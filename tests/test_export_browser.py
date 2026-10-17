"""The PDF and SVG of a changed view, pressed in headless Chromium (D-033, amended 2026-10-01).

Rytis found twice that nothing a reader changes on a page reached its PDF and SVG buttons. This
builds a module page and the standalone pages of the three kinds of page control from the
test demos (a linked view, a volcano with live calls, bars with filters), runs
``tests/browser.mjs`` under node against them, presses each button after moving a control,
and reads the files: page size in mm, embedded fonts, and the change itself. Beside them it
builds a thinned volcano, a linked view with long levels and a band map with long run names,
and measures the pages at 375, 1,280 and 1,440 px, the figures opened alone at 1,280, 1,512
and 2,560 px, and pages whose window is resized after they open. Skipped where no Chromium is
found; set ``CTX_CHROMIUM`` to one.
"""

import base64
import dataclasses
import io
import json
import os
import re
import shutil
import signal
import subprocess
from datetime import date
from pathlib import Path
from xml.etree import ElementTree

import polars as pl
import pytest
from fontTools.ttLib import TTFont
from pypdf import PdfReader
from synth import counts, demo_linked, demo_volcano
from test_viz import _painted, _rgb, filtered_bars

from ctx.io.naming import stem
from ctx.report.builder import (
    FigureText,
    Group,
    ModulePage,
    ParameterCard,
    Reproduction,
    SelfCheck,
    write,
)
from ctx.report.figures import finite
from ctx.report.text import word as text_word
from ctx.viz import scatter
from ctx.viz.aggregates import Aggregate, Lead
from ctx.viz.save import figure_block, write_figure
from ctx.viz.style import (
    ASPECT,
    BAND_COLOURS,
    DATA_COLOURS,
    LINE_WIDTHS,
    MARKER_SIZE,
    PAGE_PLOT_PX,
    PHONE_PLOT_PX,
    PLOTLY_MARKER_PX,
    SIZES,
    THEMES,
    WIDTHS_MM,
)

ROOT = Path(__file__).resolve().parents[1]
JSPDF_START_WIDTH = 0.200025
CHROMIUM = os.environ.get("CTX_CHROMIUM") or next(
    filter(None, map(shutil.which, ("chromium", "chromium-browser", "google-chrome"))), None
)
pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
SVG = "{http://www.w3.org/2000/svg}"
MM = 72 / 25.4
# Rytis, 2026-10-02: the page's subset of Liberation Sans is renamed; its PostScript names.
RENAMED = ("CtxSans-Regular", "CtxSans-Bold")
# The two paints of the linked view, indices into the data colours: orange and vermillion.
PAINTS = (0, 5)
READING = ("Each mark is one row of the demo.",)
# Rytis, 2026-10-01: opened alone, a figure is drawn 2x bigger; stated apart from style.
TWICE = 2
# plotly draws a legend's marker at most 16 px across (style.PLOTLY_LEGEND_KEY_PX).
LEGEND_MARKER_CAP_PX = 16
# Rytis, 2026-10-02: opened alone, a figure keeps twice the page's type but no wider than its
# window; laptop windows of 1,280 and 1,512 px, and one wide enough for the full width.
ALONE_WINDOWS = (1280, 1512, 2560)
# Each figure's plots among the module page's, in page order: the linked view's scatter and
# its parallel sets, the volcano, the bars.
ON_PAGE = {"volcano": (2,), "linked": (0, 1)}
# Run names in capitals, 19 characters at the longest: the review's band map whose title ran
# into its column labels (2026-10-01).
RUNS = ("RAP1_UNINDUCED_REP1", "RAP1_UNINDUCED_REP2", "RAP1_IAA_30M_REP1", "WT_REP1", "WT_REP2")
# The windows the band map is opened in and then resized to, in order (review, 2026-10-02):
# alone from a phone to a laptop, a wide screen and back, through its narrow layout and out;
# the module page from a phone to a laptop.
BANDS_ALONE_WINDOWS = (375, 1280, 2560, 375)
BANDS_PAGE_WINDOWS = (375, 1440)
# Rytis, 2026-10-07: below the phone width, Open alone draws a figure as the page does, and 2x
# from the floor up. A phone, then windows whose box falls either side of the 634 px floor.
FLOOR_WINDOWS = (375, 700, 768)

LABELS = """Array.from(document.querySelectorAll(%s)).map(function (a) {
  return [a.textContent, a.getAttribute("href")]; })"""
PAINT = """(function () {
  var fig = document.querySelector("#figure-1"), links = %(links)s;
  var before = links();
  var div = fig.querySelector(".fig-plot .js-plotly-plot");
  var sets = fig.querySelector(".fig-below .js-plotly-plot");
  var colour = function (k) { fig.querySelector("input[name=colour][value='" + k + "']").click(); };
  colour(%(first)d);
  div.emit("plotly_selected", {points: [0, 1, 2, 3, 4, 5].map(function (i) {
    return {curveNumber: 0, pointNumber: i, pointIndex: i}; })});
  colour(%(second)d);
  sets.emit("plotly_click", {points: [20, 21, 22].map(function (i) {
    return {curveNumber: 0, pointNumber: i}; })});
  return {before: before, after: links(), summary: fig.querySelector(".brush-summary").textContent,
          markers: div.data[0].marker.color, codes: Array.from(sets.data[0].line.color)};
})()"""
CALLS = """(function () {
  var fig = document.querySelector("%(figure)s"), form = fig.querySelector("form.calls");
  var links = %(links)s, before = links();
  form.elements.fc.value = "2";
  form.elements.fc.dispatchEvent(new Event("change", {bubbles: true}));
  form.elements.cutoff.value = "0.01";
  form.elements.cutoff.dispatchEvent(new Event("change", {bubbles: true}));
  var div = fig.querySelector(".js-plotly-plot"), cfg = JSON.parse(
    document.getElementById(form.getAttribute("data-controls")).textContent);
  var lines = cfg.cutoff_lines.map(function (i) { return div.data[i]; });
  return {before: before, after: links(), summary: fig.querySelector(".calls-summary").textContent,
          vertical: [lines[0].x[0], lines[1].x[0]], horizontal: lines[2].y[0],
          xtitle: div.layout.xaxis.title.text, up: div.data[cfg.traces[cfg.levels.up]].x.length,
          down: div.data[cfg.traces[cfg.levels.down]].x.length,
          roles: Array.from(fig.querySelectorAll("a[data-export]")).map(function (a) {
            return [a.getAttribute("role"), a.getAttribute("tabindex")]; })};
})()"""
# At the defaults the buttons open the figure files; the test draws the view anyway, through
# the page's own ctxExport, to compare it with them.
DRAW = """(function () {
  var fig = document.querySelector("%(figure)s");
  var cfg = JSON.parse(fig.querySelector("script[id^=ctx-controls]").textContent);
  return Promise.all(["pdf", "svg"].map(function (ext) {
    return ctxExport.draw(fig, cfg, ext).then(function (blob) {
      return new Promise(function (done) {
        var reader = new FileReader(); reader.onload = function () { done(reader.result); };
        reader.readAsDataURL(blob);
      });
    });
  }));
})()"""
# Regular and bold text through the page's own PDF path, with the fonts the page carries.
BOLD = """(function () {
  var setup = JSON.parse(document.getElementById("ctx-export").textContent);
  var svg = '<svg xmlns="http://www.w3.org/2000/svg" width="240pt" height="100pt" ' +
    'viewBox="0 0 240 100"><text x="10" y="30" style="font-family: Arial; font-size: 9px">' +
    'Regular −2 ≥ ρ</text><text x="10" y="60" style="font-family: Arial; font-size: 12px; ' +
    'font-weight: bold">Bold ρ</text></svg>';
  return ctxExport.pdf(svg, {width_pt: 240, height_pt: 100}, setup).then(function (blob) {
    return new Promise(function (done) {
      var reader = new FileReader(); reader.onload = function () { done(reader.result); };
      reader.readAsDataURL(blob);
    });
  });
})()"""
BOX = """function box(e) { var r = e.getBoundingClientRect();
  return [r.left, r.top, r.right, r.bottom]; }"""
# Every box a page plot draws, for a figure.
PLOTS = """(function () {
  %(box)s
  return Array.from(document.querySelectorAll("%(figure)s .js-plotly-plot")).map(function (d) {
    var all = function (s) { return Array.from(d.querySelectorAll(s)).map(box); };
    var hint = d.closest(".fig-plot, .fig-below").querySelector(".scroll-hint");
    var narrow = d.layout.meta && d.layout.meta.narrow;
    return {svg: box(d.querySelector("svg.main-svg")), xticks: all(".xtick text"),
            yticks: all(".ytick text"), xtitle: all(".xtitle"), dims: all(".dimlabel"),
            annotations: all(".annotation"), legendpoints: all(".legendpoints path"),
            legendlines: Array.from(d.querySelectorAll(".legendlines path")).map(function (e) {
              return [box(e), parseFloat(getComputedStyle(e).strokeWidth)]; }),
            width: d.getBoundingClientRect().width, hint: !hint.hidden,
            box: box(d.closest(".fig-plot, .fig-below")),
            edges: Array.from(d.data[0].y || []),
            narrow: narrow && narrow.restyle ? narrow.restyle[0][0].y[0] : null};
  });
})()"""
LAYOUT = """(function () {
  %(box)s
  var fig = document.querySelector("#figure-1");
  return {plot: box(fig.querySelector(".fig-plot")), reading: box(fig.querySelector(".reading")),
          below: box(fig.querySelector(".fig-below")),
          sets: fig.querySelector(".fig-below .js-plotly-plot").getBoundingClientRect().width};
})()"""
# Every text of an SVG file opened alone, in the file's own units, and the file's size.
TEXTS = """(function () {
  var svg = document.documentElement, view = svg.viewBox.baseVal;
  var origin = svg.getBoundingClientRect(), k = view.width / origin.width;
  return {page: [view.width, view.height], texts: Array.from(svg.querySelectorAll("text"))
    .filter(function (t) { return t.textContent.trim(); }).map(function (t) {
      var r = t.getBoundingClientRect();
      return [t.textContent, (r.left - origin.left) * k, (r.right - origin.left) * k,
              (r.top - origin.top) * k, (r.bottom - origin.top) * k]; })};
})()"""
FILTER = """(function () {
  var fig = document.querySelector("#figure-3"), select = fig.querySelector("select[data-column]");
  var links = %(links)s, before = links(), lead = fig.querySelector(".fig-lead");
  var was = lead.textContent;
  select.value = "LDN"; select.dispatchEvent(new Event("change", {bubbles: true}));
  var div = fig.querySelector(".js-plotly-plot");
  return {before: before, after: links(), summary: fig.querySelector(".filter-summary").textContent,
          shown: div.layout.xaxis.categoryarray, lead: [was, lead.textContent],
          spans: lead.querySelectorAll("[data-fill]").length};
})()"""
# Every plot's drawn and div width, then again once a 15 px margin narrows every card with no
# window resize, as a scrollbar appearing after the first draw does: polled until each plot is
# drawn at its div's width, or 3 s.
NARROWED = """new Promise(function (done) {
  function widths() {
    return Array.from(document.querySelectorAll(".js-plotly-plot")).map(function (d) {
      return {drawn: d.querySelector("svg.main-svg").getBoundingClientRect().width,
              div: d.getBoundingClientRect().width};
    });
  }
  var before = widths(), start = Date.now();
  document.querySelectorAll(".fig").forEach(function (f) { f.style.marginRight = "15px"; });
  (function poll() {
    var after = widths();
    var fits = after.every(function (p) { return Math.abs(p.drawn - p.div) <= 1; });
    if ((fits && Date.now() - start > 300) || Date.now() - start > 3000) {
      done({before: before, after: after});
    } else setTimeout(poll, 100);
  })();
})"""
# Every value of every plot's customdata that is null, which plotly prints as 0 under a number
# format, and every one that is "", which it prints as nothing.
HOVERED = """(function () {
  var nulls = 0, blanks = 0;
  document.querySelectorAll(".js-plotly-plot").forEach(function (d) {
    d.data.forEach(function (t) {
      (t.customdata || []).forEach(function (row) {
        (Array.isArray(row) ? row : [row]).forEach(function (v) {
          if (v === null) nulls += 1; else if (v === "") blanks += 1;
        });
      });
    });
  });
  return {nulls: nulls, blanks: blanks};
})()"""
# The demo bars' lead, whose numbers follow the filter as M1's and M3's do (D1).
BARS_LEAD = Lead(
    "{runs} runs shown; median {median} million reads, {short} below thirty.",
    (Aggregate("runs", "count"), Aggregate("median", "median", "reads"),
     Aggregate("short", "below", "reads", value=30.0)),
    note="The caption counts every run.",
)  # fmt: skip
RESET = """new Promise(function (done) {
  var fig = document.querySelector("%(figure)s");
  %(undo)s
  setTimeout(function () { done(%(links)s); }, 50);
})"""
UNDO = {
    "linked": 'fig.querySelector("[data-clear]").click();',
    "volcano": 'fig.querySelector("form.calls").reset();',
    "bars": 'var s = fig.querySelector("select[data-column]"); s.value = "";'
    's.dispatchEvent(new Event("change", {bubbles: true}));',
}
STACK = """(function () {
  var fig = document.querySelector("#figure-1"), box = function (s) {
    var r = fig.querySelector(s).getBoundingClientRect();
    return [r.left, r.top, r.right, r.bottom];
  };
  return {plot: box(".fig-plot"), reading: box(".reading"), below: box(".fig-below"),
          page: document.documentElement.clientWidth, scroll: document.documentElement.scrollWidth};
})()"""
SIZES_JS = """(function () {
  var div = document.querySelector(".fig-plot .js-plotly-plot");
  var tick = div.querySelector(".xtick text");
  return {width: div.getBoundingClientRect().width,
          tick: parseFloat(getComputedStyle(tick).fontSize),
          marker: div._fullData[0].marker.size};
})()"""


# Every plot of a page: the width plotly drew it, its div's width, its box's inner width,
# whether its scroll hint shows, the type size of its first tick label or dimension title, the
# boxes of its drawing, its row labels, the parallel sets' dimension titles, their levels and
# their bars; and the page's width and scroll width.
FIT = """(function () {
  %s
  return {page: document.documentElement.clientWidth, scroll: document.documentElement.scrollWidth,
    plots: Array.from(document.querySelectorAll(".js-plotly-plot")).map(function (d) {
      var holder = d.closest(".fig-plot, .fig-below");
      var text = d.querySelector(".xtick text, .ytick text, .dimlabel");
      var all = function (s) { return Array.from(d.querySelectorAll(s)).map(box); };
      var keys = Array.from(d.querySelectorAll(".legend .traces")).map(function (t) {
        var marks = Array.from(t.querySelectorAll(".legendpoints path, .legendlines path"));
        return [Math.max.apply(null, marks.map(function (m) { return box(m)[2]; })),
                box(t.querySelector(".legendtext"))[0]];
      });
      return {drawn: d.querySelector("svg.main-svg").getBoundingClientRect().width,
              div: d.getBoundingClientRect().width,
              box: holder.clientWidth, hint: !holder.querySelector(".scroll-hint").hidden,
              text: text && parseFloat(getComputedStyle(text).fontSize),
              svg: box(d.querySelector("svg.main-svg")), yticks: all(".ytick text"),
              titles: all(".dimlabel"), levels: all(".catlabel"), bars: all(".catrect"),
              keys: keys};
    })};
})()""".replace("%s", BOX, 1)
# FIT once plotly has drawn every plot at its div's width, or after 5 s: plotly redraws a plot
# on a window's resize after its own short delay, and a page's narrow layout swaps in as the
# div's width crosses its threshold.
SETTLED = """new Promise(function (done) {
  var start = Date.now();
  (function poll() {
    var fit = %s;
    var settled = fit.plots.every(function (p) { return Math.abs(p.drawn - p.div) <= 0.5; });
    if (settled || Date.now() - start > 5000) done(fit); else setTimeout(poll, 100);
  })();
})""".replace("%s", FIT, 1)


# SETTLED, then the theme toggled and the page measured again.
TOGGLED = (
    SETTLED
    + """.then(function (first) {
  document.getElementById("ctx-theme-toggle").click();
  return new Promise(function (done) { setTimeout(function () {
    done({first: first, toggled: %s}); }, 500); });
})""".replace("%s", FIT, 1)
)


def _links(figure: str) -> str:
    return "function () { return " + LABELS % json.dumps(f"{figure} a[data-export]") + "; }"


def _module(blocks) -> ModulePage:
    return ModulePage(
        title="M2 Demo: synthetic", purpose="Purpose.", input="Input.", method="Method.",
        parameters=(ParameterCard("--x", "A card.", "1", "Moves it."),),
        results=(Group("Figures", tuple(blocks)),), interpretation="Interpretation.",
        questions=(SelfCheck("A question.", "An answer."),), pitfalls=("A pitfall.",),
        reproduce=Reproduction("ctx demo", 1.0, 1, 1, "laptop", "none", {"python": "3.12"}),
    )  # fmt: skip


def _long_levels():
    """The demo's linked view with outer levels longer than the scatter's margins."""
    data = demo_linked()
    frame = data.frame.with_columns(
        tissue=pl.col("tissue").replace({"PB": "peripheral blood", "BM": "bone marrow aspirate"}),
        qc_status=pl.col("qc_status").replace({"pass": "passes every band",
                                                "flag": "flagged for review"}),
    )  # fmt: skip
    colours = {"passes every band": BAND_COLOURS["good"], "flagged for review": BAND_COLOURS["bad"]}
    return dataclasses.replace(data, frame=frame, colours=colours)


def _band_map():
    """Every banded and FastQC row of the registry over five runs named in long capitals."""
    from ctx.qc import figures  # noqa: PLC0415
    from ctx.qc.registry import load  # noqa: PLC0415

    registry = load(ROOT / "conf" / "qc_metrics.toml")
    metrics = [m for m in registry.metrics.values() if m.banded or m.kind == "status"]
    values = pl.DataFrame({"run_accession": list(RUNS),
                           **{m.column: ["pass" if m.kind == "status" else m.good] * len(RUNS)
                              for m in metrics}})  # fmt: skip
    return figures.band_map(values, registry, "run_accession")


def _write(out: Path, page: str, demos: dict) -> tuple[Path, dict]:
    """Each figure's files and standalone page, and a module page holding them all."""
    blocks = []
    for name, data in demos.items():
        figure_stem = stem("qc", f"{name}_demo", 1, day=date(2026, 10, 1))
        text = FigureText(title=name) if name == "volcano" else FigureText(name, reading=READING)
        if name == "bars":
            text = FigureText(name, reading=READING, follows=BARS_LEAD)
        write_figure(data, out, figure_stem, title=name, data_releases={}, text=text)
        blocks.append(figure_block(data, figure_stem, text))
    overview = stem("qc", f"{page}_overview", 1, day=date(2026, 10, 1))
    written = write(_module(blocks), out / f"{overview}.html")
    return written, {b.stem.split("_")[2]: b for b in blocks}


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """Each kind's figure files and standalone page, and a module page holding all three; a
    page each for a thinned volcano, long linked levels and a band map."""
    out = tmp_path_factory.mktemp("pages")
    demos = {"linked": demo_linked(), "volcano": demo_volcano(counts(), width="double"),
             "bars": filtered_bars()}  # fmt: skip
    page, blocks = _write(out, "demo", demos)
    pages = {"demo": page}
    with pytest.MonkeyPatch.context() as patch:
        # The page draws every 3rd of the 704 points no setting can call.
        patch.setattr(scatter, "THIN_MAX", 300)
        thin = demo_volcano(counts(), width="double")
        pages["thin"], thinned = _write(out / "thin", "thin", {"volcano": thin})
    blocks = {**blocks, "thin": thinned["volcano"]}
    pages["long"], _ = _write(out / "long", "long", {"linked": _long_levels()})
    pages["bands"], _ = _write(out / "bands", "bands", {"bandmap": _band_map()})
    pages["untested"], _ = _write(out / "untested", "untested", {"volcano": _untested()})
    return out, page, {**demos, "thin": thin}, blocks, pages


def _untested():
    """The demo volcano with every seventh hover value missing: a gene with no identifier."""
    data = demo_volcano(counts(), width="double")
    hover = data.frame.with_row_index().select(
        pl.exclude("index", "gene_id"),
        gene_id=pl.when(pl.col("index") % 7 > 0).then(pl.col("gene_id")))  # fmt: skip
    return dataclasses.replace(data, frame=hover.select(data.frame.columns))


@pytest.fixture(scope="module")
def ran(built, tmp_path_factory):
    out, page, _, blocks, pages = built
    downloads = tmp_path_factory.mktemp("downloads")
    url = page.as_uri()
    alone = (out / f"{blocks['linked'].stem}.html").as_uri()
    tasks = []
    paints = {"first": PAINTS[0], "second": PAINTS[1]}
    kinds = (("linked", 1, PAINT, paints), ("volcano", 2, CALLS, {}), ("bars", 3, FILTER, {}))
    for name, number, act, extra in kinds:
        figure = f"#figure-{number}"
        tasks.append({
            "name": name, "url": url, "figure": figure,
            "act": act % {"links": _links(figure), "figure": figure, **extra},
            "reset": RESET % {"figure": figure, "undo": UNDO[name], "links": _links(figure) + "()"},
        })  # fmt: skip
    one = "#figure-1"
    tasks += [
        {
            "name": "alone_paint",
            "url": alone,
            "figure": one,
            "act": PAINT % {"links": _links(one), "first": PAINTS[0], "second": PAINTS[1]},
        },
        {"name": "thin", "url": pages["thin"].as_uri(), "figure": one,
         "act": CALLS % {"links": _links(one), "figure": one}},
        {"name": "volcano_phone", "url": url, "width": 375, "height": 812, "figure": "#figure-2",
         "act": CALLS % {"links": _links("#figure-2"), "figure": "#figure-2"}},
        {"name": "long", "url": pages["long"].as_uri(), "figure": one,
         "act": PAINT % {"links": _links(one), **paints}},
        {"name": "bold", "url": url, "act": BOLD},
        {"name": "defaults", "url": url, "act": "Promise.all([" + ", ".join(
            DRAW % {"figure": f"#figure-{n}"} for n in (2, 3)) + "])"},
        {"name": "phone", "url": url, "width": 375, "height": 812,
         "act": "({stack: " + STACK + ", linked: " + PLOTS % {"box": BOX, "figure": one}
                + ", volcano: " + PLOTS % {"box": BOX, "figure": "#figure-2"} + "})"},
        {"name": "page", "url": url, "act": SIZES_JS},
        {"name": "alone", "url": alone, "width": 2560, "act": SIZES_JS},
        {"name": "page_volcano", "url": url, "act": PLOTS % {"box": BOX, "figure": "#figure-2"}},
        {"name": "alone_volcano", "url": (out / f"{blocks['volcano'].stem}.html").as_uri(),
         "act": PLOTS % {"box": BOX, "figure": one}},
    ]  # fmt: skip
    for width in (1280, 1440):
        tasks.append({"name": f"desktop_{width}", "url": url, "width": width,
                      "act": LAYOUT % {"box": BOX}})  # fmt: skip
    standalone = {kind: (out / f"{blocks[kind].stem}.html").as_uri() for kind in ON_PAGE}
    tasks.append({"name": "fit_page", "url": url, "width": 1440, "act": FIT})
    for kind, where in standalone.items():
        for width in ALONE_WINDOWS:
            tasks.append({"name": f"fit_{kind}_{width}", "url": where, "width": width, "act": FIT})
    tasks.append({"name": "reflow", "url": standalone["volcano"], "width": 2560, "act": FIT,
                  "resize": [{"width": 1280, "act": SETTLED}]})  # fmt: skip
    bands_stem = stem("qc", "bandmap_demo", 1, day=date(2026, 10, 1))
    bands_alone = pages["bands"].parent / f"{bands_stem}.html"
    for name, where, width in (("bands_1280", pages["bands"], 1280),
                               ("bands_375", pages["bands"], 375),
                               ("bands_alone", bands_alone, 1440),
                               ("bands_alone_768", bands_alone, 768),
                               ("bands_alone_375", bands_alone, 375)):  # fmt: skip
        tasks.append({"name": name, "url": where.as_uri(), "width": width, "height": 812,
                      "act": PLOTS % {"box": BOX, "figure": one}})  # fmt: skip
    for name, where, widths in (("bands_alone_reflow", bands_alone, BANDS_ALONE_WINDOWS),
                                ("bands_reflow", pages["bands"], BANDS_PAGE_WINDOWS)):  # fmt: skip
        steps = [{"width": width, "height": 812, "act": SETTLED} for width in widths]
        tasks.append({**steps[0], "name": name, "url": where.as_uri(), "resize": steps[1:]})
    tasks += _fit_tasks(page, pages, standalone["volcano"], bands_alone.as_uri())
    results = _browse(tasks, downloads)
    for result in results.values():
        result["paths"] = {ext: downloads / f["guid"] for ext, f in result["files"].items()}
    for ext, encoded in zip(("pdf", "svg"), results["defaults"]["acted"][0], strict=True):
        results["defaults"]["paths"][f"volcano_{ext}"] = _decoded(encoded, downloads, f"v.{ext}")
    for ext, encoded in zip(("pdf", "svg"), results["defaults"]["acted"][1], strict=True):
        results["defaults"]["paths"][f"bars_{ext}"] = _decoded(encoded, downloads, f"b.{ext}")
    results["bold"]["paths"]["pdf"] = _decoded(results["bold"]["acted"], downloads, "bold.pdf")
    return results, blocks


def _fit_tasks(page: Path, pages: dict, volcano: str, bands: str) -> list[dict]:
    """Every card narrowed with no window resize, on the demo and band map pages; the volcano
    and the band map opened alone either side of the floor; the volcano whose hover misses
    values, drawn by its live calls."""
    tasks = [{"name": f"narrowed_{name}", "url": where.as_uri(), "width": 1280, "act": NARROWED}
             for name, where in (("demo", page), ("bands", pages["bands"]))]  # fmt: skip
    tasks += [{"name": f"floor_{kind}_{width}", "url": where, "width": width, "height": 812,
               "act": TOGGLED}
              for kind, where in (("volcano", volcano), ("bandmap", bands))
              for width in FLOOR_WINDOWS]  # fmt: skip
    one = "#figure-1"
    tasks.append({"name": "alone_volcano_phone", "url": volcano, "width": 375, "height": 812,
                  "figure": one, "act": CALLS % {"links": _links(one), "figure": one}})  # fmt: skip
    return [*tasks, {"name": "untested", "url": pages["untested"].as_uri(), "act": HOVERED}]


def _browse(tasks: list[dict], downloads: Path, chromium: str | None = None,
            timeout: float = 900) -> dict[str, dict]:  # fmt: skip
    """Run tasks in headless Chromium through tests/browser.mjs; results by task name.

    The driver traces each step to stderr, so one past ``timeout`` names the step it stopped
    in. It runs in a session of its own, and past ``timeout`` the whole session is killed:
    killing the driver alone left its Chromium running, four of them for days (K.13).
    """
    plan = {"chromium": chromium or CHROMIUM, "downloads": str(downloads), "tasks": tasks,
            "trace": True}  # fmt: skip
    node = shutil.which("node")
    assert node, "nodejs is pinned in the ctx-dev environment; run tests with pixi"
    with subprocess.Popen([node, str(ROOT / "tests" / "browser.mjs")], stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                          start_new_session=True) as proc:  # fmt: skip
        try:
            stdout, stderr = proc.communicate(json.dumps(plan), timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            _, stderr = proc.communicate()
            raise AssertionError(f"browser driver past {timeout:g} s:\n{stderr[-2000:]}") from None
    # The driver's own error, which a bare exit status hides.
    assert proc.returncode == 0, stderr[-4000:]
    return {r["name"]: r for r in json.loads(stdout)}


def _decoded(data_url: str, folder: Path, name: str) -> Path:
    path = folder / name
    path.write_bytes(base64.b64decode(data_url.split(",", 1)[1]))
    return path


@pytest.fixture(scope="module")
def texts(ran, tmp_path_factory):
    """Every text of every SVG the buttons gave, measured in Chromium in the file's own pt."""
    results, _ = ran
    folder = tmp_path_factory.mktemp("svgs")
    tasks = []
    for name, result in results.items():
        for key, path in result["paths"].items():
            if key.endswith("svg"):
                copy = folder / f"{name}_{key}.svg"
                copy.write_bytes(path.read_bytes())
                tasks.append({"name": f"{name}_{key}", "url": copy.as_uri(), "act": TEXTS})
    return {name: r["acted"] for name, r in _browse(tasks, folder).items()}


def _svg(path: Path) -> ElementTree.Element:
    return ElementTree.parse(path).getroot()


def _hex(value: str) -> str:
    """A colour plotly writes, rgb(r, g, b), black or #RRGGBB, as #RRGGBB."""
    if found := re.fullmatch(r"rgb\((\d+), ?(\d+), ?(\d+)\)", value.strip()):
        return "#" + "".join(f"{int(v):02X}" for v in found.groups())
    return {"black": "#000000"}.get(value.strip(), value.strip().upper())


def _paint(element: ElementTree.Element, what: str) -> str | None:
    """An element's fill or stroke, from its style or its attribute."""
    found = re.search(rf"(?:^|;)\s*{what}:\s*([^;]+)", element.get("style", ""))
    value = found.group(1) if found else element.get(what)
    return None if value is None else _hex(value)


def _fills(root: ElementTree.Element, cls: str) -> list[str]:
    """The fill of every element whose class starts with ``cls``, as #RRGGBB."""
    found = (_paint(e, "fill") for e in root.iter() if e.get("class", "").startswith(cls))
    return [f for f in found if f]


# plotly's legend scrollbar, drawn 0 by 0 when the legend fits; no reader sees it.
UNSEEN = ("scrollbar",)


def _stray_colours(root: ElementTree.Element) -> set[str]:
    """Every fill and stroke the SVG paints that no style token names."""
    theme = THEMES["white"]
    allowed = {*DATA_COLOURS, theme.ink, theme.grid, theme.reference, "NONE"}
    painted = {_paint(e, what) for e in root.iter() if e.get("class") not in UNSEEN
               for what in ("fill", "stroke")}  # fmt: skip
    return painted - allowed - {None}


def _pdf(path: Path) -> tuple[tuple[float, float], set[str], str]:
    """Page size in mm, the fonts embedded with their file, and the text."""
    page = PdfReader(path).pages[0]
    size = (float(page.mediabox.width) / MM, float(page.mediabox.height) / MM)
    fonts = set()
    for entry in page["/Resources"]["/Font"].values():
        font = entry.get_object()
        for descendant in font["/DescendantFonts"]:
            descriptor = descendant.get_object()["/FontDescriptor"].get_object()
            assert "/FontFile2" in descriptor, "a font named but not embedded"
        fonts.add(str(font["/BaseFont"]).lstrip("/"))
    return size, fonts, " ".join(page.extract_text().split())


def _programs(path: Path) -> dict[str, dict[int, set[str]]]:
    """Each embedded font program's name records by name ID, by the PDF's name for the font."""
    page, found = PdfReader(path).pages[0], {}
    for entry in page["/Resources"]["/Font"].values():
        font = entry.get_object()
        for descendant in font["/DescendantFonts"]:
            descriptor = descendant.get_object()["/FontDescriptor"].get_object()
            data = descriptor["/FontFile2"].get_object().get_data()
            names: dict[int, set[str]] = {}
            for record in TTFont(io.BytesIO(data))["name"].names:
                names.setdefault(record.nameID, set()).add(record.toUnicode())
            found[str(font["/BaseFont"]).lstrip("/")] = names
    return found


def _near(a: tuple[float, ...], b: tuple[float, ...]) -> bool:
    """Equal as the export writes a colour, to four decimals per channel. jsPDF wrote two, up
    to 0.0047 off, which a tolerance of 0.005 let through (review, 2026-10-01)."""
    return all(abs(x - y) <= 0.00005 + 1e-9 for x, y in zip(a, b, strict=True))


def _expected_mm(data) -> tuple[float, float]:
    width, height = data.inches()
    return width * 25.4, height * 25.4


@pytest.mark.parametrize("name", ["linked", "volcano", "bars"])
def test_buttons_give_the_files_at_the_defaults_and_say_when_they_draw(ran, name):
    """At the defaults PDF and SVG link the figure files; moved, they name the view and link
    no file, so no middle click or "Save link as" reaches the default file under the view's
    label; undone, the files again. Mutations: the switch never called from a control's
    wiring; the href kept while moved."""
    results, blocks = ran
    result, figure_stem = results[name], blocks[name].stem
    files = [
        f"{figure_stem}_{theme}.{ext}" for ext in ("pdf", "svg") for theme in ("white", "black")
    ]
    reference = [
        [ext.upper(), f] for ext, f in zip(("pdf",) * 2 + ("svg",) * 2, files, strict=True)
    ]
    assert result["acted"]["before"] == reference
    assert [label for label, _ in result["acted"]["after"]] == ["PDF of this view"] * 2 + [
        "SVG of this view"] * 2  # fmt: skip
    assert [href for _, href in result["acted"]["after"]] == [None] * 4
    assert result["reset"] == reference
    assert {f["name"] for f in result["files"].values()} == {f"{figure_stem}_view.pdf",
                                                             f"{figure_stem}_view.svg"}  # fmt: skip


@pytest.mark.parametrize("name", ["linked", "volcano", "bars"])
def test_files_take_the_figure_files_size_and_embed_ctx_sans(ran, built, name):
    """85 or 170 mm wide as the figure's own files, Ctx Sans embedded and no other font,
    9 pt text and 10 pt titles, and the SVG sized in pt as matplotlib's is; the caption ends
    with the live line and the note that the file is the page's drawing (the plan's question
    1, Rytis 2026-10-01). Mutations: the face left to svg2pdf (Helvetica), the SVG left in
    plotly's px, the note left out of the caption."""
    results, _ = ran
    data = built[2][name]
    width, height = _expected_mm(data)
    assert width == WIDTHS_MM[data.width]
    size, fonts, text = _pdf(results[name]["paths"]["pdf"])
    assert size == pytest.approx((width, height), abs=0.01)
    assert fonts == {"CtxSans-Regular"}
    colours, sizes, widths, page_fill = _painted(results[name]["paths"]["pdf"])
    assert sizes <= {SIZES["tick"], SIZES["axis_label"]} and not page_fill
    # jsPDF opens every page at its own 0.200025 pt before a stroke sets another.
    assert widths - {JSPDF_START_WIDTH} <= set(LINE_WIDTHS.values())
    theme = THEMES["white"]
    allowed = [_rgb(c) for c in (*DATA_COLOURS, theme.ink, theme.grid, theme.reference)]
    assert all(any(_near(c, a) for a in allowed) for c in colours), colours
    root = _svg(results[name]["paths"]["svg"])
    assert _stray_colours(root) == set()
    assert float(root.get("width").removesuffix("pt")) == pytest.approx(width * MM)
    assert float(root.get("height").removesuffix("pt")) == pytest.approx(height * MM)
    assert results[name]["acted"]["summary"] in text
    note = text_word("ui", "export.note")
    assert text.endswith(note) and text.count(note) == 1


def test_bold_text_embeds_the_renamed_faces_from_the_subset(ran):
    """The page's fonts, cut to its characters and renamed (Rytis 2026-10-02), still draw:
    regular and bold text through the page's PDF path embed CtxSans-Regular and CtxSans-Bold,
    each under its own name; each embedded program's own names are Ctx Sans's, with no
    Liberation in its family, unique, full or PostScript name; and the minus sign, the
    greater-than-or-equal sign and rho come back as text. Mutations: the PostScript names not
    set, so both faces were embedded as one name; the subset not renamed, so the PDF named
    CtxSans-Regular a program whose own names were Liberation Sans's."""
    path = ran[0]["bold"]["paths"]["pdf"]
    _, fonts, text = _pdf(path)
    assert fonts == set(RENAMED)
    for name, records in _programs(path).items():
        assert (records[1], records[6]) == ({"Ctx Sans"}, {name})
        assert not any("Liberation" in v for i in (1, 3, 4, 6) for v in records[i]), records
        assert "Copyright (c) 2012 Red Hat, Inc." in next(iter(records[0]))
    assert "Regular −2 ≥ ρ" in text and "Bold ρ" in text


def test_painted_view_draws_both_views_in_the_paint(ran, built):
    """The painted samples in their paints in both views, in the SVG and the PDF.
    Mutation: the parcats restyled on the scatter's plot."""
    results, _ = ran
    acted = results["linked"]["acted"]
    paints = [DATA_COLOURS[k] for k in PAINTS]
    assert acted["markers"][:6] == [paints[0]] * 6 and acted["markers"][20:23] == [paints[1]] * 3
    assert (
        acted["codes"][:6] == [PAINTS[0] + 1] * 6 and acted["codes"][20:23] == [PAINTS[1] + 1] * 3
    )
    root = _svg(results["linked"]["paths"]["svg"])
    points = _fills(root, "point")
    assert [points.count(p) for p in paints] == [acted["markers"].count(p) for p in paints]
    bands = set(_fills(root, "path"))
    assert set(paints) <= bands, "the parallel sets carry no paint"
    colours, *_ = _painted(results["linked"]["paths"]["pdf"])
    assert all(any(_near(c, _rgb(p)) for c in colours) for p in paints)
    height = WIDTHS_MM["single"] * ASPECT["linked"]
    assert _pdf(results["linked"]["paths"]["pdf"])[0][1] == pytest.approx(height, abs=0.01)


def _axis_map(root: ElementTree.Element, axis: str) -> tuple[float, float]:
    """(px per unit, px of zero) of an axis, from its numbered tick labels."""
    ticks = []
    for group in root.iter(f"{SVG}g"):
        if group.get("class") != f"{axis}tick":
            continue
        label = group.find(f"{SVG}text")
        moved = re.search(r"translate\(([-\d.]+),([-\d.]+)\)", label.get("transform", ""))
        value = label.text.replace("−", "-") if label is not None and label.text else ""
        if moved and re.fullmatch(r"-?\d+(\.\d+)?", value):
            ticks.append((float(value), float(moved.group(1 if axis == "x" else 2))))
    (v0, p0), (v1, p1) = ticks[0], ticks[-1]
    scale = (p1 - p0) / (v1 - v0)
    return scale, p0 - v0 * scale


def test_moved_cutoffs_draw_the_calls_and_lines_at_the_new_values(ran):
    """Up and down at the new rule, and the cutoff lines at ±2 and at the new p line, read
    off the SVG's own axes. Mutation: the export drawn from the figure as first loaded."""
    results, _ = ran
    acted = results["volcano"]["acted"]
    assert acted["vertical"] == [2, -2]
    root = _svg(results["volcano"]["paths"]["svg"])
    fills = _fills(root, "point")
    up, down = "#D55E00", "#0072B2"
    assert (fills.count(up), fills.count(down)) == (acted["up"], acted["down"])
    assert f"Up in case: {acted['up']}." in acted["summary"]
    lines = [re.findall(r"[-\d.]+", p.get("d")) for p in root.iter(f"{SVG}path")
             if p.get("class") == "js-line" and "L" in p.get("d", "")]  # fmt: skip
    vertical = sorted(float(d[0]) for d in lines if d[0] == d[2])
    (horizontal,) = [float(d[1]) for d in lines if d[1] == d[3]]
    x_scale, _ = _axis_map(root, "x")
    assert (vertical[1] - vertical[0]) / 2 / x_scale == pytest.approx(2.0, abs=0.01)
    # Both in the plot's own frame: the vertical lines run up from y = 0.
    bottom = max(float(d[1]) for d in lines if d[0] == d[2])
    y_scale, _ = _axis_map(root, "y")
    assert (horizontal - bottom) / y_scale == pytest.approx(acted["horizontal"], abs=0.02)


def _xtitle(path: Path) -> list[list[str]]:
    """Each x-axis title of an SVG file, as its lines: plotly writes a line per tspan."""
    titles = []
    for text in _svg(path).iter(f"{SVG}text"):
        if "xtitle" in text.get("class", "").split():
            lines = [s.text or "" for s in text.iter(f"{SVG}tspan") if s.get("class") == "line"]
            titles.append(lines or ["".join(text.itertext())])
    return titles


def test_a_view_moved_on_a_phone_keeps_the_files_x_title(ran):
    """The page wraps an x-axis title wider than a phone's plot, but a view's files are the
    figure files' width: the volcano moved at 375 px exports the x title the page draws at
    1,280 px, one line as built. Mutations: the export drawing the title as the page last
    wrapped it, which broke "Fold change (log2 fold change)" after "fold" in the 170 mm files
    (review, 2026-10-02); the page never keeping the title as built."""
    phone, desktop = ran[0]["volcano_phone"], ran[0]["volcano"]
    assert "<br>" in phone["acted"]["xtitle"] and "<br>" not in desktop["acted"]["xtitle"]
    (on_phone,), (on_desktop,) = (_xtitle(r["paths"]["svg"]) for r in (phone, desktop))
    assert on_phone == on_desktop == [desktop["acted"]["xtitle"]]


def test_filtered_view_draws_the_bars_shown(ran):
    """The bars of the chosen level only, by name on the axis. Mutation: the filter reset
    before the export reads the plot."""
    results, _ = ran
    acted = results["bars"]["acted"]
    root = _svg(results["bars"]["paths"]["svg"])
    names = [t.text for g in root.iter(f"{SVG}g") if g.get("class") == "xtick"
             for t in g.iter(f"{SVG}text")]  # fmt: skip
    assert names == acted["shown"] and 0 < len(names) < 40
    bars = [p for g in root.iter(f"{SVG}g") if g.get("class") == "point" for p in g]
    assert len(bars) == len(names)


def test_a_filter_refills_its_lead_as_python_fills_it(ran, built):
    """D1: choosing LDN refills the lead's numbers, in their spans, with ``ctx.viz.aggregates``
    over the LDN rows, its words kept; the live line adds that the caption counts every run,
    which the export appends to the caption. Mutation: the refill removed from the page script
    (the lead kept every bar's numbers)."""
    acted = ran[0]["bars"]["acted"]
    rows = [finite(r) for r in built[3]["bars"].rows]
    every, ldn = BARS_LEAD.text(rows), BARS_LEAD.text([r for r in rows if r["cell_type"] == "LDN"])
    assert every != ldn and acted["lead"] == [every, ldn] and acted["spans"] == 3
    assert acted["summary"].endswith(" " + BARS_LEAD.note)
    text = " ".join(PdfReader(ran[0]["bars"]["paths"]["pdf"]).pages[0].extract_text().split())
    assert BARS_LEAD.note in text


def test_a_narrowed_card_redraws_every_plot_at_its_width(ran):
    """D11: a scrollbar that appears after the first draw narrows every card with no window
    resize, and plotly redraws only on one: QC band maps were laid out 582 px wide in a 567 px
    div. Every card narrowed 15 px, every plot of the demo page and the band map page is drawn
    at its div's width within 1 px. Mutation: the page script's ResizeObserver removed (each
    plot kept the width it was drawn at, 15 px wider than its div)."""
    for name in ("demo", "bands"):
        acted = ran[0][f"narrowed_{name}"]["acted"]
        pairs = zip(acted["before"], acted["after"], strict=True)
        narrowed = [b["div"] - a["div"] for b, a in pairs]
        assert any(n == pytest.approx(15, abs=1) for n in narrowed), (name, narrowed)
        assert all(abs(p["drawn"] - p["div"]) <= 1 for p in acted["after"]), (name, acted)


@pytest.mark.parametrize("kind", ["volcano", "bandmap"])
def test_alone_below_the_floor_is_drawn_as_the_page_draws_it(ran, kind):
    """D17, Rytis 2026-10-07: opened alone where the box is narrower than twice a phone's plot
    (a 375 px phone, a 700 px window), the volcano and the band map are drawn at their box's
    width with the module page's 20 px ticks and no scroll hint; from the floor up (a 768 px
    window) at twice the type, as wide as their box. The page never scrolls, every legend key
    ends before its name, and the theme toggle keeps the type each is drawn at. Mutations: the
    floor kept below the phone width (634 px in a 317 px box, ticks at 40 px, under the hint);
    the legend not drawn anew at the page's size (each key stood 15 px into its name); the
    toggle swapping in the templates at twice the size there."""
    page = ran[0]["fit_page"]["acted"]["plots"][0]["text"]
    floor, regimes = TWICE * PHONE_PLOT_PX, []
    for width in FLOOR_WINDOWS:
        acted = ran[0][f"floor_{kind}_{width}"]["acted"]
        for fit in (acted["first"], acted["toggled"]):
            (plot,) = fit["plots"]
            below = plot["box"] < floor
            assert plot["text"] == pytest.approx(page if below else TWICE * page), width
            assert plot["drawn"] == pytest.approx(plot["box"], abs=0.5), width
            assert not plot["hint"] and fit["scroll"] <= fit["page"], width
            assert plot["keys"] and all(key < name for key, name in plot["keys"]), plot["keys"]
        regimes.append(below)
    assert regimes == [True, True, False]
    assert page == pytest.approx(20)


def test_live_calls_print_a_missing_value_as_nothing(ran):
    """Package C's hover rule on the page: a volcano over filtered genes, whose untested genes
    have no adjusted p, draws its live calls on load from its rows, and no customdata value is
    null, which plotly prints as 0 under a number format. Mutation: the live calls' customdata
    of raw values."""
    acted = ran[0]["untested"]["acted"]
    assert acted["nulls"] == 0 and acted["blanks"] > 0


def test_no_request_and_no_error(ran):
    """Offline: the page asks for nothing but itself and the files it saves, and logs no error.
    Mutations: a fetch of an https URL planted in the page script; an exception thrown in the
    export's draw()."""
    results, _ = ran
    for result in results.values():
        assert [
            u for u in result["requests"] if not u.startswith(("file:", "blob:", "data:"))
        ] == []
        assert result["errors"] == [], result["name"]


def test_linked_view_stacks_at_phone_width(ran):
    """At 375 px the scatter, its reading guide and the parallel sets stack, each inside the
    page, which does not scroll sideways. Mutation: the parallel sets inside the plot column."""
    acted = ran[0]["phone"]["acted"]["stack"]
    plot, reading, below = acted["plot"], acted["reading"], acted["below"]
    assert acted["scroll"] <= acted["page"] == 375
    assert plot[3] <= reading[1] + 1 and reading[3] <= below[1] + 1
    assert all(box[0] >= 0 and box[2] <= acted["page"] for box in (plot, reading, below))


def test_alone_draws_the_page_figure_at_twice_its_size(ran):
    """Measured in the browser: the standalone scatter's tick text and markers twice the
    module page's, and its width twice the page's at full layout (style.PAGE_PLOT_PX) in a
    2,560 px window, wide enough for it (Rytis 2026-10-02); the
    volcano's legend keys twice as long and as thick, its marks twice as large up to plotly's
    16 px cap. Mutations: ALONE_SCALE 1.5; the legend's key length left out of the template,
    so the zoom never reached it."""
    page, alone = ran[0]["page"]["acted"], ran[0]["alone"]["acted"]
    assert alone["tick"] == pytest.approx(TWICE * page["tick"])
    assert (page["marker"], alone["marker"]) == pytest.approx(
        (PLOTLY_MARKER_PX, TWICE * PLOTLY_MARKER_PX))  # fmt: skip
    assert alone["width"] == pytest.approx(TWICE * PAGE_PLOT_PX["beside"])
    (small,), (large,) = ran[0]["page_volcano"]["acted"], ran[0]["alone_volcano"]["acted"]
    for (a, wa), (b, wb) in zip(small["legendlines"], large["legendlines"], strict=True):
        assert (b[2] - b[0], wb) == pytest.approx((TWICE * (a[2] - a[0]), TWICE * wa))
    for a, b in zip(small["legendpoints"], large["legendpoints"], strict=True):
        assert b[2] - b[0] == pytest.approx(min(TWICE * (a[2] - a[0]), LEGEND_MARKER_CAP_PX))
    assert small["legendpoints"][0][2] - small["legendpoints"][0][0] == pytest.approx(
        PLOTLY_MARKER_PX)  # fmt: skip


@pytest.mark.parametrize("name", ["linked", "alone_paint", "volcano", "alone_volcano_phone"])
def test_exported_marks_are_the_figure_files_size(ran, name):
    """A point is MARKER_SIZE pt across in the export, from the module page, from the page
    opened alone at twice the size, and opened alone on a phone, where it is drawn at the
    page's size. Mutations: the page's zoom not divided back out; the zoom left at 2 where a
    figure opened alone is drawn at the page's size."""
    root = _svg(ran[0][name]["paths"]["svg"])
    radii = {float(re.match(r"M([\d.]+),0A", p.get("d")).group(1)) for p in root.iter(f"{SVG}path")
             if p.get("class", "").startswith("point")}  # fmt: skip
    assert radii == {MARKER_SIZE / 2}


def _intersect(a: list[float], b: list[float]) -> bool:
    return a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]


def _inside(inner: list[float], outer: list[float], slack: float = 0.5) -> bool:
    return all((inner[0] >= outer[0] - slack, inner[2] <= outer[2] + slack,
                inner[1] >= outer[1] - slack, inner[3] <= outer[3] + slack))  # fmt: skip


@pytest.mark.parametrize("width", [1280, 1440])
def test_linked_view_spans_the_card_under_the_guide_at_desktop_width(ran, width):
    """The parallel sets start where the scatter does and end where the reading guide does,
    and at 1,440 px they are PAGE_PLOT_PX["across"] wide. Mutation: the parallel sets moved
    into the plot column, beside the guide, which passed a test of the HTML's order alone
    (review, 2026-10-01)."""
    acted = ran[0][f"desktop_{width}"]["acted"]
    plot, reading, below = acted["plot"], acted["reading"], acted["below"]
    assert plot[2] < reading[0], "the guide sits beside the scatter"
    assert below[0] == pytest.approx(plot[0]) and below[2] == pytest.approx(reading[2])
    assert below[1] >= max(plot[3], reading[3])
    if width == 1440:
        assert acted["sets"] == pytest.approx(PAGE_PLOT_PX["across"])


def test_at_phone_width_nothing_runs_into_anything(ran):
    """At 375 px the parallel sets keep their titles apart and scroll sideways, saying so; the
    scatter's and the volcano's x titles stay inside their plots. Mutations: the parallel
    sets' least width dropped, so four titles overlapped by up to 55 px; the page's title fit
    switched off, which clipped "PC1, ... (VST units)" at the screen's edge."""
    acted = ran[0]["phone"]["acted"]
    scatter_plot, sets = acted["linked"]
    dims = [d for d in sets["dims"] if d[2] > d[0]]
    assert len(dims) == 5
    assert not any(_intersect(a, b) for i, a in enumerate(dims) for b in dims[i + 1 :])
    assert sets["hint"] and sets["width"] > sets["box"][2] - sets["box"][0]
    assert not scatter_plot["hint"]
    for plot in (scatter_plot, acted["volcano"][0]):
        (title,) = plot["xtitle"]
        assert _inside(title, plot["svg"]), (title, plot["svg"])


@pytest.mark.parametrize(
    "name", ["bands_1280", "bands_375", "bands_alone", "bands_alone_768", "bands_alone_375"]
)
def test_band_map_title_and_labels_stay_clear(ran, name):
    """Run names in long capitals: the column title sits below every column label, on the
    page at 1,280 and 375 px and opened alone at 1,440, 768 and 375 px; every row label starts
    inside the plot; the rows are a phone's, each as tall as its wrapped label, where the plot
    is narrower than its labels need, and the static's elsewhere. Mutations: the title's drop
    estimated at the average glyph (it overlapped the labels by 18 px); the narrow row labels
    never swapped in (8 of 11 started left of the plot); opened alone, the narrow layout's
    width and height left at the page's, so at 768 px 7 of 11 row labels started up to 198 px
    left of the plot (2026-10-02); the phone's band edges never restyled (merge, 2026-10-07)."""
    (plot,) = ran[0][name]["acted"]
    phone = name in ("bands_375", "bands_alone_768", "bands_alone_375")
    assert plot["narrow"] and (plot["edges"] == plot["narrow"]) == phone
    (title,) = [a for a in plot["annotations"]]
    assert plot["xticks"] and not any(_intersect(title, t) for t in plot["xticks"])
    assert title[1] >= max(t[3] for t in plot["xticks"])
    assert all(_inside(t, plot["svg"]) for t in plot["yticks"])
    assert _inside(title, plot["svg"])


def test_alone_on_a_phone_fits_the_screen_at_the_pages_size(ran):
    """Opened alone in a 375 px window, the band map is drawn at its box's width with no hint,
    at the page's size, its phone layout in: Rytis 2026-10-07, below the phone width Open alone
    fits the screen at the page's size; every row label stays inside the plot
    (``test_band_map_title_and_labels_stay_clear``). Mutation: the floor kept, which drew it
    634 px wide in its 317 px box under the hint."""
    (plot,) = ran[0]["bands_alone_375"]["acted"]
    assert plot["width"] == pytest.approx(plot["box"][2] - plot["box"][0], abs=0.5)
    assert plot["width"] < TWICE * PHONE_PLOT_PX and not plot["hint"]


@pytest.mark.parametrize(("name", "windows"), [("bands_alone_reflow", BANDS_ALONE_WINDOWS),
                                               ("bands_reflow", BANDS_PAGE_WINDOWS)])  # fmt: skip
def test_band_map_reflows_as_the_window_changes(ran, name, windows):
    """Rytis, 2026-10-02: the plot reflows on resize. The band map, whose page script swaps its
    narrow layout in and out, opened alone at 375 px and resized to 1,280, 2,560 and 375 px, is
    drawn at each step at twice the page's width or its box's, whichever is less, and below
    twice a phone's plot at its box's at the page's size (Rytis 2026-10-07), which the steps
    reach in turn; on the module page, opened at 375 px and widened to 1,440, at its box's
    width; every row label stays inside the plot. Mutations: the narrow layout's relayout sent
    without autosize and a null width, so plotly pinned the width it drew at: alone, 634 px in
    a 1,167 px box at 1,280 with 7 of 11 row labels up to 219 px left of the plot, and on the
    module page 317 px in a 686 px box at 1,440 (review, 2026-10-02); autosize sent without the
    null width, which the height's own edit pins again."""
    result = ran[0][name]
    plots = [step["plots"][0] for step in (result["acted"], *result["resized"])]
    assert len(plots) == len(windows)
    floor, cap = TWICE * PHONE_PLOT_PX, TWICE * PAGE_PLOT_PX["beside"]
    alone = name == "bands_alone_reflow"
    for plot in plots:
        expected = min(plot["box"], cap) if alone else plot["box"]
        assert plot["drawn"] == pytest.approx(expected, abs=0.5), (windows, plot["box"])
        assert plot["yticks"] and all(_inside(t, plot["svg"]) for t in plot["yticks"])
    if alone:
        regimes = ["page" if p["box"] < floor else "cap" if p["box"] > cap else "box"
                   for p in plots]  # fmt: skip
        assert regimes == ["page", "box", "cap", "page"]


def _words(text: str) -> set[str]:
    """The words of a text that are not numbers: tick values differ between the renderers."""
    return {
        w
        for w in re.findall(r"[^\s,;:()]+", text.replace("−", "-"))
        if not re.fullmatch(r"-?[\d.]+\.?", w)
    }


def _reference_words(path: Path) -> set[str]:
    return {w for t in _svg(path).iter(f"{SVG}text") for w in _words("".join(t.itertext()))}


def _dashes(path: Path) -> set[str]:
    data = PdfReader(path).pages[0].get_contents().get_data().decode("latin-1")
    return {
        " ".join(f"{float(v):g}" for v in d.split())
        for d in re.findall(r"\[([\d. ]*)\] [\d.]+ d", data)
        if d.strip()
    }


def _points(path: Path) -> int:
    """Marks drawn: a scatter's points and a bar chart's bars, each of class "point"."""
    return sum(1 for e in _svg(path).iter() if "point" in e.get("class", "").split())


@pytest.mark.parametrize("name", ["volcano", "bars"])
def test_at_the_defaults_the_view_is_the_figure_file(ran, built, name):
    """Drawn at the defaults, a view's file holds every word, label, dash and point the
    reference file does: only what the reader moves may differ. Mutations: plotly's own
    dashes ("dot" came out dashed); the gene labels left out; the volcano's thinned points
    left out."""
    out, _, demos, blocks, _ = built
    paths = ran[0]["defaults"]["paths"]
    reference = out / f"{blocks[name].stem}_white"
    words = _words(" ".join(PdfReader(paths[f"{name}_pdf"]).pages[0].extract_text().split()))
    assert _reference_words(reference.with_suffix(".svg")) <= words
    assert _dashes(paths[f"{name}_pdf"]) == _dashes(reference.with_suffix(".pdf"))
    drawn = demos[name].frame.height if name == "volcano" else len(demos[name]._categories())
    assert _points(paths[f"{name}_svg"]) == drawn


def test_a_thinned_volcano_draws_every_point_and_its_files_caption(ran, built):
    """The page draws every 3rd of the points no setting can call; moved, its PDF and SVG draw
    all of them, as the figure file does, and the caption is the figure file's. Mutations: the
    thinned points left out of the export (it drew the page's 4,186 of 8,016 in the review);
    the page's thinning sentence carried into the files."""
    result, data = ran[0]["thin"], built[2]["thin"]
    assert "every 3rd" in built[3]["thin"].caption
    assert len(built[3]["thin"].rows) < data.frame.height
    assert _points(result["paths"]["svg"]) == data.frame.height
    text = " ".join(PdfReader(result["paths"]["pdf"]).pages[0].extract_text().split())
    assert data.caption in text and "This page draws" not in text
    assert result["acted"]["summary"] in text


def test_a_moved_volcano_keeps_the_figure_files_labels(ran, built):
    """Moving the thresholds keeps the gene names the figure file prints. Mutation: the labels
    left out of the export."""
    out, _, demos, blocks, _ = built
    labels = {label["text"] for label in blocks["volcano"].controls["export"]["labels"]}
    reference = _reference_words(out / f"{blocks['volcano'].stem}_white.svg")
    assert labels and labels <= reference
    words = _words(" ".join(PdfReader(ran[0]["volcano"]["paths"]["pdf"]).pages[0]
                            .extract_text().split()))  # fmt: skip
    assert labels <= words


def test_legend_keys_are_opaque_as_in_the_figure_files(ran):
    """A key is drawn at full strength, its points at the page's opacity. Mutation: the key
    left at the points' 0.5."""
    root = _svg(ran[0]["volcano"]["paths"]["svg"])
    keys = [e for g in root.iter(f"{SVG}g") if g.get("class") == "legendpoints" for e in g
            if _paint(e, "fill") in {"#D55E00", "#0072B2"}]  # fmt: skip
    assert len(keys) == 2
    assert all(re.search(r"opacity:\s*1\b", e.get("style", "")) for e in keys)


def test_every_text_of_every_file_is_on_its_page(texts):
    """No word of any file the buttons gave leaves its page: a linked view with outer levels
    longer than the scatter's margins included. Mutation: the parallel sets' margins left out
    of the export, so "bone marrow aspirate" started 16 pt left of the page and "passes every
    band" ended 36 pt right of it (review, 2026-10-01)."""
    assert {"long_svg", "linked_svg", "volcano_svg", "bars_svg", "thin_svg"} <= set(texts)
    for name, measured in texts.items():
        width, height = measured["page"]
        page = [0.0, 0.0, width, height]
        off = [t for t in measured["texts"] if not _inside(t[1:], page)]
        assert off == [], name
    long_words = " ".join(t[0] for t in texts["long_svg"]["texts"])
    assert "bone marrow aspirate" in long_words and "passes every band" in long_words


@pytest.mark.parametrize("kind", ["volcano", "linked"])
@pytest.mark.parametrize("width", [1280, 1512])
def test_alone_keeps_twice_the_type_within_the_window(ran, kind, width):
    """Rytis, 2026-10-02: opened alone, a figure keeps twice the page's type, marks and lines,
    but its plot is no wider than its card, so the page does not scroll sideways. At 1,280 and
    1,512 px windows the double-width volcano and the linked view's scatter are drawn within
    their cards, with no scroll hint, their tick labels at twice the module page's size; on the
    module page no plot scrolls. Mutations: the width left at twice the page's (1,372 px in a
    1,167 px card at 1,280); the type left at the page's on the standalone page."""
    fit, page = ran[0][f"fit_{kind}_{width}"]["acted"], ran[0]["fit_page"]["acted"]
    assert fit["scroll"] <= fit["page"]
    plot, on_page = fit["plots"][0], page["plots"][ON_PAGE[kind][0]]
    assert plot["drawn"] <= plot["box"] + 0.5 and not plot["hint"]
    assert plot["text"] == pytest.approx(TWICE * on_page["text"])
    assert not any(p["hint"] for p in page["plots"]) and page["scroll"] <= page["page"]


def test_alone_takes_twice_the_page_width_in_a_wide_window(ran):
    """At a 2,560 px window each plot opened alone is twice as wide as the module page draws it
    at full layout (1,440 px): the volcano and the scatter 2 x 686 px, the parallel sets 2 x
    1,022 px, with no scroll hint. Mutation: the plot left at its box's width, 2,447 px."""
    page = ran[0]["fit_page"]["acted"]["plots"]
    beside, across = PAGE_PLOT_PX["beside"], PAGE_PLOT_PX["across"]
    assert [round(page[i]["drawn"]) for i in (0, 1, 2)] == [beside, across, beside]
    for kind, numbers in ON_PAGE.items():
        fit = ran[0][f"fit_{kind}_2560"]["acted"]
        for plot, i in zip(fit["plots"], numbers, strict=True):
            assert plot["drawn"] == pytest.approx(TWICE * page[i]["drawn"]) and not plot["hint"]


def test_alone_reflows_when_the_window_narrows(ran):
    """Opened at 2,560 px and narrowed to 1,280, the volcano is redrawn at its card's width,
    narrower than before, with no scroll hint, and the page does not scroll. Mutations:
    plotly's responsive redraw switched off, which kept the plot 1,372 px wide in a 1,167 px
    card; the hint checked only when its box resizes, which left it shown over the redrawn
    plot."""
    before, (after,) = ran[0]["reflow"]["acted"], ran[0]["reflow"]["resized"]
    ((wide,), (narrow,)) = before["plots"], after["plots"]
    assert wide["drawn"] == pytest.approx(TWICE * PAGE_PLOT_PX["beside"])
    assert narrow["drawn"] == pytest.approx(narrow["box"], abs=0.5)
    assert narrow["drawn"] < wide["drawn"]
    assert after["scroll"] <= after["page"] and not narrow["hint"]


def _drawn_text(boxes: list[list[float]]) -> list[list[float]]:
    """The boxes with some width: plotly writes an empty dimension title under every level."""
    return [b for b in boxes if b[2] > b[0]]


def _apart(sets: dict) -> bool:
    """No dimension title over another, and no level over another level, a title or a bar.

    A title sits on its own bar, centred over it, and its line box reaches the bar's top.
    """
    titles, levels = _drawn_text(sets["titles"]), _drawn_text(sets["levels"])
    pairs = [(a, b) for i, a in enumerate(titles) for b in titles[i + 1 :]]
    pairs += [
        (a, b) for i, a in enumerate(levels) for b in [*levels[i + 1 :], *titles, *sets["bars"]]
    ]
    return not any(_intersect(a, b) for a, b in pairs)


def test_alone_parallel_sets_fit_where_their_least_width_allows(ran):
    """The linked view's parallel sets opened alone, at twice the page's type, keep a least
    width of 1,392 px for the demo, which holds twice the page's 16 px gap between a level and
    the bar before it: plotly's bars and insets do not grow with the type. At 1,512 and 2,560
    px they fit their card, with no title or level over another or over a bar. At 1,280 px the
    card holds 1,167: they keep their least width and scroll inside the card under the hint,
    and the page does not scroll. That is the rule until Rytis decides P2Q.13: decision 2 asks
    them to fit, and drawn at the card's width a level runs over a bar at 1,280 but at 1,440
    keeps 15.9 px (D-024). Mutations: the least width doubled whole with the zoom, 1,572 px,
    which scrolled them at 1,512; the least width dropped alone, which drew them 1,167 px wide
    at 1,280; the hint never shown."""
    needs = demo_linked().below().layout.meta["least_width"]
    least = max(TWICE * scaled + fixed for scaled, fixed in needs)
    data = demo_linked()
    shown = (len(data.dimensions), sum(data.frame[c].n_unique() for c in data.dimensions))
    for width in (1512, 2560):
        fit = ran[0][f"fit_linked_{width}"]["acted"]
        sets = fit["plots"][1]
        assert (
            sets["drawn"] <= sets["box"] + 0.5 and not sets["hint"] and fit["scroll"] <= fit["page"]
        )
        drawn = (len(_drawn_text(sets["titles"])), len(_drawn_text(sets["levels"])))
        assert drawn == shown and _apart(sets)
    fit = ran[0]["fit_linked_1280"]["acted"]
    sets = fit["plots"][1]
    assert sets["drawn"] == pytest.approx(least, abs=0.5) and sets["drawn"] > sets["box"]
    assert sets["hint"] and fit["scroll"] <= fit["page"]
    assert _apart(sets)
