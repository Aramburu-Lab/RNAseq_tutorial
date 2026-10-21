"""A view exported on a phone draws its names as built, in Chromium (TODO.md L.7).

Below a phone's width the page wraps long legend entries and bar labels with ``<br>`` (the
fit's ``meta.narrow``); a view's PDF and SVG are the figure files' width, so they draw the
names as built whatever the screen. Three moved figures on one module page, exported at 375
and 1,280 px: bars filtered, their reference lines' legend names long; horizontal bars, by a
view, their labels long; a PCA, by a view, its legend Step-1's level names. Skipped where no
Chromium is found; set ``CTX_CHROMIUM``.
"""

from dataclasses import replace
from datetime import date
from pathlib import Path
from xml.etree import ElementTree

import polars as pl
import pytest
from test_export_browser import CHROMIUM, SVG, _browse, _module
from test_legend_browser import LINES, LONG
from test_scatter_legend_browser import _long_pca
from test_viz import filtered_bars

from ctx.io.naming import stem
from ctx.report.builder import FigureText, write
from ctx.report.views import viewed
from ctx.viz.bar import BarData
from ctx.viz.save import figure_block, write_figure
from ctx.viz.style import category_colours

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
DAY = date(2026, 10, 21)
WIDTHS = (375, 1280)
# Each control moved as a reader moves it, then the names the page draws, once settled.
MOVE = """(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(1500);
  const fig = document.querySelector("%(figure)s"), d = fig.querySelector(".js-plotly-plot");
  const filter = fig.querySelector("select[data-column]"), view = fig.querySelector("form.view");
  if (filter) {
    filter.value = "LDN"; filter.dispatchEvent(new Event("change", {bubbles: true}));
  } else {
    view.elements.view.value = "1";
    view.elements.view.dispatchEvent(new Event("change", {bubbles: true}));
  }
  await wait(1500);
  return {moved: !!fig.querySelector("a[data-export][role=button]"),
          legend: Array.from(d.querySelectorAll(".legendtext")).map(t => t.textContent),
          names: d.data.map(t => t.name).concat((d.layout.shapes || []).map(s => s.name)),
          ticks: (d.layout.yaxis || {}).ticktext || null};
})()"""


def _labels(n: tuple[int, ...]) -> BarData:
    """Horizontal stacked bars of two contrasts named as long as M5's genes-called labels."""
    frame = pl.DataFrame({"contrast": [*LONG, *LONG], "direction": ["up"] * 2 + ["down"] * 2,
                          "n": list(n)})  # fmt: skip
    return BarData(
        frame=frame, columns={"contrast": ("none", "Contrast."), "direction": ("none", "Way."),
                              "n": ("genes", "Genes called.")},
        category="contrast", value="n", category_label="Contrast", value_label="Genes called",
        series="direction", colours=category_colours(["up", "down"]), caption="Counts, no test.",
        horizontal=True, stacked=True)  # fmt: skip


def _views() -> dict[str, list]:
    """Each viewed figure's options, (name, contract, words), the build first."""
    words = FigureText("Genes called", reading=("Each bar is one made-up contrast.",))
    other_words = FigureText("Samples", reading=("Each dot is one made-up sample.",))
    pca = _long_pca()
    other = replace(_long_pca(), frame=pca.frame.with_columns(pl.col("pc1") * -1),
                    caption="The same samples mirrored; no test.")  # fmt: skip
    return {
        "labels": [("build", _labels((500, 20, 200, 10)), words),
                   ("other", replace(_labels((50, 200, 20, 100)), caption="Others; no test."),
                    words)],
        "pca": [("build", pca, other_words), ("mirrored", other, other_words)],
    }  # fmt: skip


@pytest.fixture(scope="module")
def seen(tmp_path_factory):
    """(each figure's names as built, by figure; the files each width exported, by task)."""
    out = tmp_path_factory.mktemp("phone_export")
    bars = replace(filtered_bars(), lines=LINES)
    text = FigureText("Bars", reading=("Each bar is one made-up run.",))
    write_figure(bars, out, stem("de", "bars_phone", 1, day=DAY), title="Bars",
                 data_releases={}, text=text)  # fmt: skip
    blocks = [figure_block(bars, stem("de", "bars_phone", 1, day=DAY), text)]
    built = {"bars": bars.interactive()}
    for name, options in _views().items():
        figure_stem = stem("de", f"{name}_phone", 1, day=DAY)
        data, words = options[0][1], options[0][2]
        write_figure(data, out, figure_stem, title=name, data_releases={}, text=words)
        blocks.append(viewed(figure_block(data, figure_stem, words), options, "Which"))
        built[name] = options[1][1].interactive()
    page = write(_module(blocks), out / f"{stem('de', 'phone_overview', 1, day=DAY)}.html")
    downloads = tmp_path_factory.mktemp("downloads")
    tasks = [{"name": f"{name}_{width}", "url": page.as_uri(), "width": width, "height": 900,
              "figure": f"#figure-{i}", "act": MOVE % {"figure": f"#figure-{i}"}}
             for i, name in enumerate(built, 1) for width in WIDTHS]  # fmt: skip
    results = _browse(tasks, downloads)
    for result in results.values():
        result["svg"] = downloads / result["files"]["svg"]["guid"]
        result["pdf"] = downloads / result["files"]["pdf"]["guid"]
    return built, results


def _texts(path: Path, cls: str) -> list[list[str]]:
    """Each text of class ``cls`` in an SVG file, as its lines: plotly writes a line per tspan."""
    out = []
    for text in ElementTree.parse(path).getroot().iter(f"{SVG}text"):
        if cls in text.get("class", "").split():
            lines = [s.text or "" for s in text.iter(f"{SVG}tspan") if s.get("class") == "line"]
            out.append(lines or ["".join(text.itertext())])
    return out


def _ticks(path: Path) -> list[list[str]]:
    """Each y tick label of an SVG file, as its lines."""
    out = []
    for g in ElementTree.parse(path).getroot().iter(f"{SVG}g"):
        if g.get("class") == "ytick":
            for text in g.iter(f"{SVG}text"):
                lines = [s.text or "" for s in text.iter(f"{SVG}tspan") if s.get("class") == "line"]
                out.append(lines or ["".join(text.itertext())])
    return out


def _legend_names(figure) -> set[str]:
    traces = {t.name for t in figure.data if t.name and t.showlegend is not False}
    return traces | {s.name for s in figure.layout.shapes if s.name and s.showlegend}


@pytest.mark.parametrize("name", ["bars", "labels", "pca"])
def test_a_view_exported_on_a_phone_draws_its_names_as_built(seen, name):
    """L.7: at 375 px the page wraps the names (``<br>`` in its plot), and the PDF and SVG of
    the moved view draw each legend entry and bar label as built, one line each, the texts the
    same view exports at 1,280 px.

    Mutation caught: the export drawing the names as the phone's layout last set them."""
    built, results = seen
    phone, desktop = results[f"{name}_375"], results[f"{name}_1280"]
    for result in (phone, desktop):
        assert result["errors"] == [] and result["acted"]["moved"], result["acted"]
        assert result["pdf"].read_bytes().startswith(b"%PDF")
    acted = phone["acted"]
    wrapped = [n for n in acted["names"] + (acted["ticks"] or []) if n and "<br>" in n]
    assert wrapped, acted  # the phone's layout is in force on the page
    legend, ticks = (_texts(phone["svg"], "legendtext"), _ticks(phone["svg"]))
    assert legend == _texts(desktop["svg"], "legendtext") and ticks == _ticks(desktop["svg"])
    assert all(len(lines) == 1 for lines in legend + ticks), (legend, ticks)
    assert {lines[0] for lines in legend} == _legend_names(built[name])
    if name == "labels":
        assert sorted(lines[0] for lines in ticks) == sorted(LONG)
