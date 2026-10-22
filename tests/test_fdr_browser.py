"""M6's FDR cutoff (P2S.5, D-109) in Chromium, on the pipeline fixture's M6 page: the live
ORA, the GSEA bars and one set's NES in every contrast, each at three cutoffs and back to the
build's, at 1,280 and 375 px in both themes, then each figure's PDF and SVG at a cutoff.

What plotly drew is checked against the contracts ``ctx function`` builds at that cutoff and
the words its page functions write there. The twins themselves are tested under node in
``test_fdr_page.py``. Each test names the mutations it catches.
"""

import html
import math
import re
from datetime import date

import polars as pl
import pytest
from conftest import DE_SPEC, _ran, once
from synth import DE_GROUPS, design_counts, write_gather_outputs, write_gene_sets
from test_export_browser import CHROMIUM, _browse

from ctx.dataset import load
from ctx.de.page import every
from ctx.function import live, run
from ctx.function.page import gsea_text, nes_text, ora_text, set_words
from ctx.report import text

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
SIZES = {"wide_light": (1280, "light"), "wide_dark": (1280, "dark"),
         "phone_light": (375, "light"), "phone_dark": (375, "dark")}  # fmt: skip
CUTS = ("0.1", "0.01", "0.001", "0.05")
# Each figure the cutoff marks, by the form that holds its select.
FIGURES = {"ora": "form.ora", "gsea": "form.view:has(select[name=fdr])", "nes": "form.fdr"}
# The same figures as the driver finds them to press their buttons (no :has inside :has).
CARDS = {
    "ora": ".fig:has(form.ora)",
    "gsea": ".fig:has(form.view select[name=fdr])",
    "nes": ".fig:has(form.fdr)",
}
READ = """
  const read = (fig) => {
    const d = fig.querySelector('.js-plotly-plot');
    const words = s => (fig.querySelector(s) || {}).textContent || '';
    return {traces: d.data.filter(t => t.type === 'bar').map(t => ({name: t.name,
              y: Array.from(t.y), x: Array.from(t.x), colour: t.marker.color})),
            shapes: (d.layout.shapes || []).map(s => ({x0: s.x0, x1: s.x1, name: s.name})),
            legend: d.layout.showlegend, lead: words('.fig-lead'), caption: words('figcaption'),
            reading: Array.from(fig.querySelectorAll('.reading li')).map(l => l.textContent),
            live: words('[aria-live]'), moved: !!fig.querySelector('a[data-export][role=button]'),
            ticks: Array.from(d.querySelectorAll('.yaxislayer-above .ytick text')).map(
              t => t.getBoundingClientRect().left)};
  };"""
WALK = (
    """(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(2000);
  READ
  const out = {figures: {}};
  for (const [name, sel] of Object.entries(FIGURES)) {
    const form = document.querySelector(sel), fig = form.closest('.fig');
    fig.scrollIntoView({block: 'start', behavior: 'instant'});
    await wait(150);
    const seen = {loaded: read(fig), cuts: {}};
    for (const cut of CUTS) {
      form.elements.fdr.value = cut;
      form.elements.fdr.dispatchEvent(new Event('change', {bubbles: true}));
      await wait(250); seen.cuts[cut] = read(fig);
    }
    out.figures[name] = seen;
  }
  out.box = document.documentElement.clientWidth; out.scroll = document.documentElement.scrollWidth;
  return out;
})()""".replace("READ", READ)
    .replace("FIGURES", repr(FIGURES))
    .replace("CUTS", repr(list(CUTS)))
)
# One figure set to 0.01, whose PDF and SVG buttons the driver then presses.
AT = """(async () => {
  await new Promise(r => setTimeout(r, 2000));
  const form = document.querySelector(SEL);
  form.elements.fdr.value = '0.01';
  form.elements.fdr.dispatchEvent(new Event('change', {bubbles: true}));
  await new Promise(r => setTimeout(r, 300));
  return form.closest('.fig').id;
})()"""


def _root(root):
    """de and function on the synthetic design, with twelve more sets holding 1 to 12 of
    UP_SET's genes among null ones: their FDRs fall between the presets, so each cutoff moves
    some set across the mark. The pipeline fixture's sets sit at FDR 0 or far above 0.1."""
    write_gather_outputs(root, "desynth", design_counts(DE_GROUPS, batch_lfc=2.0), date.today())
    gmt = root / load(DE_SPEC).resource("sets").dest
    write_gene_sets(gmt)
    up = next(line.split("\t")[2:] for line in gmt.read_text().splitlines()
              if line.startswith("UP_SET\t"))  # fmt: skip
    null = [f"G{i:05d}" for i in range(400, 900)]
    mixed = ["\t".join([f"MIXED_{k:02d}", "synthetic", *up[:k], *null[k * 30:k * 30 + 40 - k]])
             for k in range(1, 13)]  # fmt: skip
    gmt.write_text(gmt.read_text() + "\n".join(mixed) + "\n")
    for module in ("de", "function"):
        _ran(root, module, "--dataset", "desynth", "--dataset-file", str(DE_SPEC))
    return root


@pytest.fixture(scope="module")
def m6(tmp_path_factory):
    """(the M6 folder, the walks by size, and the page)."""
    root = once(tmp_path_factory, "fdr_cutoff", _root)
    folder = root / "outputs" / "06_function"
    (page,) = folder.glob("*_function_desynth_overview_v1.html")
    tasks = [{"name": name, "url": page.as_uri(), "width": width, "height": 900, "scheme": scheme,
              "act": WALK} for name, (width, scheme) in SIZES.items()]  # fmt: skip
    return folder, _browse(tasks, tmp_path_factory.mktemp("downloads")), page


def _newest(folder, what: str) -> pl.DataFrame:
    (path,) = folder.glob(f"*_desynth_{what}_v1.parquet")
    return pl.read_parquet(path)


def _words():
    spec = load(DE_SPEC)
    return spec, every(spec.contrasts, spec.columns)


def _drawn(traces: list[dict]) -> list[tuple[str, tuple[str, ...], str]]:
    return [(t["name"], tuple(map(str, t["y"])), t["colour"]) for t in traces]


def _plotted(contract) -> list[tuple[str, tuple[str, ...], str]]:
    return [(t.name, tuple(map(str, t.y)), t.marker.color) for t in contract.interactive().data]


@pytest.mark.parametrize("size", list(SIZES))
def test_the_ora_redraws_each_cutoff_as_live_bars_draws_it(m6, size):
    """The live ORA at each cutoff draws ``live.bars``' sets at that FDR, each in its family's
    colour or grey, the dashed line at -log10 of it named for it, the lead, caption, reading
    line and live line in its words; back at 0.05, the build's.

    Mutations caught: the cutoff read nowhere (the build's bars at every cutoff); the line left
    at 0.05; the legend named for 0.05; the lead left as built; the export switch left off."""
    folder, results, _ = m6
    seen = results[size]
    assert seen["errors"] == []
    acted = seen["acted"]
    assert acted["scroll"] <= acted["box"], (acted["scroll"], acted["box"])
    spec, said = _words()
    featured = spec.de.featured
    table = _newest(folder, f"{featured}_ora")
    background = len(live.symbols(_newest(folder.parents[0] / "05_de", f"{featured}_results")))
    walk = acted["figures"]["ora"]
    families = text.raw("function", "live.families")
    counts = []
    for cut in CUTS:
        alpha, page = float(cut), walk["cuts"][cut]
        bars = live.bars(table, alpha=alpha)
        drawn = {(t["name"], label) for t in page["traces"] for label in t["y"]}
        key = {True: "live.legend.passing", False: "live.legend.failing"}
        want = {(text.fill("function", key[p], family=families[d], alpha=alpha), set_words(s))
                for d, s, p in zip(bars["direction"], bars["gene_set"], bars["passes"],
                                   strict=True)}  # fmt: skip
        assert drawn == want, cut
        (line,) = page["shapes"]
        assert (
            math.isclose(line["x0"], -math.log10(alpha), rel_tol=1e-12) and line["x0"] == line["x1"]
        )
        assert line["name"] == f"FDR {cut}"
        up, down = (int((table.filter(pl.col("direction") == d)["fdr"] < alpha).sum())
                    for d in ("up", "down"))  # fmt: skip
        words = ora_text(said[featured], background, up, down, int((~bars["passes"]).sum()), alpha)
        assert page["lead"] == words.lead and page["reading"][-1] == words.reading[-1], cut
        assert page["caption"] == run._ora_caption(said[featured], background, alpha)  # noqa: SLF001
        assert f"set results pass FDR < {cut}, BH per family" in page["live"]
        assert page["moved"] == (alpha != live.ALPHA), cut
        counts.append(int(bars["passes"].sum()))
    loaded, back = walk["loaded"], walk["cuts"]["0.05"]
    assert not loaded["moved"] and (back["lead"], back["caption"]) == (
        loaded["lead"],
        loaded["caption"],
    )
    assert sorted(_drawn(back["traces"])) == sorted(_drawn(loaded["traces"]))
    assert len(set(counts)) > 1, "the fixture must pass more sets at one cutoff than another"


@pytest.mark.parametrize("size", list(SIZES))
def test_the_gsea_and_nes_bars_mark_each_cutoff_as_the_build_marks_it(m6, size):
    """The GSEA bars and one set's NES at each cutoff: the traces ``_gsea_bars`` and
    ``_featured_bars`` draw there, the GSEA lead and both reading lines in their words, the
    live line naming the cutoff; back at 0.05 the figure as built. No label left of the screen.

    Mutations caught: <= in the twin (a set at FDR 0.05 highlighted); the lead counting the
    bars drawn; the reading line left as built; the live line left out of the export's caption
    line; a cutoff moved back leaving the last marks."""
    folder, results, _ = m6
    acted = results[size]["acted"]
    assert results[size]["errors"] == []
    spec, said = _words()
    featured = spec.de.featured
    gsea = _newest(folder, f"{featured}_gsea")
    tables = {
        n: _newest(folder, f"{n}_gsea")
        for n in said
        if list(folder.glob(f"*_desynth_{n}_gsea_v1.parquet"))
    }
    moves = set()
    for name in ("gsea", "nes"):
        walk = acted["figures"][name]
        for cut in CUTS:
            alpha, page = float(cut), walk["cuts"][cut]
            if name == "gsea":
                contract = run._gsea_bars(said[featured], gsea, alpha)  # noqa: SLF001
                shown = (int((contract.frame["fdr"] < alpha).sum()), contract.frame.height)
                words = gsea_text(said[featured], int((gsea["fdr"] < alpha).sum()), gsea.height,
                                  shown, alpha)  # fmt: skip
                if alpha != live.ALPHA:
                    assert page["lead"] == words.lead, cut
            else:
                contract = run._featured_bars(spec.function.featured_set, tables, said, alpha)  # noqa: SLF001
                words = nes_text(spec.function.featured_set, len(tables), alpha)
            assert sorted(_drawn(page["traces"])) == sorted(_plotted(contract)), (name, cut)
            assert page["reading"][-1] == words.reading[-1], (name, cut)
            assert page["moved"] == (alpha != live.ALPHA)
            assert (f"Sets marked at FDR < {cut};" in page["live"]) == (alpha != live.ALPHA)
            assert all(left >= 0 for left in page["ticks"]), (name, cut, page["ticks"])
            moves.add((name, tuple(len(t["y"]) for t in page["traces"])))
        back, loaded = walk["cuts"]["0.05"], walk["loaded"]
        assert _drawn(back["traces"]) == _drawn(loaded["traces"]), name
        assert (back["lead"], back["reading"]) == (loaded["lead"], loaded["reading"]), name
    # The fixture's premise: some cutoff moves a bar between the marks on each figure.
    assert {n for n, _ in moves} == {"gsea", "nes"} and len(moves) > 2, moves  # noqa: PLR2004


@pytest.mark.parametrize("name", list(FIGURES))
def test_each_figure_exports_the_view_at_its_cutoff(m6, name, tmp_path):
    """At FDR 0.01 each figure's PDF and SVG are the view: the legend names 0.01, the ORA's
    line too and its caption the cutoff, and the live line names it.

    Mutations caught: the export ignoring the cutoff (the build's files, or a moved view drawn
    with the build's caption); the ORA's export caption left at 0.05."""
    _, _, page = m6
    task = {
        "name": "export",
        "url": page.as_uri(),
        "width": 1280,
        "height": 900,
        "act": AT.replace("SEL", repr(FIGURES[name])),
        "figure": CARDS[name],
    }
    downloads = tmp_path / "downloads"
    # A button left on the figure files downloads no view, and the driver stops waiting for it:
    # its own browser, so that stop fails this test and leaves the walks standing.
    try:
        (seen,) = _browse([task], downloads).values()
    except AssertionError as stopped:
        pytest.fail(f"no view downloaded: {str(stopped)[-300:]}")
    assert seen["errors"] == []
    files = seen["files"]
    assert files["svg"]["name"].endswith("_view.svg") and files["pdf"]["name"].endswith("_view.pdf")
    svg = (downloads / files["svg"]["guid"]).read_text()
    legend = [html.unescape(s) for s in re.findall(r'class="legendtext"[^>]*>([^<]*)', svg)]
    words = html.unescape(" ".join(re.findall(r">([^<>]+)</text>", svg)))
    if name == "ora":
        assert "FDR 0.01" in legend and any("FDR < 0.01 (BH)" in n for n in legend), legend
        assert "FDR = 0.01" in words and "FDR = 0.05" not in words
        assert "pass FDR < 0.01" in words
    else:
        assert "FDR < 0.01" in legend and "FDR < 0.05" not in legend, legend
        assert "Sets marked at FDR < 0.01" in words
    assert (downloads / files["pdf"]["guid"]).read_bytes().startswith(b"%PDF")
