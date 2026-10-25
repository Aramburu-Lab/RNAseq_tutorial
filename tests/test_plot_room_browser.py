"""The plot is the largest part of its card, in Chromium (Rytis 2026-10-23, D-117).

Rytis: "The biggest part in the card should be the plot, not the combination of the legend
and the 'How to read' blocks." On the core review pages the plot's drawing area was smaller
than its legend's plus its reading guide's in 28 of 88 cards, at 1,280 and 1,920 px alike. A
module page of the figure types a page draws, a legend-heavy scatter among them, each with a
reading guide, is opened at 375, 1,280 and 1,920 px. Skipped where no Chromium is found; set
``CTX_CHROMIUM``.
"""

import tempfile
from pathlib import Path

import polars as pl
import pytest
from synth import counts, demo_bar, demo_curves, demo_heatmap, demo_surface, demo_volcano
from test_export_browser import CHROMIUM, _browse, _module

from ctx.report.builder import FigureText, write
from ctx.viz.save import figure_block
from ctx.viz.scatter import ScatterData
from ctx.viz.style import PAGE_PLOT_ASPECT, PAGE_PLOT_PX, PHONE_PX, category_colours

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
WIDTHS = (PHONE_PX, 1280, 1920)
# Each visible card's plot drawing (plotly's _size), its legends, its reading guide and the
# card, in px and px²; whether its plot holds its own aspect (a heatmap's cells); and the page.
CARDS = """new Promise(r => setTimeout(r, 3000)).then(() => {
  const opened = Array.from(document.querySelectorAll(".reading details")).map(d => d.open);
  const area = e => {
    if (!e) return 0;
    const b = e.getBoundingClientRect();
    return b.width * b.height;
  };
  const cards = Array.from(document.querySelectorAll("figure.fig[id]")).filter(f => f.offsetParent)
    .map(f => {
      const gd = f.querySelector(".fig-plot .js-plotly-plot"), s = gd._fullLayout._size;
      const legend = Array.from(gd.querySelectorAll(".legend")).reduce((a, l) => a + area(l), 0);
      const plot = gd.querySelector(".nsewdrag").getBoundingClientRect();
      const keys = Array.from(gd.querySelectorAll(".legend, .colorbar")).map(l => {
        const b = l.getBoundingClientRect(); return [b.left, b.top, b.right, b.bottom]; });
      return {title: f.querySelector(".fig-title").textContent, card: area(f),
              plot: s.w * s.h, plotw: s.w, ploth: s.h, legend: legend,
              reading: area(f.querySelector(".reading")), keys: keys,
              box: f.querySelector(".fig-plot").clientWidth,
              area: [plot.left, plot.top, plot.right, plot.bottom],
              held: gd.data.some(t => t.type === "heatmap")};
    });
  const html = document.documentElement;
  return {cards: cards, box: html.clientWidth, scroll: html.scrollWidth, opened: opened};
})"""
# Every guide opened as the page loads, before its plots settle, then CARDS (D-119).
OPEN = 'document.querySelectorAll(".reading details").forEach(d => { d.open = true; }); ' + CARDS
LONG = (
    "tumour-associated macrophages",
    "tissue-resident alveolar macrophages",
    "classical monocytes from blood",
    "non-classical monocytes from blood",
    "low-density neutrophils",
    "normal-density neutrophils",
    "bone-marrow precursors",
)
READING = ("Each dot is one sample, placed by its two scores.",
           "Colour is the cell type the sample sheet gives; shape is the batch.",
           "Samples of one type sit together where the scores separate the types.",
           "A sample far from its own type's cloud is worth a second look at its QC.",
           "The axes are in the units of the scores, not in reads or genes.")  # fmt: skip


def _crowded() -> ScatterData:
    """Seven colour levels with long names by two shapes, 14 legend entries: the legend of
    Step-1's mitochondrial share scatter (13 entries, 136,756 px² against a 95,480 px² plot)."""
    flags = ("not flagged on either", "flagged on either share")
    rows = [(f"run{i:03d}", LONG[i % 7], flags[i % 2], float(i % 13), float(i % 11))
            for i in range(70)]  # fmt: skip
    frame = pl.DataFrame(rows, schema=["run", "cell", "batch", "x", "y"], orient="row")
    columns = {"run": ("none", "Run."), "cell": ("none", "Cell type."), "batch": ("none", "Batch."),
               "x": ("%", "Mitochondrial share."), "y": ("%", "rRNA share.")}  # fmt: skip
    return ScatterData(frame=frame, columns=columns, x="x", y="y", x_label="Mitochondrial share",
                       y_label="Nuclear rRNA share", colour="cell",
                       colours=category_colours(list(LONG)), symbol="batch",
                       caption="Synthetic runs; no test.", hover=("run",))  # fmt: skip


def room_page(folder: Path) -> Path:
    """A module page of a legend-heavy scatter, curves, bars, a volcano, a surface and a
    heatmap, each with a reading guide but the volcano, which takes its own from its calls."""
    synthetic = counts()
    figures = [("crowded", _crowded()), ("curves", demo_curves()), ("bars", demo_bar(synthetic)),
               ("surface", demo_surface()), ("heatmap", demo_heatmap())]  # fmt: skip
    blocks = [figure_block(data, f"2026-10-23_demo_{name}_v1", FigureText(name.title(),
              "A lead.", READING)) for name, data in figures]  # fmt: skip
    blocks.append(figure_block(demo_volcano(synthetic), "2026-10-23_demo_volcano_v1",
                               FigureText("Volcano", "A lead.")))  # fmt: skip
    return write(_module(blocks), folder / "2026-10-23_demo_room_v1.html")


def measure(page: Path, widths=WIDTHS, act: str = CARDS) -> dict[int, dict]:
    """``act`` on ``page`` at each width, by width."""
    tasks = [{"name": str(w), "url": page.as_uri(), "width": w, "height": 1000, "act": act}
             for w in widths]  # fmt: skip
    with tempfile.TemporaryDirectory() as folder:
        seen = _browse(tasks, Path(folder))
    return {int(name): result["acted"] for name, result in seen.items()}


@pytest.fixture(scope="module")
def measured(tmp_path_factory) -> dict[int, dict]:
    """The room page measured at every width."""
    return measure(room_page(tmp_path_factory.mktemp("room")))


def test_reading_guides_load_folded_and_opening_one_keeps_the_plot(measured, tmp_path):
    """Every guide loads closed; opened, it is counted folded, so the plot keeps its size.

    Catches: a guide that loads open; the fit counting an opened guide's list (D-119)."""
    assert measured[1280]["opened"] == [False] * 6, measured[1280]["opened"]
    opened = measure(room_page(tmp_path), (1280,), OPEN)[1280]
    assert opened["opened"] == [True] * 6
    for closed, open_ in zip(measured[1280]["cards"], opened["cards"], strict=True):
        assert abs(open_["plot"] - closed["plot"]) <= 0.01 * closed["plot"], (closed, open_)


def test_the_plot_outweighs_its_legend_and_reading_guide(measured):
    """At 1,280 and 1,920 px every card's plot is larger than its legend and guide together,
    and every plot that holds no aspect of its own is at least 9/16 of its width tall.

    Catches: the reading guide back beside the plot; the 16:9 floor removed from the fit."""
    for width in (1280, 1920):
        cards = measured[width]["cards"]
        assert len(cards) == 6, width
        for c in cards:
            assert c["plot"] > c["legend"] + c["reading"], (width, c)
            if not c["held"]:
                assert c["ploth"] >= PAGE_PLOT_ASPECT * c["plotw"] - 1, (width, c)


def test_legends_sit_above_or_below_the_plot(measured):
    """No legend or colour bar of a page plot starts right of the plot's left edge and beside
    it: each lies wholly above or below the plot's drawing area.

    Catches: a legend or colour bar back at the plot's side (the surface's colour bar was)."""
    for width in (1280, 1920):
        for c in measured[width]["cards"]:
            left, top, right, bottom = c["area"]
            for key in c["keys"]:
                assert key[3] <= top + 1 or key[1] >= bottom - 1, (width, c["title"], key)


def test_a_wide_screen_gives_every_plot_more_room(measured):
    """Each card's plot draws larger at 1,920 px than at 1,280 px, in a box wider than the
    1,072 px column gave it (PAGE_PLOT_PX).

    Catches: the content column fixed at one width again (a 1,920 px window then still drew
    larger plots than 1,280 px, whose contents rail leaves 968 px)."""
    for narrow, wide in zip(measured[1280]["cards"], measured[1920]["cards"], strict=True):
        assert wide["plot"] > narrow["plot"], (narrow, wide)
        assert wide["box"] > PAGE_PLOT_PX["across"] + 1, wide


def test_a_phone_never_scrolls_sideways(measured):
    """At 375 px the page is no wider than its window."""
    phone = measured[PHONE_PX]
    assert phone["box"] == PHONE_PX and phone["scroll"] <= phone["box"], phone
