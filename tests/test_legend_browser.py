"""A bar figure's legend and labels at a phone's width in Chromium: each inside the plot's box.

M1's depth legend ran past its plot at 375 px: "ENCODE, older projects: 20 million" with its
key took 355 px of a 317 px plot (work-order item 2's leftover, measured 2026-10-14 on
Step-1's page: the entry's right edge at 384 px in a 375 px window). Skipped where no Chromium
is found; set ``CTX_CHROMIUM``.
"""

from dataclasses import replace
from datetime import date

import polars as pl
import pytest
from synth import counts, demo_bar
from test_export_browser import CHROMIUM, _browse, _module

from ctx.io.naming import stem
from ctx.report.builder import FigureText, write
from ctx.viz.axes import phone_legend
from ctx.viz.bar import BarData
from ctx.viz.save import figure_block, write_figure
from ctx.viz.style import category_colours

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
LINES = ((20.0, "ENCODE, older projects: 20 million"), (30.0, "ENCODE: 30 million aligned reads"))
READ = """new Promise(r => setTimeout(r, 1500)).then(() => {
  const d = document.querySelector('.js-plotly-plot'), box = d.getBoundingClientRect();
  return {div: [box.left, box.right], view: document.documentElement.clientWidth,
          items: Array.from(d.querySelectorAll('.legend .traces')).map(t => {
            const b = t.getBoundingClientRect(); return [t.textContent, b.left, b.right]; })};
})"""


@pytest.fixture(scope="module")
def ran(tmp_path_factory):
    out = tmp_path_factory.mktemp("legend")
    data = replace(demo_bar(counts()), hline=None, lines=LINES)
    figure_stem = stem("gather", "legend_demo", 1, day=date(2026, 10, 14))
    text = FigureText("Depth", reading=("Each bar is one made-up group.",))
    write_figure(data, out, figure_stem, title="Depth", data_releases={}, text=text)
    page = out / f"{figure_stem}.html"
    tasks = [{"name": str(w), "url": page.as_uri(), "width": w, "height": 900, "act": READ}
             for w in (375, 1280)]  # fmt: skip
    return data, _browse(tasks, tmp_path_factory.mktemp("downloads"))


# Each y tick label's box against the plot's svg, on a module page, once the page settles.
TICKS = """new Promise(r => setTimeout(r, 2000)).then(() => {
  const d = document.querySelector('.js-plotly-plot');
  const box = d.querySelector('svg').getBoundingClientRect();
  return {svg: [box.left, box.right], view: document.documentElement.clientWidth,
          scroll: document.documentElement.scrollWidth,
          ticks: Array.from(d.querySelectorAll('.ytick text')).map(t => {
            const b = t.getBoundingClientRect();
            return [t.textContent, t.querySelectorAll('tspan.line').length, b.left, b.right]; })};
})"""
LONG = ("Neutrophils, PDAC against steady state, no batch term",
        "Neutrophils, PDAC against steady state, batch in the design")  # fmt: skip


def test_long_bar_labels_wrap_inside_the_plot_on_a_phone(tmp_path):
    """E2b: at 375 px M5's genes-called bars set their contrast names left of the plot, where
    every one ran off its edge. Each label now ends inside the plot's svg at 375 and 1,280 px,
    wrapped onto lines at 375; the figure files keep the labels as built.

    Mutations caught: no narrow labels; labels wrapped wider than half a phone's plot. A PDF
    of a moved view prints the labels as built: tests/test_phone_export_browser.py."""
    frame = pl.DataFrame({"contrast": [*LONG, *LONG], "direction": ["up"] * 2 + ["down"] * 2,
                          "n": [500, 20, 200, 10]})  # fmt: skip
    data = BarData(
        frame=frame, columns={"contrast": ("none", "Contrast."), "direction": ("none", "Way."),
                              "n": ("genes", "Genes called.")},
        category="contrast", value="n", category_label="Contrast", value_label="Genes called",
        series="direction", colours=category_colours(["up", "down"]), caption="Counts, no test.",
        horizontal=True, stacked=True)  # fmt: skip
    built = data.interactive()
    assert list(built.layout.yaxis.ticktext) == list(LONG)
    figure_stem = stem("de", "labels_demo", 1, day=date(2026, 10, 19))
    text = FigureText("Genes called", reading=("Each bar is one made-up contrast.",))
    write_figure(data, tmp_path, figure_stem, title="Genes called", data_releases={}, text=text)
    page = write(
        _module([figure_block(data, figure_stem, text)]),
        tmp_path / f"{stem('de', 'labels_overview', 1, day=date(2026, 10, 19))}.html",
    )
    tasks = [{"name": str(w), "url": page.as_uri(), "width": w, "height": 900, "act": TICKS}
             for w in (375, 1280)]  # fmt: skip
    seen = _browse(tasks, tmp_path / "downloads")
    for width, result in seen.items():
        acted = result["acted"]
        assert result["errors"] == [] and acted["scroll"] == acted["view"], (width, acted)
        assert len(acted["ticks"]) == len(LONG), acted
        left, right = acted["svg"]
        assert all(left - 0.5 <= a and b <= right + 0.5 for _, _, a, b in acted["ticks"]), acted
    assert all(lines > 1 for _, lines, _, _ in seen["375"]["acted"]["ticks"])


def test_a_phone_wraps_each_legend_entry_inside_the_plot(ran):
    """At 375 px every entry, the lines' included, ends inside the plot's box and the window;
    at 1,280 px the entries keep one line each.

    Mutations caught: no phone legend; the lines' shape names left unwrapped; the wrap room
    without the key's width."""
    data, results = ran
    assert phone_legend(data.interactive()) is not None
    for width in ("375", "1280"):
        seen = results[width]
        assert seen["errors"] == []
        acted = seen["acted"]
        names = {n for n, _, _ in acted["items"]}
        assert (
            {n for _, n in LINES} <= names
            if width == "1280"
            else names.isdisjoint({n for _, n in LINES})
        )  # wrapped: the words run on without the line break
        right = min(acted["div"][1], acted["view"])
        assert all(r <= right + 0.5 for _, _, r in acted["items"]), acted
