"""A bar figure's legend at a phone's width in Chromium: every entry inside the plot's box.

M1's depth legend ran past its plot at 375 px: "ENCODE, older projects: 20 million" with its
key took 355 px of a 317 px plot (work-order item 2's leftover, measured 2026-10-14 on
Step-1's page: the entry's right edge at 384 px in a 375 px window). Skipped where no Chromium
is found; set ``CTX_CHROMIUM``.
"""

from dataclasses import replace
from datetime import date

import pytest
from synth import counts, demo_bar
from test_export_browser import CHROMIUM, _browse

from ctx.io.naming import stem
from ctx.report.builder import FigureText
from ctx.viz.bar import phone_legend
from ctx.viz.save import write_figure

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
