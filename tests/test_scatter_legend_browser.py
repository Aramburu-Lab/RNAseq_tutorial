"""A scatter's legend at a phone's width in Chromium: each entry inside its plot's svg.

At 375 px M5's PCA legend ran past its plot on Step-1: an entry's right edge at 407 px in a
317 px svg (TODO.md L.2), a condition and a read-length batch named together. The PCA here
carries names as long, beside a volcano and an MA plot on one module page. Skipped where no
Chromium is found; set ``CTX_CHROMIUM``.
"""

from dataclasses import replace
from datetime import date

import polars as pl
import pytest
from synth import counts, demo_ma, demo_pca, demo_volcano
from test_export_browser import CHROMIUM, _browse, _module

from ctx.io.naming import stem
from ctx.report.builder import FigureText, write
from ctx.viz.save import figure_block, write_figure
from ctx.viz.style import category_colours

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
# Step-1's PCA levels: a condition, and the read-length batch its points are shaped by.
CELLS = {"type_a": "steady_state", "type_b": "hsct_early"}
BATCHES = {"control": "Reads ≤ 76 nt (min 35 nt)", "case": "Reads > 76 nt (max 500 nt)"}
# Each figure's legend entries against its svg, and the page against the window, once settled.
READ = """new Promise(r => setTimeout(r, 2500)).then(() => ({
  view: document.documentElement.clientWidth, scroll: document.documentElement.scrollWidth,
  figures: Array.from(document.querySelectorAll('.js-plotly-plot')).map(d => {
    const box = d.querySelector('svg').getBoundingClientRect();
    return {svg: [box.left, box.right],
            names: d.data.filter(t => t.showlegend !== false).map(t => t.name),
            items: Array.from(d.querySelectorAll('.legend .traces')).map(t => {
              const b = t.getBoundingClientRect();
              return [t.textContent, t.querySelectorAll('tspan.line').length, b.left, b.right];
            })};
  })}))"""


def _long_pca():
    """The demo PCA with Step-1's level names: coloured by condition, shaped by batch."""
    data = demo_pca(counts())
    frame = data.frame.with_columns(
        cell_type=pl.col("cell_type").replace(CELLS),
        condition=pl.col("condition").replace(BATCHES),
    )
    return replace(data, frame=frame, colours=category_colours(sorted(CELLS.values())))


@pytest.fixture(scope="module")
def seen(tmp_path_factory):
    out = tmp_path_factory.mktemp("scatter_legend")
    synthetic, day = counts(), date(2026, 10, 21)
    blocks, built = [], []
    for name, data in (("pca", _long_pca()), ("volcano", demo_volcano(synthetic)),
                       ("ma", demo_ma(synthetic))):  # fmt: skip
        figure_stem = stem("de", f"{name}_legend", 1, day=day)
        reading = () if data.calls else ("Each dot is one made-up row.",)
        text = FigureText(name, reading=reading)
        write_figure(data, out, figure_stem, title=name, data_releases={}, text=text)
        blocks.append(figure_block(data, figure_stem, text))
        traces = data.interactive().data
        built.append([t.name for t in traces if t.showlegend is not False and t.name])
    page = write(_module(blocks), out / f"{stem('de', 'legend_overview', 1, day=day)}.html")
    tasks = [{"name": f"{w}_{s}", "url": page.as_uri(), "width": w, "height": 900, "scheme": s,
              "act": READ} for w in (375, 1280) for s in ("light", "dark")]  # fmt: skip
    return built, _browse(tasks, tmp_path_factory.mktemp("downloads"))


@pytest.mark.parametrize("scheme", ["light", "dark"])
def test_a_phone_keeps_every_scatter_legend_entry_inside_its_svg(seen, scheme):
    """L.2: at 375 px every entry of the PCA's, the volcano's and the MA plot's legends ends
    inside its plot's svg, and the page does not scroll sideways.

    Mutations caught: no phone legend on a scatter (entries to 439 px in a 317 px svg); the
    wrap room without the entry's end gap (an entry 1.4 px past the svg)."""
    built, results = seen
    result = results[f"375_{scheme}"]
    acted = result["acted"]
    assert result["errors"] == [] and acted["scroll"] == acted["view"], acted
    assert len(acted["figures"]) == len(built)
    for figure in acted["figures"]:
        left, right = figure["svg"]
        assert figure["items"], figure
        assert all(left - 0.5 <= a and b <= right + 0.5 for _, _, a, b in figure["items"]), figure


@pytest.mark.parametrize("scheme", ["light", "dark"])
def test_a_desktop_keeps_each_scatter_legend_as_built(seen, scheme):
    """At 1,280 px every legend entry is the figure's own name on one line, as built."""
    built, results = seen
    result = results[f"1280_{scheme}"]
    acted = result["acted"]
    assert result["errors"] == [] and acted["scroll"] == acted["view"], acted
    for names, figure in zip(built, acted["figures"], strict=True):
        assert figure["names"][: len(names)] == names, figure
        assert sorted(t for t, _, _, _ in figure["items"]) == sorted(names), figure
        assert all(lines <= 1 for _, lines, _, _ in figure["items"]), figure
