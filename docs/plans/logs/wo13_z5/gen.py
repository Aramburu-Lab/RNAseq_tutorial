"""The scatter-legend test's module page, written by whichever src is on the path."""

import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, "tests")
from synth import counts, demo_ma, demo_volcano  # noqa: E402
from test_export_browser import _module  # noqa: E402
from test_scatter_legend_browser import _long_pca  # noqa: E402

from ctx.io.naming import stem  # noqa: E402
from ctx.report.builder import FigureText, write  # noqa: E402
from ctx.viz.save import figure_block, write_figure  # noqa: E402

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
synthetic, day = counts(), date(2026, 10, 21)
blocks = []
for name, data in (
    ("pca", _long_pca()),
    ("volcano", demo_volcano(synthetic)),
    ("ma", demo_ma(synthetic)),
):
    figure_stem = stem("de", f"{name}_legend", 1, day=day)
    reading = () if data.calls else ("Each dot is one made-up row.",)
    text = FigureText(name, reading=reading)
    write_figure(data, out, figure_stem, title=name, data_releases={}, text=text)
    blocks.append(figure_block(data, figure_stem, text))
print(write(_module(blocks), out / f"{stem('de', 'legend_overview', 1, day=day)}.html"))
