"""M4's page budget: a sample-by-sample heatmap thins whole samples, a scatter keeps its lesson.

Mutations each test was seen to fail on are named in its docstring (2026-10-08).
"""

import numpy as np
import polars as pl

from ctx.quant import figures as fig
from ctx.viz.style import HEAT

ROWS = HEAT["max_page_cells"]


def sheet(n: int) -> pl.DataFrame:
    """n samples of two cell types and two batches."""
    return pl.DataFrame({"sample": [f"s{i:03d}" for i in range(n)],
                         "cell_type": ["a", "b"] * (n // 2), "batch": ["x"] * n})  # fmt: skip


def heatmap(n: int):
    """F3's contract over a made-up rho of n samples."""
    rng = np.random.default_rng(1)
    noise = rng.normal(0, 0.05, (n, n))
    rho = np.clip(0.8 + (noise + noise.T) / 2, -1, 1)
    np.fill_diagonal(rho, 1.0)
    return fig.correlation_heatmap(rho, sheet(n), "caption")


def test_heatmap_thins_whole_samples():
    """Above the cell budget the page draws every k-th sample, its row and column whole.

    Mutation: the page keeps a thinned sample's row but not its column (the column filter
    dropped in ``_paged``): cells name samples the page does not draw.
    """
    data = heatmap(120)
    assert data.frame.height == 120 * 119 // 2 > ROWS
    assert data.stride() == 2
    page = data.page_table()
    kept = set(data.rows_in_order()[::2])
    assert set(page["sample_row"]) | set(page["sample_column"]) <= kept
    assert page.height == 60 * 59 // 2 <= ROWS
    assert "draws 60 of the 120 rows" in data.caption_text()
    assert data.caption == "caption"


def test_stride_counts_cells():
    """The stride counts the frame's cells, so a lower triangle under the budget is whole.

    Mutation: the stride counting units squared, as a full matrix would: 100 samples would
    then be thinned although their 4,950 cells fit.
    """
    data = heatmap(100)
    assert data.frame.height == 4_950 <= ROWS
    assert data.stride() == 1
    assert data.caption_text() == "caption"
    assert heatmap(200).stride() == 2  # 100 samples, 4,950 cells; units squared would say 3
