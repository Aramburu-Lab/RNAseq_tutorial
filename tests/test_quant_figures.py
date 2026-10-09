"""M4's page budget: a sample-by-sample heatmap thins whole samples, a scatter keeps its lesson.

Mutations each test was seen to fail on are named in its docstring (2026-10-08).
"""

import numpy as np
import polars as pl
import pytest

from ctx.quant import figures as fig
from ctx.viz.style import HEAT, SPENT_GREYS, contrast, strip_colours

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


def test_cells_are_square():
    """A sample-by-sample heatmap draws square cells in the files and on the page, thinned too.

    Mutations: the static's equal aspect dropped (``ax.set_aspect`` in ``static``): the 120
    samples' cells stand taller than wide; ``_paged`` resetting ``thin``: the page's 60 drawn
    samples fall back to 60 capped rows of 4.5 mm and stand taller than wide, as Step-1's
    105 did, 2.9 times (2026-10-09); the track bands not lifted in ``_wide``: plotly holds
    no aspect and the page draws cells wider than tall past 686 px.
    """
    data = heatmap(120)
    paged = data._paged(data.stride())
    assert abs(data.cell_aspect() - 1) < 0.01
    assert abs(paged.cell_aspect() - 1) < 0.01
    # 105 drawn, as Step-1's page: rows unlabelled, so plotly must hold the aspect itself.
    page = heatmap(210)._paged(2).interactive().layout.yaxis
    assert page.scaleanchor == "x"
    assert abs(page.scaleratio - 1) < 0.01


def test_track_names_keep_their_row():
    """A square cell is shorter than a label, so a track's row keeps the 4.5 mm a label needs.

    Mutation: every static row one unit tall (``_static_edges`` without the thin case): 120
    columns give 1.4 mm rows and the two track names overprint.
    """
    data = heatmap(120)
    edges = data._static_edges()
    assert (edges[1] - edges[0]) * data._pitch_mm() >= HEAT["row_mm"]
    assert edges[-1] - edges[3] == 120  # two tracks and the gap, then one unit per sample


def test_seven_cell_types_take_seven_colours_and_the_batches_none_of_them():
    """Step-1's 7 cell types each take a colour, and the batch strip repeats none of them.

    Mutation: each strip coloured from the first data colour (``category_colours`` per strip,
    as before 2026-10-09): the two batches take the first two cell types' colours.
    """
    n = 14
    frame = pl.DataFrame({"sample": [f"s{i:02d}" for i in range(n)],
                          "cell_type": [f"c{i % 7}" for i in range(n)],
                          "batch": ["short", "long"] * (n // 2)})  # fmt: skip
    rng = np.random.default_rng(2)
    rho = np.clip(0.8 + rng.normal(0, 0.05, (n, n)), -1, 1)
    rho = (rho + rho.T) / 2
    np.fill_diagonal(rho, 1.0)
    cells, batches = fig.correlation_heatmap(rho, frame, "caption").tracks
    assert len(set(cells.colours.values())) == 7
    assert not set(cells.colours.values()) & set(batches.colours.values())
    assert len(set(batches.colours.values())) == 2


def test_strip_colours_refuse_what_greys_cannot_tell_apart():
    """A strip of three levels with no data colour left is refused, never given repeats.

    Mutation: the spent-colours branch taking any strip (its length check dropped): three
    levels zipped onto two greys, one level left without a colour.
    """
    spent = strip_colours([list("abcdefg"), ["x", "y"]])[1]
    assert spent == dict(zip("xy", SPENT_GREYS, strict=True))
    with pytest.raises(ValueError, match="3 strip levels; 0 data colours left"):
        strip_colours([list("abcdefg"), ["x", "y", "z"]])
    assert contrast(*SPENT_GREYS) >= 3  # the 3:1 graphics threshold, between the two levels
