"""A paired contrast's box draws one line per donor between its groups, both renderers (E.1e).

The breast cohort's tumour against normal pairs 113 donors; here five donors, one with a second
sample in one group and two with no value in the other, which no line can stand for. Each test
names the mutations it catches.
"""

from dataclasses import replace

import numpy as np
import polars as pl
import pytest

from ctx.viz.box import OFFSETS, box_data, paired_lines
from ctx.viz.style import THEMES

GENES = ("ESR1", "KRT5")
# donor: (normal, tumour) per gene; d3 has two tumour samples, d4 a missing normal value.
VALUES = {
    "d1": (1.0, 4.0),
    "d2": (2.0, 5.0),
    "d3": (3.0, 6.0),
    "d4": (None, 7.0),
    "d5": (float("nan"), 8.0),
}  # d4 and d5: no normal value, missing or not a number


def paired_frame() -> pl.DataFrame:
    """Two genes, four donors, the second gene shifted by 10; d3 sampled twice as tumour."""
    rows = []
    for g, gene in enumerate(GENES):
        for donor, (normal, tumour) in VALUES.items():
            rows.append((gene, "normal", f"{donor}n", donor,
                         None if normal is None else normal + 10 * g))  # fmt: skip
            rows.append((gene, "tumour", f"{donor}t", donor, tumour + 10 * g))
        rows.append((gene, "tumour", "d3t2", "d3", 9.0 + 10 * g))
    return pl.DataFrame(rows, schema=["feature", "group", "sample", "donor", "value"], orient="row")


@pytest.fixture
def paired():
    """The contract with ``pairs`` naming the donor column."""
    data = box_data(paired_frame(), columns={"value": ("log2 reads", "x")}, feature_label="Gene",
                    value_label="Expression", groups=("normal", "tumour"))  # fmt: skip
    return replace(data, pairs="donor")


def test_a_line_joins_each_donor_with_one_value_in_each_group(paired):
    """Mutations caught: a donor with two samples in one group joined to one of them; a donor
    with a missing value joined; lines drawn for an unpaired contrast; the groups swapped."""
    lines = paired_lines(paired)
    assert lines == {"ESR1": [("d1", 1.0, 4.0), ("d2", 2.0, 5.0)],
                     "KRT5": [("d1", 11.0, 14.0), ("d2", 12.0, 15.0)]}  # fmt: skip
    assert paired_lines(replace(paired, pairs="")) == {}
    with pytest.raises(ValueError, match="pairs column"):
        replace(paired, pairs="individual")


def test_the_files_draw_each_line_from_point_to_point_under_the_points(paired):
    """Mutations caught: the lines left out of the files; drawn from the feature's centre, not
    the boxes' offsets; drawn over the points; in ink, which the page cannot draw on both
    themes."""
    from matplotlib.collections import LineCollection  # noqa: PLC0415 - heavy

    for theme in THEMES.values():
        ax = paired.static(theme).axes[0]
        (found,) = [c for c in ax.collections if isinstance(c, LineCollection)]
        segments = [s.tolist() for s in found.get_segments()]
        assert segments == [[[-0.2, 1.0], [0.2, 4.0]], [[-0.2, 2.0], [0.2, 5.0]],
                            [[0.8, 11.0], [1.2, 14.0]], [[0.8, 12.0], [1.2, 15.0]]]  # fmt: skip
        points = [c for c in ax.collections if not isinstance(c, LineCollection)]
        assert all(found.get_zorder() < p.get_zorder() for p in points)
        assert np.allclose(found.get_colors()[0][:3], [0x8C / 255] * 3)
    assert OFFSETS == (-0.2, 0.2)


def test_the_page_draws_the_same_lines_and_boxes_where_the_files_do(paired):
    """Mutations caught: the page's lines left out; boxes grouped by plotly, whose offsets miss
    the lines' ends; the lines drawn over the points; hover naming a position, not the gene."""
    fig = paired.interactive()
    line, reference, tested = fig.data
    assert line.mode == "lines" and line.line.color == "#8C8C8C"
    x = [v for v in line.x if v is not None]
    y = [v for v in line.y if v is not None]
    assert x == [-0.2, 0.2, -0.2, 0.2, 0.8, 1.2, 0.8, 1.2]
    assert y == [1.0, 4.0, 2.0, 5.0, 11.0, 14.0, 12.0, 15.0]
    assert sorted(set(reference.x)) == [-0.2, 0.8] and sorted(set(tested.x)) == [0.2, 1.2]
    assert fig.layout.boxmode == "overlay" and reference.width == pytest.approx(0.4)
    assert list(fig.layout.xaxis.ticktext) == list(GENES)
    assert "%{customdata[" in reference.hovertemplate.split("<br>")[0]
    assert {d[-1] for d in reference.customdata} == set(GENES)
    unpaired = replace(paired, pairs="").interactive()
    assert [t.type for t in unpaired.data] == ["box", "box"]


def test_a_wrapped_name_is_drawn_and_measured_line_by_line(paired):
    """Mutations caught: a name's break left to plotly, which prints it as nothing; the least
    width measured on the whole name (genecard's OwnBox did both until E.1e)."""
    from ctx.viz.box import least_width  # noqa: PLC0415
    from ctx.viz.metrics import text_width  # noqa: PLC0415
    from ctx.viz.style import BOX_GAP_PX, SIZES, on_page  # noqa: PLC0415

    long = "Neutrophils, G-CSF\ndonors against\nsteady state"
    frame = paired.frame.with_columns(pl.col("feature").replace({"ESR1": long}))
    data = box_data(frame, columns={"value": ("log2 reads", "x")}, feature_label="Contrast",
                    value_label="Expression", groups=("normal", "tumour"))  # fmt: skip
    assert data.interactive().layout.xaxis.ticktext[0] == long.replace("\n", "<br>")
    tick = on_page(SIZES["tick"])
    ((scaled, _),) = least_width(data, (0.0, 20.0))
    assert scaled < 2 * (text_width(long, tick) + BOX_GAP_PX)
    assert scaled >= 2 * (max(text_width(line, tick) for line in long.split("\n")) + BOX_GAP_PX)
