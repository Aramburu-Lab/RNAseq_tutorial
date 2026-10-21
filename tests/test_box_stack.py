"""V.3: a box figure stacks its p and r labels and its names in rows once one row overlaps.

Measured 2026-10-21 at 170 mm with ``p=3.7e-05`` labels: 10 features passed, 11 gave 10
overlaps and ``write_figure`` refused, though a dataset file may list any number of box genes.
Names overlapping in one row were thinned instead, every other one hidden: ten 12-character
names lost five. matplotlib's ``get_xticklabels`` returns only the visible names, so counting
them counts what thinning left.

V.4: the size a page's export of a view takes, its rows, and the legend's lift.
"""

import dataclasses
import math
from datetime import date

import polars as pl
import pytest

from ctx.io.naming import stem
from ctx.viz.box import Stack, box_data, stacking
from ctx.viz.save import layout_problems, thin_tick_labels, write_figure
from ctx.viz.style import STACK_GAP_PT, THEMES


def boxes(n: int, p: float = 3.7e-05, name: str = "GENE{:02d}"):
    """``n`` features of eight samples a group, each label printing ``p``."""
    rows = [(name.format(i), g, f"{g}{j}", float(j + (20 if g == "case" else 0) + i))
            for i in range(n) for g in ("control", "case") for j in range(8)]  # fmt: skip
    frame = pl.DataFrame(rows, schema=["feature", "group", "sample", "value"], orient="row")
    columns = {"feature": ("none", "Gene."), "group": ("none", "Group."),
               "sample": ("none", "Sample."), "value": ("log2 CPM", "Expression.")}  # fmt: skip
    data = box_data(frame, columns=columns, feature_label="Gene",
                    value_label="Expression", groups=("control", "case"))  # fmt: skip
    tests = tuple(dataclasses.replace(t, p_adj=p) for t in data.caption.tests)
    return dataclasses.replace(data, caption=dataclasses.replace(data.caption, tests=tests))


CROWDED = {"11": (11, "GENE{:02d}"), "15": (15, "GENE{:02d}"), "20": (20, "GENE{:02d}"),
           "10_long_names": (10, "LOC1053728{:02d}")}  # fmt: skip


@pytest.mark.parametrize("case", list(CROWDED))
def test_every_label_and_name_is_drawn_apart_when_one_row_cannot_hold_them(case, tmp_path):
    """Every p and r label and every name drawn, none overlapping, none thinned away, every
    label inside the axes, and the file written. Mutations caught: the label rows held at one
    (``labels=1`` in ``stacking``), which leaves 10 overlaps at 11 features; the name rows held
    at one, which thins ten long names to five."""
    n, name = CROWDED[case]
    for theme in THEMES.values():
        fig = boxes(n, name=name).static(theme)
        thin_tick_labels(fig)
        ax = fig.axes[0]
        assert len(ax.get_xticklabels()) == n
        assert layout_problems(fig) == []
        frame = ax.get_window_extent()
        notes = [t for t in ax.texts if t.get_text().startswith("p=")]
        assert len(notes) == n
        assert all(frame.x0 <= t.get_window_extent().x0 and t.get_window_extent().x1 <= frame.x1
                   for t in notes)  # fmt: skip
    write_figure(boxes(n, name=name), tmp_path, stem("de", f"box_{n}", 1, day=date(2026, 10, 21)),
                 title="t", data_releases={})  # fmt: skip


def test_ten_features_are_drawn_in_one_row_as_before():
    """At ten features one row holds every label: no stacking, so the figure is the one drawn
    before stacking existed (its PDF and SVG bytes were compared on fixtures, 2026-10-21)."""
    fig = boxes(10).static(THEMES["white"])
    ax = fig.axes[0]
    assert stacking(fig, ax, [t for t in ax.texts if t.get_text().startswith("p=")], Stack()) \
        == Stack()  # fmt: skip


def _top_label_and_legend(fig) -> tuple[float, float, float]:
    """(highest label's top, legend's bottom, axes' top) in pt, on the drawn figure."""
    fig.canvas.draw()
    ax, to_pt = fig.axes[0], 72 / fig.dpi
    notes = [t.get_window_extent() for t in ax.texts if t.get_text().startswith("p=")]
    return (max(b.y1 for b in notes) * to_pt, ax.get_legend().get_window_extent().y0 * to_pt,
            ax.get_window_extent().y1 * to_pt)  # fmt: skip


STACKED = {**CROWDED, "30": (30, "GENE{:02d}")}
# (label rows, legend lifted) at p=3.7e-05, measured 2026-10-21.
ROWS = {"11": (2, False), "15": (2, False), "20": (3, True), "30": (5, True),
        "10_long_names": (1, False)}  # fmt: skip


@pytest.mark.parametrize("case", ["4", *STACKED])
def test_the_size_the_export_takes_is_the_size_the_files_draw(case):
    """V.4: ``inches()``, the size a page's PDF and SVG of a view take, is the drawn figure's.
    Mutation caught: ``inches()`` left at the one-row size (3.68 in against 3.95 at 20)."""
    n, name = STACKED.get(case, (4, "GENE{:02d}"))
    data = boxes(n, name=name)
    assert data.inches() == tuple(data.static(THEMES["white"]).get_size_inches())


@pytest.mark.parametrize("case", list(STACKED))
def test_the_legend_rises_only_as_far_as_the_highest_label_needs(case):
    """V.4: the legend clears the highest label by ``STACK_GAP_PT``, or stays on the axes where
    the headroom holds every row, and the axes keep their one-row height. At 11 features the
    second row sits inside the axes, so the page keeps its one-row 3.68 in (V.3 drew 3.97).
    Mutation caught: the legend lifted by every extra row (11 and 15 lifted 20.6 pt for none)."""
    n, name = STACKED[case]
    data = boxes(n, name=name)
    top, legend, axes = _top_label_and_legend(data.static(THEMES["white"]))
    if data.stack.legend_pt:
        assert legend - top == pytest.approx(STACK_GAP_PT, abs=0.02)
    else:
        assert legend == pytest.approx(axes, abs=0.01) and top <= axes - STACK_GAP_PT
    plain, wide = boxes(4).static(THEMES["white"]), data.static(THEMES["white"])
    assert wide.axes[0].get_window_extent().height == pytest.approx(
        plain.axes[0].get_window_extent().height, abs=0.01)  # fmt: skip
    assert (data.stack.labels, data.stack.legend_pt > 0) == ROWS[case]


@pytest.mark.parametrize("case", list(STACKED))
def test_the_export_takes_the_rows_of_the_files(case):
    """V.4: the page's export puts label i in row ``i % labels`` of the files' ``Stack``, every
    label shown, the names upright. Mutation caught: the export left with nothing, which draws
    the page's phone layout, its labels hidden."""
    n, name = STACKED[case]
    data = boxes(n, name=name)
    layers, stack = data.export_layers(), data.stack
    relayout = layers["relayout"]
    pitch = relayout["annotations[1].yshift"] if stack.labels > 1 else 1.0
    order = {t.feature: j for j, t in enumerate(data.caption.tests)}
    features = data.frame["feature"].unique(maintain_order=True).to_list()
    for i, f in enumerate(features):
        at = f"annotations[{order[f]}]"
        assert relayout[f"{at}.visible"] is True
        assert relayout[f"{at}.yshift"] == pytest.approx(i % stack.labels * pitch)
    assert (stack.labels == 1 or pitch >= stack.label_pt) and relayout["xaxis.tickangle"] == 0
    assert ("lift" in layers) == (stack.labels > 1)


def test_a_missing_value_leaves_the_labels_over_the_boxes():
    """V.4: numpy's max over a column holding a NaN is NaN, and the labels went to y = NaN, never
    drawn. Mutation caught: the top taken over every value, missing ones included."""
    frame = boxes(3).frame.with_columns(
        value=pl.when(pl.int_range(pl.len()) == 0).then(float("nan")).otherwise(pl.col("value"))
    )
    data = dataclasses.replace(boxes(3), frame=frame)
    ax = data.static(THEMES["white"]).axes[0]
    ys = [t.get_position()[1] for t in ax.texts if t.get_text().startswith("p=")]
    assert len(ys) == 3 and all(math.isfinite(y) for y in ys)  # noqa: PLR2004 - three genes
    assert all(math.isfinite(a.y) for a in data.interactive().layout.annotations)
