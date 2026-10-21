"""V.3: a box figure stacks its p and r labels and its names in rows once one row overlaps.

Measured 2026-10-21 at 170 mm with ``p=3.7e-05`` labels: 10 features passed, 11 gave 10
overlaps and ``write_figure`` refused, though a dataset file may list any number of box genes.
Names overlapping in one row were thinned instead, every other one hidden: ten 12-character
names lost five. matplotlib's ``get_xticklabels`` returns only the visible names, so counting
them counts what thinning left.
"""

import dataclasses
from datetime import date

import polars as pl
import pytest

from ctx.io.naming import stem
from ctx.viz.box import Stack, box_data, stacking
from ctx.viz.save import layout_problems, thin_tick_labels, write_figure
from ctx.viz.style import THEMES


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
