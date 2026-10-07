"""Every numeric hover value has a number format, from one rule for every contract (item 2, D9).

Each contract formatted its hover on its own or not at all: Step-1's M6 printed an FDR as
0.13016621743036838 and a rank as 1.23e+4. ``axes.hover_format`` is the rule; a missing value
is "", which plotly prints as nothing, where under a number format it printed 0 (plotly.js
4.1.1, measured in Chrome).
"""

import dataclasses
import re

import polars as pl
import pytest
from synth import (
    counts,
    demo_bar,
    demo_bars_at_scale,
    demo_box,
    demo_bubble,
    demo_curves,
    demo_density,
    demo_linked,
    demo_ma,
    demo_pca,
    demo_strip,
    demo_volcano,
)

from ctx.function import lite
from ctx.function import run as m6
from ctx.viz.axes import FLOAT_HOVER, INT_HOVER, hover_format, hover_value

FIELD = re.compile(r"%\{([^}:]+)(:[^}]*)?\}")
# The numbers each demo gains on hover: a count, and a measurement missing on the first row.
EXTRA = {"reads": ("reads", "A count."), "share": ("fraction", "A measurement.")}


def _numeric(value: object) -> bool:
    return isinstance(value, int | float) and not isinstance(value, bool)


def _with_numbers(data):
    """The demo with an integer and a float hover column more, so the customdata of every
    contract carries numbers, a missing one among them."""
    frame = (
        data.frame.with_row_index("_i")
        .with_columns(
            reads=pl.lit(12_345, pl.Int64),
            share=pl.when(pl.col("_i") == 0).then(None).otherwise(0.13016621743036838),
        )
        .drop("_i")
    )
    return dataclasses.replace(data, frame=frame, columns={**data.columns, **EXTRA},
                               hover=(*data.hover, *EXTRA))  # fmt: skip


def _running():
    positions = (10, 2_000, 12_345)
    result = lite.Enriched("SET", 3, 0.5, 2.0, 0.01, 0.0, positions,
                           tuple(0.001 * i for i in range(15_000)))  # fmt: skip
    return m6._running("contrast", result)


def _gsea():
    rows = [(f"S{i}", 40 + i, 0.5, 2.0, 0.001, 0.13016621743036838 / (i + 1)) for i in range(5)]
    return m6._gsea_bars("c", pl.DataFrame(rows, schema=list(m6.GSEA_COLUMNS), orient="row"))


def _ora():
    rows = [("up", f"S{i}", 40, 20, 4.5, 1e-6, 1e-5, "") for i in range(3)]
    return m6._ora_bars("c", pl.DataFrame(rows, schema=list(m6.ORA_COLUMNS), orient="row"), 900)


DEMOS = {
    "volcano": lambda: _with_numbers(demo_volcano(counts())),
    "ma": lambda: _with_numbers(demo_ma(counts())),
    "pca": lambda: _with_numbers(demo_pca(counts())),
    "bubble": lambda: _with_numbers(demo_bubble()),
    "bar": lambda: _with_numbers(demo_bar(counts())),
    "bar across": lambda: _with_numbers(demo_bar(counts(), horizontal=True)),
    "bars at scale": lambda: _with_numbers(demo_bars_at_scale(False)),
    "box": lambda: _with_numbers(demo_box(counts())),
    "curves": lambda: _with_numbers(demo_curves()),
    "density": lambda: _with_numbers(demo_density()),
    "linked": lambda: _with_numbers(demo_linked()),
    "strip": lambda: _with_numbers(demo_strip()),
    "running score": _running,
    "gsea": _gsea,
    "ora": _ora,
}


def test_the_rule():
    assert [hover_format(t) for t in (pl.Int32(), pl.UInt64(), pl.Float32(), pl.Float64())] == [
        INT_HOVER, INT_HOVER, FLOAT_HOVER, FLOAT_HOVER]  # fmt: skip
    assert hover_format(pl.String()) == hover_format(pl.Boolean()) == ""
    assert [hover_value(v) for v in (None, float("nan"), 0.0, 12, "g1")] == ["", "", 0.0, 12, "g1"]


@pytest.mark.parametrize("name", list(DEMOS))
def test_every_numeric_hover_has_a_format(name):
    """Every hover field of every trace: an x, y or customdata column holding a number has a
    format, and a formatted column holds no missing value plotly would print as 0. Mutations:
    a contract's customdata left unformatted; a missing value left as null."""
    fig, checked = DEMOS[name]().interactive(), 0
    for trace in fig.data:
        for field, fmt in FIELD.findall(getattr(trace, "hovertemplate", None) or ""):
            if field in ("x", "y"):
                values = list(getattr(trace, field))
            else:
                index = int(re.fullmatch(r"customdata\[(\d+)\]", field).group(1))
                values = [row[index] for row in trace.customdata]
            if any(_numeric(v) for v in values):
                checked += 1
                assert fmt, (name, trace.name, field)
                assert None not in values, (name, trace.name, field)
    assert checked >= 2, f"{name}: no numeric hover was checked"


def test_a_rank_prints_whole():
    """M6's running score: rank 12,345 printed as 1.23e+4 while ranks were floats."""
    fig = _running().interactive()
    (points,) = [t for t in fig.data if t.mode == "markers"]
    assert points.hovertemplate.startswith(f"Rank: %{{x{INT_HOVER}}}")
    assert list(points.x) == [10, 2_000, 12_345]
