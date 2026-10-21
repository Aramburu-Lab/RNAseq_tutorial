"""The one definition of an axis title, a hover value and a phone's legend, shared by figures."""

import math
import re
from collections.abc import Sequence
from typing import Any

import polars as pl

from ctx.viz.annotate import wrap_title
from ctx.viz.metrics import text_width, wrap_to
from ctx.viz.style import (
    PHONE_PLOT_PX,
    PLOTLY_LEGEND_ENTRY_PX,
    PLOTLY_MARGIN_PX,
    SIZES,
    on_page,
)

# Units that are no unit: the title is the label alone.
NO_UNIT = frozenset({"none"})
# Plotly never wraps an axis title; an x title longer than this runs past a 375 px phone's
# figure (measured on M3's scatters: 29 to 31 characters crossed the edge by 5 to 24 px), so
# its unit takes a line of its own on the page. A y title takes the figure files' lines.
SCREEN_TITLE_CHARS = 26
# Hover number formats: significant digits for measurements, thousands for counts.
FLOAT_HOVER, INT_HOVER = ":.3~g", ":,d"


def _words(text: str) -> set[str]:
    """Lower-case words without a plural s, so "genes" and "Genes called" share "gene"."""
    return {w.removesuffix("s") for w in re.findall(r"[a-z0-9]+", text.lower())}


def axis_title(label: str, unit: str) -> str:
    """``label (unit)``, or the label alone when it already names the unit.

    A count such as "Genes (genes)" repeats itself: the label still names what is counted,
    so every axis keeps its unit in words.

    Args:
        label: What the axis shows.
        unit: The unit from the column dictionary.

    Returns:
        The title the figure prints.
    """
    words = _words(unit)
    # A unit with no word in it, such as %, is never already named by the label.
    if unit in NO_UNIT or (words and words <= _words(label)):
        return label
    return f"{label} ({unit})"


def screen_title(title: str) -> str:
    """An interactive x-axis title: ``axis_title``'s, its unit on a second line when it is long.

    Args:
        title: The title ``axis_title`` gave.

    Returns:
        The title with a line break before the unit's opening parenthesis, or unchanged when
        it is short or has no unit.
    """
    label, bracket, unit = title.rpartition(" (")
    if len(title) <= SCREEN_TITLE_CHARS or not bracket or not label:
        return title
    return f"{label}<br>({unit}"


def page_y_title(title: str, room: float) -> str:
    """A y-axis title as the page prints it: the lines the figure files print, from one call.

    The page's own PDF and SVG of a moved view copy the page's y title, so a page that broke
    it elsewhere would carry that break into them too.

    Args:
        title: The title ``axis_title`` gave.
        room: The files' drawn axis height in pt, ``annotate.files_room`` (P2S.8).

    Returns:
        ``wrap_title``'s lines, joined by plotly's line break.
    """
    return "<br>".join(wrap_title(title, room).split("\n"))


def hover_format(dtype: pl.DataType) -> str:
    """Plotly's number format for a hover value of this type: none for text.

    A whole number prints with thousands separators (12,345), a measurement to three
    significant digits; unformatted, plotly printed an FDR as 0.13016621743036838.

    Args:
        dtype: The column's polars type.

    Returns:
        A ``:``-led d3 format, or "" for a value printed as it is.
    """
    if dtype.is_integer():
        return INT_HOVER
    return FLOAT_HOVER if dtype.is_float() else ""


def hover_value(value: Any) -> Any:
    """A value as plotly's customdata carries it: a missing one (None or NaN) as "".

    Under a number format plotly.js 4.1.1 prints a missing value (null) as 0 (measured in
    Chrome), and "" as nothing, which is what is there.

    Args:
        value: One hover value.

    Returns:
        The value, or "" for a missing one.
    """
    return "" if value is None or (isinstance(value, float) and math.isnan(value)) else value


def hover_rows(frame: pl.DataFrame, columns: Sequence[str]) -> list[tuple[Any, ...]]:
    """Each row's hover values for plotly's customdata, through ``hover_value``.

    Args:
        frame: The rows a trace draws.
        columns: The hover columns, in customdata order.

    Returns:
        One tuple per row.
    """
    return [tuple(hover_value(v) for v in row) for row in frame.select(columns).rows()]


def hover_field(field: str, dtype: pl.DataType) -> str:
    """A plotly hover template field, ``%{field}`` with the format ``hover_format`` gives.

    Args:
        field: The plotly name, such as ``x`` or ``customdata[2]``.
        dtype: The type of the values it shows.

    Returns:
        The field as a hover template prints it.
    """
    return f"%{{{field}{hover_format(dtype)}}}"


def phone_legend(fig: Any) -> dict[str, Any] | None:
    """A phone's legend, its entries wrapped to a phone's plot, where an entry is wider.

    At 375 px M1's depth legend ran 38 px past its plot and off the screen: "ENCODE, older
    projects: 20 million" is 355 px with its key (work-order item 2's leftover), and M5's PCA
    legend reached 407 px in a 317 px svg on Step-1 (L.2), so bars and scatters both take
    it. The page script swaps the wrapped names in below ``below`` px (``_page.js.j2``,
    fit): traces by restyle, the lines' shapes by relayout. An entry takes its symbol and a
    gap beside its text (``PLOTLY_LEGEND_ENTRY_PX``).

    Returns:
        The narrow layout, or None when every entry fits a phone's plot.
    """
    size = on_page(SIZES["legend"])
    traces = [(i, t.name) for i, t in enumerate(fig.data) if t.name and t.showlegend is not False]
    shapes = [(i, s.name) for i, s in enumerate(fig.layout.shapes) if s.name and s.showlegend]
    widest = max((text_width(n, size) for _, n in traces + shapes), default=0.0)
    room = PHONE_PLOT_PX - PLOTLY_LEGEND_ENTRY_PX
    if widest <= room:
        return None
    wrapped = {n: "<br>".join(wrap_to(n, room, size)) for _, n in traces + shapes}
    return {
        "below": widest + PLOTLY_LEGEND_ENTRY_PX + PLOTLY_MARGIN_PX["l"] + PLOTLY_MARGIN_PX["r"],
        "relayout": {f"shapes[{i}].name": wrapped[n] for i, n in shapes},
        "restyle": [[{"name": [wrapped[n] for _, n in traces]}, [i for i, _ in traces]]]
        if traces else [],
    }  # fmt: skip
