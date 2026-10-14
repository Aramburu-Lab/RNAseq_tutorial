"""The one definition of an axis title and a hover value, shared by both renderings."""

import math
import re
from collections.abc import Sequence
from typing import Any

import polars as pl

from ctx.viz.annotate import wrap_title

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
