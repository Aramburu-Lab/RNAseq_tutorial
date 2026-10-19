"""The page helpers of ``ctx.viz.heatmap``: its layouts' words and sizes, hover text and colours.

Apart from ``heatmap.py`` so that file stays under the line limit; every name here serves
``HeatmapData.interactive`` and nothing else.
"""

import itertools
from collections.abc import Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING

from ctx.viz.metrics import wrap_to
from ctx.viz.style import (
    PLOTLY_HOVER_PX,
    PLOTLY_NO_LEVEL,
    PLOTLY_TICK_PX,
    SIZES,
    on_page,
    tint,
)

if TYPE_CHECKING:
    from ctx.viz.heatmap import Swatch

# An axis drawn bare: no line, ticks or grid; labels clear of the cells; margins grown to hold
# them.
BARE = {"showline": False, "ticks": "", "ticklabelstandoff": PLOTLY_TICK_PX, "zeroline": False,
         "showgrid": False, "automargin": True}  # fmt: skip


@dataclass(frozen=True)
class Words:
    """What one page layout prints, each text with ``<br>`` where it breaks.

    Row labels by grid row, legend entries and the column title.
    """

    ticks: Mapping[int, str]
    legend: list[str]
    title: str


@dataclass(frozen=True)
class Wide:
    """The static's layout on the page (``HeatmapData._wide``).

    Attributes:
        width: The plot's width it is laid out for, in px: a module page's, or the least its
            labels need where that is more.
        band: A row's height there, in px: the unit of ``edges``.
        ratio: A row's height over a column's width there.
        below: The narrowest plot it is drawn in, in px.
        held: Whether plotly holds ``ratio`` at every width it is drawn at.
        edges: The bands' edges top down, in rows.
    """

    width: float
    band: float
    ratio: float
    below: int
    held: bool
    edges: list[float]


def page_colour(swatch: "Swatch") -> str:
    """A swatch's colour on the page; "nothing to show" is a grey that reads in both themes."""
    colour, alpha = (swatch.colour, 1.0) if swatch.colour else PLOTLY_NO_LEVEL
    return tint(colour, alpha * swatch.alpha)


def band_centres(edges: list[float]) -> list[float]:
    """The middle of each band between consecutive ``edges``."""
    return [(top + bottom) / 2 for top, bottom in itertools.pairwise(edges)]


def named(label: str, value: str) -> str:
    """``label: value``; a label written as a sentence loses its full stop first."""
    return f"{label.removesuffix('.')}: {value}"


def hover_lines(lines: list[str]) -> str:
    """Hover lines, each wrapped to a phone's plot, which a wider hover box ran off."""
    return "<br>".join(fit_px(line, PLOTLY_HOVER_PX) for line in lines)


def fit_px(text: str | list[str], width_px: float, size: float = SIZES["tick"]) -> str:
    """Text wrapped at spaces to lines of at most ``width_px`` on the page, in plotly's ``<br>``.

    Measured as ``ctx.viz.metrics`` measures every wrapped text. A word longer than a line
    stays whole: split, it would read as two words. A list is the text's words already, as
    ``metrics.words`` holds a parenthesised unit whole.
    """
    lines = wrap_to(text, width_px, on_page(size))
    return "<br>".join(lines) if lines else (text if isinstance(text, str) else "")


def line_count(text: str) -> int:
    """How many lines a plotly text runs to."""
    return text.count("<br>") + 1
