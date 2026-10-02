"""The one definition of an axis title, shared by the static and the interactive figures."""

import re

# Units that are no unit: the title is the label alone.
NO_UNIT = frozenset({"none"})
# Plotly never wraps an axis title; one longer than this runs past a 375 px phone's figure
# (measured on M3's scatters: 29 to 31 characters crossed the edge by 5 to 24 px), so its unit
# takes a line of its own on the page. The PDF and SVG wrap by measured width instead.
SCREEN_TITLE_CHARS = 26


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
    """An interactive axis title: ``axis_title``'s, its unit on a second line when it is long.

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
