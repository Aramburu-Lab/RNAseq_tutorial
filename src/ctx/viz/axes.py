"""The one definition of an axis title, shared by the static and the interactive figures."""

import re

# Units that are no unit: the title is the label alone.
NO_UNIT = frozenset({"none"})


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
    if unit in NO_UNIT or _words(unit) <= _words(label):
        return label
    return f"{label} ({unit})"
