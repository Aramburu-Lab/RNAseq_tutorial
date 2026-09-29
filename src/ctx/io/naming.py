"""Output names: ``YYYY-MM-DD_<module>_<description>_v<N>``, built and checked in one place.

A figure appends one theme token after the version (``_white`` or ``_black``), and a data
file's companion appends ``_README.md``. Nothing else may follow the version.
"""

import re
from datetime import date
from pathlib import Path

# (token, NN): the NN numbers the per-module output directory, outputs/NN_<token>/.
MODULES: tuple[tuple[str, int], ...] = (
    ("setup", 0),
    ("gather", 1),
    ("qc", 2),
    ("process", 3),
    ("quant", 4),
    ("de", 5),
    ("function", 6),
    ("proteo", 7),
    ("ncrna", 8),
    ("report", 9),
    ("lite_build", 10),
)
THEMES = ("white", "black")
# Theme words may only appear as the trailing token, or `volcano_white_v1` would pass.
_RESERVED = {*THEMES, "light", "dark", "readme"}
_TOKENS = "|".join(sorted((t for t, _ in MODULES), key=len, reverse=True))
_NAME = re.compile(
    rf"(?P<day>\d{{4}}-\d{{2}}-\d{{2}})_(?P<module>{_TOKENS})_"
    r"(?P<desc>[a-z0-9]+(?:_[a-z0-9]+)*)_v(?P<version>[1-9]\d*)"
    rf"(?:_(?P<theme>{'|'.join(THEMES)}))?"
    r"(?:(?P<readme>_README\.md)|\.(?P<ext>[a-z0-9]+))"
)


def check(name: str) -> bool:
    """Say whether a file name follows the output naming rule.

    Args:
        name: A bare file name, without directories.

    Returns:
        True when the name is legal.
    """
    match = _NAME.fullmatch(name)
    if match is None or (match["readme"] and match["theme"]):
        return False
    if _RESERVED & set(match["desc"].split("_")):
        return False
    try:
        date.fromisoformat(match["day"])
    except ValueError:
        return False
    return True


def stem(module: str, description: str, version: int, *, day: date) -> str:
    """Build the stem every output of one step shares.

    Args:
        module: A token from ``MODULES``, e.g. ``"de"``.
        description: Lowercase words joined by underscores, e.g. ``"volcano"``.
        version: Positive integer, bumped when the content changes meaning.
        day: The date the output was made.

    Returns:
        The stem, e.g. ``"2026-09-26_de_volcano_v1"``.

    Raises:
        ValueError: If the parts do not form a legal name.
    """
    built = f"{day.isoformat()}_{module}_{description}_v{version}"
    if not check(f"{built}.csv"):
        raise ValueError(f"not a legal output stem: {built!r}")
    return built


def figure_files(figure_stem: str) -> list[str]:
    """Every file one figure writes (D16): PDF and SVG per theme, one HTML, the data and README.

    Args:
        figure_stem: A stem from ``stem()``.

    Returns:
        The file names, in a fixed order.
    """
    rendered = [f"{figure_stem}_{theme}.{ext}" for ext in ("pdf", "svg") for theme in THEMES]
    return [*rendered, f"{figure_stem}.html", f"{figure_stem}.csv", f"{figure_stem}_README.md"]


# Tool-native output folders keep the tool's own file names, as nf-core's `results/` does
# (AGENTS.md): MultiQC's report and tables under `outputs/02_qc/multiqc/` (D-028).
NATIVE_DIRS = frozenset({"multiqc"})


def lint(directory: Path) -> tuple[int, list[str]]:
    """Check every non-hidden file under a directory, outside tool-native folders.

    Args:
        directory: The output directory to scan, recursively.

    Returns:
        The number of files scanned, and the relative paths whose names break the rule.
    """
    files = [
        p
        for p in sorted(directory.rglob("*"))
        if p.is_file()
        and not any(part.startswith(".") for part in p.relative_to(directory).parts)
        and not NATIVE_DIRS & set(p.relative_to(directory).parts[:-1])
    ]
    return len(files), [str(p.relative_to(directory)) for p in files if not check(p.name)]
