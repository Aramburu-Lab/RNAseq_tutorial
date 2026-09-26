"""Write every file of one figure (D16): PDF and SVG per theme, the data and its README."""

from collections.abc import Mapping
from pathlib import Path

from ctx.io.naming import check
from ctx.io.readme_companion import validate, write_table
from ctx.viz import box
from ctx.viz.style import THEMES


def write_figure(
    data: box.BoxData, directory: Path, figure_stem: str, *, data_releases: Mapping[str, str]
) -> list[Path]:
    """Validate everything, render each theme once and save it as PDF and SVG, then the data.

    PDF and SVG come from the same Figure object, so they cannot differ in content. The page
    size is never trimmed (no ``bbox_inches="tight"``), so the width stays 85 or 170 mm.

    Args:
        data: The figure's data contract.
        directory: Output directory, created if missing.
        figure_stem: A stem from ``ctx.io.naming.stem()``.
        data_releases: External releases the data derive from; empty when none.

    Returns:
        The written files, in ``ctx.io.naming.figure_files()`` order minus the HTML.

    Raises:
        ValueError: On an illegal stem or an undescribed column, before any file is written.
    """
    if not check(f"{figure_stem}.pdf"):
        raise ValueError(f"not a legal output stem: {figure_stem!r}")
    validate(data.frame, figure_stem, data.columns, ("csv",))
    directory.mkdir(parents=True, exist_ok=True)
    rendered = {}
    for theme in THEMES.values():
        fig = box.static(data, theme)
        for ext in ("pdf", "svg"):
            rendered[(theme.name, ext)] = directory / f"{figure_stem}_{theme.name}.{ext}"
            fig.savefig(rendered[(theme.name, ext)])
    ordered = [rendered[(t, ext)] for ext in ("pdf", "svg") for t in THEMES]
    table = write_table(
        data.frame, directory, figure_stem, columns=data.columns, data_releases=data_releases,
        formats=("csv",),
    )  # fmt: skip
    return [*ordered, *table]
