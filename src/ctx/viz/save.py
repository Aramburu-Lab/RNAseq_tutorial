"""Write every file of one figure (D16): PDF and SVG per theme, one HTML page, data and README."""

from collections.abc import Mapping
from pathlib import Path
from typing import Any, Protocol

import polars as pl

from ctx.io.naming import check
from ctx.io.readme_companion import validate, write_table
from ctx.report.builder import FigureBlock, FigurePage, write
from ctx.viz.style import THEMES, Theme


class Figure(Protocol):
    """What every figure contract offers: one data object, two renderers, one caption."""

    def table(self) -> pl.DataFrame:
        """Every number the figure prints, one row per drawn element."""

    def table_columns(self) -> Mapping[str, tuple[str, str]]:
        """(unit, meaning) for every column of ``table()``."""

    def static(self, theme: Theme) -> Any:
        """The matplotlib Figure for one theme."""

    def interactive(self) -> Any:
        """The plotly Figure."""

    def caption_text(self) -> str:
        """The caption the page prints under the figure."""


def write_figure(
    data: Figure,
    directory: Path,
    figure_stem: str,
    *,
    title: str,
    data_releases: Mapping[str, str],
) -> list[Path]:
    """Validate everything, render each theme once as PDF and SVG, then the page and the data.

    PDF and SVG come from the same Figure object, so they cannot differ in content. The page
    size is never trimmed (no ``bbox_inches="tight"``), so the width stays 85 or 170 mm. The
    CSV carries every number the figure prints.

    Args:
        data: The figure's data contract.
        directory: Output directory, created if missing.
        figure_stem: A stem from ``ctx.io.naming.stem()``.
        title: Title of the HTML page.
        data_releases: External releases the data derive from; empty when none.

    Returns:
        The written files, in ``ctx.io.naming.figure_files()`` order.

    Raises:
        ValueError: On an illegal stem or an undescribed column, before any file is written.
    """
    if not check(f"{figure_stem}.pdf"):
        raise ValueError(f"not a legal output stem: {figure_stem!r}")
    rows = data.table()
    columns = dict(data.table_columns())
    validate(rows, figure_stem, columns, ("csv",))
    directory.mkdir(parents=True, exist_ok=True)
    rendered = {}
    for theme in THEMES.values():
        fig = data.static(theme)
        for ext in ("pdf", "svg"):
            rendered[(theme.name, ext)] = directory / f"{figure_stem}_{theme.name}.{ext}"
            fig.savefig(rendered[(theme.name, ext)])
    block = FigureBlock(
        data.interactive(), data.caption_text(), rows.to_dicts(), f"{figure_stem}.csv"
    )
    page = write(FigurePage(title, block), directory / f"{figure_stem}.html")
    table = write_table(
        rows, directory, figure_stem, columns=columns, data_releases=data_releases, formats=("csv",)
    )
    return [*(rendered[(t, ext)] for ext in ("pdf", "svg") for t in THEMES), page, *table]
