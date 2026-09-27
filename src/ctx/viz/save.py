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

    def page_table(self) -> pl.DataFrame:
        """The rows the HTML page draws and embeds; the CSV beside it holds ``table()``."""

    def static(self, theme: Theme) -> Any:
        """The matplotlib Figure for one theme."""

    def interactive(self) -> Any:
        """The plotly Figure."""

    def caption_text(self) -> str:
        """The caption the page prints under the figure."""


def drawn_texts(fig: Any) -> list[Any]:
    """Visible, non-empty texts, without the labels of ticks outside the view limits.

    matplotlib keeps a label for every tick its locator proposes but draws only those inside
    the view interval; an undrawn label sits wherever its tick would be.
    """
    from matplotlib.text import Text  # noqa: PLC0415 - heavy; imported after the profile

    hidden = set()
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            low, high = sorted(axis.get_view_interval())
            span = high - low
            for tick, loc in zip(axis.get_major_ticks(), axis.get_majorticklocs(), strict=False):
                if not low - 1e-9 * span <= loc <= high + 1e-9 * span:
                    hidden |= {id(tick.label1), id(tick.label2)}
    return [
        t
        for t in fig.findobj(Text)
        if t.get_visible() and t.get_text().strip() and id(t) not in hidden
    ]


def layout_problems(fig: Any) -> list[str]:
    """Texts that leave the page or overlap another text, measured on the drawn figure."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg  # noqa: PLC0415 - heavy

    renderer = FigureCanvasAgg(fig).get_renderer()  # type: ignore[no-untyped-call]
    fig.draw(renderer)
    page = fig.bbox
    boxes = [(t.get_text(), t.get_window_extent(renderer)) for t in drawn_texts(fig)]
    problems = [
        f"off the page: {name!r}"
        for name, b in boxes
        if b.x0 < page.x0 - 1 or b.x1 > page.x1 + 1 or b.y0 < page.y0 - 1 or b.y1 > page.y1 + 1
    ]
    problems += [
        f"overlap: {a!r} and {b!r}"
        for i, (a, box_a) in enumerate(boxes)
        for b, box_b in boxes[i + 1 :]
        if box_a.overlaps(box_b)
    ]
    return problems


def figure_block(data: Figure, figure_stem: str) -> FigureBlock:
    """The interactive figure, caption and data rows, for a figure page or a module page."""
    return FigureBlock(
        data.interactive(), data.caption_text(), data.page_table().to_dicts(), f"{figure_stem}.csv"
    )


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
        ValueError: On an illegal stem or an undescribed column, before any file is written;
            or if a text leaves the page or overlaps another.
    """
    if not check(f"{figure_stem}.pdf"):
        raise ValueError(f"not a legal output stem: {figure_stem!r}")
    rows = data.table()
    columns = dict(data.table_columns())
    validate(rows, figure_stem, columns, ("csv",))
    figures = {theme.name: data.static(theme) for theme in THEMES.values()}
    for name, fig in figures.items():
        if problems := layout_problems(fig):
            raise ValueError(f"{figure_stem} ({name}): " + "; ".join(problems[:5]))
    directory.mkdir(parents=True, exist_ok=True)
    rendered = {}
    for name, fig in figures.items():
        for ext in ("pdf", "svg"):
            rendered[(name, ext)] = directory / f"{figure_stem}_{name}.{ext}"
            fig.savefig(rendered[(name, ext)])
    block = figure_block(data, figure_stem)
    page = write(FigurePage(title, block), directory / f"{figure_stem}.html")
    table = write_table(
        rows, directory, figure_stem, columns=columns, data_releases=data_releases, formats=("csv",)
    )
    return [*(rendered[(t, ext)] for ext in ("pdf", "svg") for t in THEMES), page, *table]
