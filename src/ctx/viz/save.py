"""Write every file of one figure (D16): PDF and SVG per theme, one HTML page, data and README."""

from collections.abc import Iterable, Iterator, Mapping
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import replace
from pathlib import Path
from typing import Any, Protocol

import polars as pl

from ctx.figures import panel, record
from ctx.io.naming import check
from ctx.io.readme_companion import validate, write_table
from ctx.report import export, glossary
from ctx.report.builder import FigureBlock, FigurePage, FigureText, write
from ctx.report.figures import finite
from ctx.report.marks import MARK
from ctx.viz.bar import BarData
from ctx.viz.contour import SurfaceData
from ctx.viz.curve import CurveData
from ctx.viz.km import KmData
from ctx.viz.linked import LinkedData
from ctx.viz.scatter import ScatterData
from ctx.viz.style import THEMES, Theme

# A figure with no title, lead or reading guide.
NO_TEXT = FigureText()
# The contracts that can carry page controls and draw their changed view (D-033, D-051).
CONTROLLED = (ScatterData, LinkedData, BarData, CurveData, SurfaceData, KmData)
# The Reproduce panel of each figure this build wrote with a spec, by stem (Phase 2D). A stem
# names no root, so a second build in one process would find the first's panels before its
# own files exist; ``one_build`` gives each build its own (Z.4).
_PANELS: ContextVar[dict[str, dict[str, str]] | None] = ContextVar("panels", default=None)
# The panels of figures written outside any ``one_build``: one CLI command is one process.
_PROCESS: dict[str, dict[str, str]] = {}


def _panels() -> dict[str, dict[str, str]]:
    """The panels of the build in progress."""
    found = _PANELS.get()
    return _PROCESS if found is None else found


@contextmanager
def one_build() -> Iterator[None]:
    """Hold the panels written inside to this build: none before, none kept after."""
    token = _PANELS.set({})
    try:
        yield
    finally:
        _PANELS.reset(token)


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
    """Visible, non-empty texts, without the tick labels matplotlib does not draw.

    matplotlib keeps a label for every tick its locator proposes but draws only those inside
    the view interval; an undrawn label sits wherever its tick would be. An axis switched off
    draws none of its labels, though each still reports itself visible: on Step-1's 210
    samples the parallel sets' undrawn 0 to 1 labels sat past the page (2026-09-30).
    """
    from matplotlib.text import Text  # noqa: PLC0415 - heavy; imported after the profile

    hidden = set()
    for ax in fig.axes:
        if not ax.axison:
            for axis in (ax.xaxis, ax.yaxis):
                labels = (*axis.get_ticklabels(), *axis.get_ticklabels(minor=True))
                hidden |= {id(t) for t in labels}
            continue
        for axis in (ax.xaxis, ax.yaxis):
            low, high = sorted(axis.get_view_interval())
            span = high - low
            ticks = [
                *zip(axis.get_major_ticks(), axis.get_majorticklocs(), strict=False),
                *zip(axis.get_minor_ticks(), axis.get_minorticklocs(), strict=False),
            ]
            for tick, loc in ticks:
                if not low - 1e-9 * span <= loc <= high + 1e-9 * span:
                    hidden |= {id(tick.label1), id(tick.label2)}
    return [
        t
        for t in fig.findobj(Text)
        if t.get_visible() and t.get_text().strip() and id(t) not in hidden
    ]


def thin_tick_labels(fig: Any) -> None:
    """Show every k-th tick label of an axis whose labels collide, with the smallest k that fits.

    A tick locator knows an axis's length, not its labels' size: 990,000 to 1,000,000 reads in
    steps of 2,000 on a short value axis stacked the labels on each other (2026-09-30).
    """
    from matplotlib.backends.backend_agg import FigureCanvasAgg  # noqa: PLC0415 - heavy

    renderer = FigureCanvasAgg(fig).get_renderer()  # type: ignore[no-untyped-call]
    fig.draw(renderer)
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            labels = [t for t in axis.get_ticklabels() if t.get_visible() and t.get_text()]
            boxes = [t.get_window_extent(renderer) for t in labels]
            step = 1
            while step < len(labels) and any(
                a.overlaps(b) for a, b in zip(boxes[::step], boxes[step::step], strict=False)
            ):
                step += 1
            for i, label in enumerate(labels):
                label.set_visible(i % step == 0)


def _refuse_marks(where: str, texts: Iterable[str]) -> None:
    """Refuse a glossary mark anywhere in a figure but its caption.

    Only the caption under the plot draws a mark as a term with its tooltip, and
    ``place_caption`` prints it as its display text in the PDF and SVG. A title, an axis
    label or a hover text would print the brackets in the plotly figure (D-048).

    Raises:
        ValueError: Naming the first marked text.
    """
    if found := next((t for t in texts if MARK.search(t)), None):
        raise ValueError(f"{where}: a glossary mark outside the caption: {found!r}")


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


# What makes two renders of one figure differ though nothing drawn did (656, Phase 2D): the
# PDF's CreationDate and the SVG's dc:date, set from the clock, and the SVG's clip-path and
# glyph ids, hashed with a random salt unless svg.hashsalt is set.
NO_DATE = {"pdf": {"CreationDate": None}, "svg": {"Date": None}}
SVG_SALT = "ctx"


def save_reproducibly(fig: Any, path: Path) -> None:
    """Save a PDF or SVG with no creation time and fixed ids, so a rerun is byte-equal."""
    import matplotlib as mpl  # noqa: PLC0415 - heavy; imported after the profile

    with mpl.rc_context({"svg.hashsalt": SVG_SALT}):
        fig.savefig(path, metadata=NO_DATE[path.suffix.lstrip(".")])


def figure_block(data: Figure, figure_stem: str, text: FigureText = NO_TEXT) -> FigureBlock:
    """The interactive figure, caption, data rows and words, for a figure or a module page.

    A figure with page controls carries them (a scatter's live calls, a linked view's brush,
    a bar's filters); a scatter with live calling also takes, unless ``text`` brings its own,
    the "How to read" lines at the default setting. A lead that follows its rows
    (``text.follows``) is filled here over the rows the page embeds, the rows its filter
    refills it from, and a figure with filters carries it for the page script.

    Raises:
        ValueError: If the figure's title or plotly text holds a glossary mark, or a scatter
            with live calls brings its own reading.
    """
    controls = data.controls() if isinstance(data, CONTROLLED) else None
    if text.follows is not None:
        rows = [finite(r) for r in data.page_table().to_dicts()]
        if controls is not None and controls["kind"] == "filter":
            controls = {**controls, "lead": text.follows.controls(rows)}
        text = replace(text, lead=text.follows.text(rows), follows=None)
    if isinstance(data, ScatterData) and data.calls:
        # The page script refills these lines by position from the figure's own templates;
        # lines from elsewhere would be overwritten, or leave a template with no line.
        if text.reading:
            raise ValueError("a figure with live calls takes its reading from calls.reading")
        text = replace(text, reading=data.reading())
    figure = data.interactive()
    below = data.below() if isinstance(data, LinkedData) else None
    drawn = [figure, *([below] if below is not None else [])]
    texts = (t for f in drawn for t in export.strings(f.to_dict()))
    _refuse_marks(figure_stem, [text.title, *texts])
    if controls is not None and isinstance(data, CONTROLLED):
        # A moved control makes the PDF and SVG buttons draw the view as set (D-033), with
        # the figure files' caption and what they draw that the page's plot lacks.
        layers = data.export_layers()
        controls = {**controls, "export": export.settings(data.inches(), data.caption,
                                                          figure_stem, layers)}  # fmt: skip
    return FigureBlock(
        figure, data.caption_text(), data.page_table().to_dicts(), f"{figure_stem}.csv",
        text=text, stem=figure_stem, controls=controls, below=below,
        reproduce=_panels().get(figure_stem),
    )  # fmt: skip


def write_figure(
    data: Figure,
    directory: Path,
    figure_stem: str,
    *,
    title: str,
    data_releases: Mapping[str, str],
    text: FigureText = NO_TEXT,
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
        text: Title, lead and reading guide; ``back`` links the page to its module page.

    Returns:
        The written files, in ``ctx.io.naming.figure_files()`` order.

    Raises:
        ValueError: Before any file is written: on an illegal stem, an undescribed column, a
            glossary mark outside the caption, or a text that leaves the page or overlaps
            another.
        GlossaryError: Before any file is written, if the caption marks a term the glossary
            does not hold.
    """
    from matplotlib.text import Text  # noqa: PLC0415 - heavy; imported after the profile

    if not check(f"{figure_stem}.pdf"):
        raise ValueError(f"not a legal output stem: {figure_stem!r}")
    rows = data.table()
    columns = dict(data.table_columns())
    validate(rows, figure_stem, columns, ("csv",))
    block = figure_block(data, figure_stem, text)
    _refuse_marks(figure_stem, [title])
    glossary.used([block.caption])
    figures = {theme.name: data.static(theme) for theme in THEMES.values()}
    for name, fig in figures.items():
        _refuse_marks(f"{figure_stem} ({name})", (t.get_text() for t in fig.findobj(Text)))
        thin_tick_labels(fig)
        if problems := layout_problems(fig):
            raise ValueError(f"{figure_stem} ({name}): " + "; ".join(problems[:5]))
    directory.mkdir(parents=True, exist_ok=True)
    rendered = {}
    for name, fig in figures.items():
        for ext in ("pdf", "svg"):
            rendered[(name, ext)] = directory / f"{figure_stem}_{name}.{ext}"
            save_reproducibly(fig, rendered[(name, ext)])
    table = write_table(
        rows, directory, figure_stem, columns=columns, data_releases=data_releases, formats=("csv",)
    )
    # The spec hashes the data table, so it is found once the table is written (Phase 2D).
    if (spec := record.describe(figure_stem, directory, table[0])) is not None:
        with table[-1].open("a") as readme:
            readme.write("\n" + record.readme(spec))
        # The module builds its page's block in a second call, figure_block, by this stem.
        _panels()[figure_stem] = panel.block(spec)
        block = replace(block, reproduce=_panels()[figure_stem])
    page = write(FigurePage(title, block), directory / f"{figure_stem}.html")
    return [*(rendered[(t, ext)] for ext in ("pdf", "svg") for t in THEMES), page, *table]
