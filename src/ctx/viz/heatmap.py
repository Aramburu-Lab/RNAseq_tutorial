"""Heatmap: every cell holds a level, every level has a swatch; vector cells, a swatch legend.

Two uses share the contract. A QC band map colours samples × metrics good, usable, bad or no
band. An expression map bins z-scores into a stepped diverging scale, two data colours at
three opacities (``diverging_levels``), orders rows and columns by average-linkage clustering
on 1 − Pearson (``cluster_order``), draws category tracks above the columns and a statistic
beside each row. ``static()`` and ``interactive()`` read the same ``HeatmapData``.
"""

import itertools
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import Any

import numpy as np
import polars as pl

from ctx.viz import metrics
from ctx.viz.annotate import place_caption
from ctx.viz.axes import axis_title
from ctx.viz.metrics import text_width, wrap_to
from ctx.viz.style import (
    BAR_LABEL_ROTATION,
    HEAT,
    HEAT_ALPHAS,
    LINE_SPACING,
    LINE_WIDTHS,
    MM_PER_INCH,
    PAGE_PLOT_PX,
    PHONE_PLOT_PX,
    PHONE_PX,
    PLOTLY_CAP_SHIFT,
    PLOTLY_HEAT_ROW_PX,
    PLOTLY_HOVER_PX,
    PLOTLY_LABEL_BOX_EM,
    PLOTLY_LEGEND_GAP_PX,
    PLOTLY_LEGEND_LINE_PX,
    PLOTLY_LEGEND_PAD_PX,
    PLOTLY_LEGEND_SWATCH_PX,
    PLOTLY_LEGEND_SYMBOL_PX,
    PLOTLY_MARGIN_PX,
    PLOTLY_NO_LEVEL,
    PLOTLY_PUSH_PAD_PX,
    PLOTLY_TICK_PX,
    PLOTLY_TITLE_STANDOFF_PX,
    ROLES,
    SIZES,
    THEMES,
    WIDTHS_MM,
    Theme,
    apply_matplotlib,
    on_page,
    plotly_steps,
    plotly_templates,
    tint,
)

# A blank row between the tracks and the cells, so a track never reads as a data row.
GAP = 1


@dataclass(frozen=True)
class Swatch:
    """A level's colour and opacity; None is the theme's grid grey, for "nothing to show"."""

    level: str
    colour: str | None
    alpha: float = 1.0


@dataclass(frozen=True)
class Track:
    """A category drawn as one strip of cells above the heatmap, one cell per column."""

    name: str
    values: Mapping[str, str]
    colours: Mapping[str, str]


def _suffix(n: int) -> str:
    """The ordinal suffix of n: 1st, 2nd, 3rd, 4th, 11th."""
    return "th" if n % 100 in (11, 12, 13) else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")


def cluster_order(matrix: np.ndarray, labels: Sequence[str]) -> tuple[str, ...]:
    """Row labels in the leaf order of average-linkage clustering on 1 − Pearson.

    A constant row has no correlation with anything; such rows follow the clustered ones, in
    their given order.
    """
    from scipy.cluster.hierarchy import leaves_list, linkage  # noqa: PLC0415 - heavy

    varying = [i for i in range(len(labels)) if np.std(matrix[i]) > 0]
    flat = [labels[i] for i in range(len(labels)) if i not in set(varying)]
    if len(varying) < 2:  # noqa: PLR2004 - linkage needs two rows
        return tuple(labels[i] for i in varying) + tuple(flat)
    tree = linkage(matrix[varying], method="average", metric="correlation")
    return tuple(labels[varying[i]] for i in leaves_list(tree)) + tuple(flat)


def diverging_levels(
    values: Sequence[float | None], edges: tuple[float, float, float] = (0.5, 1.0, 2.0)
) -> tuple[list[str | None], tuple[Swatch, ...]]:
    """Bin values into a stepped diverging scale: three steps each side of a blank middle.

    Below zero is the reference colour, above it the comparison colour; each step further
    from zero is more opaque (``HEAT_ALPHAS``). A missing value has no level.

    Returns:
        Each value's level, and the swatches from most negative to most positive.
    """
    e1, e2, e3 = edges
    names = [f"< -{e3:g}", f"-{e3:g} to -{e2:g}", f"-{e2:g} to -{e1:g}", f"-{e1:g} to {e1:g}",
             f"{e1:g} to {e2:g}", f"{e2:g} to {e3:g}", f"≥ {e3:g}"]  # fmt: skip
    swatches = (
        *(Swatch(names[i], ROLES["reference"], HEAT_ALPHAS[2 - i]) for i in range(3)),
        Swatch(names[3], None),
        *(Swatch(names[4 + i], ROLES["comparison"], HEAT_ALPHAS[i]) for i in range(3)),
    )
    bounds = (-e3, -e2, -e1, e1, e2, e3)
    levels = [None if v is None or np.isnan(v) else names[int(np.searchsorted(bounds, v, "right"))]
              for v in values]  # fmt: skip
    return levels, swatches


@dataclass(frozen=True)
class HeatmapData:
    """The contract both renderers read.

    Attributes:
        frame: Long table, one row per cell: the row, the column, the value and its level.
        columns: (unit, meaning) for every column of ``frame``.
        row: Column naming a cell's row, e.g. the gene or the sample.
        column: Column naming a cell's column.
        value: Column with the number a cell stands for (may be null).
        level: Column with the cell's level, a key of ``swatches`` (null draws nothing).
        row_label: What a row is.
        column_label: What a column is, as hover names it.
        swatches: Every level's colour, in legend order.
        caption: The caption, unwrapped.
        row_order: Rows top to bottom; empty keeps first appearance.
        column_order: Columns left to right; empty keeps first appearance.
        tracks: Categories drawn above the columns.
        row_stats: Row to the number printed at its right, e.g. an adjusted p.
        row_stats_label: What that number is.
        row_stats_format: Its format spec on the figure; the CSV keeps the number.
        width: The column; every figure takes ``FIGURE_WIDTH``, the double column.
        column_title: The column axis's title, where it says more than ``column_label``:
            what a column is, how many there are and in what order.
        row_display: Row to the label its axis prints, with a line break where it wraps, as
            a band map's metric over its unit and cut-offs; the row itself where absent.
        thin: Rows and columns name the same units (a sample against a sample): above
            ``HEAT["max_page_cells"]`` cells, the page draws every k-th unit in drawn order,
            its row and its column whole, with the cells between them. The files draw every
            cell; the CSV holds every cell.
    """

    frame: pl.DataFrame
    columns: Mapping[str, tuple[str, str]]
    row: str
    column: str
    value: str
    level: str
    row_label: str
    column_label: str
    swatches: tuple[Swatch, ...]
    caption: str
    row_order: tuple[str, ...] = ()
    column_order: tuple[str, ...] = ()
    tracks: tuple[Track, ...] = ()
    row_stats: Mapping[str, float] | None = None
    row_stats_label: str = ""
    row_stats_format: str = ".2g"
    width: str = "double"
    column_title: str = ""
    row_display: Mapping[str, str] | None = None
    thin: bool = False

    def __post_init__(self) -> None:
        """Refuse a contract the renderers could only draw wrongly."""
        if missing := [c for c in (self.row, self.column, self.value, self.level)
                       if c not in self.frame.columns]:  # fmt: skip
            raise ValueError(f"missing columns {missing}")
        levels = set(self.frame[self.level].drop_nulls())
        if stray := levels - {s.level for s in self.swatches}:
            raise ValueError(f"levels with no swatch: {sorted(stray)}")
        if self.frame.select(self.row, self.column).is_duplicated().any():
            raise ValueError("a cell appears twice")
        for track in self.tracks:
            if gaps := set(self.columns_in_order()) - set(track.values):
                raise ValueError(f"track {track.name!r} has no value for {sorted(gaps)[:5]}")
            if stray := set(track.values.values()) - set(track.colours):
                raise ValueError(f"track {track.name!r}: levels {sorted(stray)} have no colour")

    def rows_in_order(self) -> list[str]:
        """Rows top to bottom."""
        seen = self.frame[self.row].unique(maintain_order=True)
        return list(self.row_order) or [str(v) for v in seen]

    def columns_in_order(self) -> list[str]:
        """Columns left to right."""
        seen = self.frame[self.column].unique(maintain_order=True)
        return list(self.column_order) or [str(v) for v in seen]

    def table(self) -> pl.DataFrame:
        """Every cell, with its row's statistic and its column's track levels."""
        frame = self.frame.with_columns(pl.col(self.row).cast(pl.String),
                                        pl.col(self.column).cast(pl.String))  # fmt: skip
        if self.row_stats is not None:
            frame = frame.with_columns(
                row_stat=pl.col(self.row).replace_strict(
                    dict(self.row_stats), default=None, return_dtype=pl.Float64
                )
            )
        return frame.with_columns([
            pl.col(self.column).replace_strict(dict(t.values)).alias(f"track_{t.name}")
            for t in self.tracks
        ])  # fmt: skip

    def table_columns(self) -> Mapping[str, tuple[str, str]]:
        """(unit, meaning) for every column of ``table()``."""
        extra = {f"track_{t.name}": ("none", f"{t.name} of the column.") for t in self.tracks}
        if self.row_stats is not None:
            meaning = f"{self.row_stats_label}, printed beside the row."
            extra["row_stat"] = ("dimensionless", meaning)
        return {**self.columns, **extra}

    def page_table(self) -> pl.DataFrame:
        """The cells the page draws: every cell, or the thinned units' cells."""
        stride = self.stride()
        return self.table() if stride == 1 else self._paged(stride).table()

    def caption_text(self) -> str:
        """The caption, then the page's thinning sentence when it draws fewer units.

        The figure files draw every cell, so their caption is ``caption`` alone.
        """
        stride = self.stride()
        if stride == 1:
            return self.caption
        from ctx.report import text  # noqa: PLC0415 - report imports viz

        units = self.rows_in_order()
        sentence = text.fill("ui", "thinned.heatmap", drawn=len(units[::stride]),
                             units=len(units), nth=f"{stride}{_suffix(stride)}",
                             total=self.frame.height)  # fmt: skip
        return f"{self.caption} {sentence}"

    def stride(self) -> int:
        """Every how many units the page draws: 1 draws all; only a ``thin`` heatmap thins.

        The stride counts cells, the frame's rows: the smallest k whose units' cells fit
        ``HEAT["max_page_cells"]``.
        """
        if not self.thin or self.frame.height <= HEAT["max_page_cells"]:
            return 1
        stride = 2
        while self._paged(stride).frame.height > HEAT["max_page_cells"]:
            stride += 1
        return stride

    def _paged(self, stride: int) -> "HeatmapData":
        """The contract of every ``stride``-th unit in drawn order, its rows and columns whole."""
        kept = self.rows_in_order()[::stride]
        frame = self.frame.filter(pl.col(self.row).cast(pl.String).is_in(kept)
                                  & pl.col(self.column).cast(pl.String).is_in(kept))  # fmt: skip
        order = tuple(kept)
        # Still ``thin``: its units are one kind on both axes, so its cells stay square; it
        # fits the budget, so it is never thinned again.
        return replace(self, frame=frame, row_order=order, column_order=order)

    def _grid(self) -> tuple[np.ndarray, list[str]]:
        """Swatch index per cell, tracks on top, the gap, then the rows; NaN draws nothing.

        Track levels get swatch indices after the data swatches, so one colour list serves.
        """
        rows, cols = self.rows_in_order(), self.columns_in_order()
        index = {s.level: i for i, s in enumerate(self.swatches)}
        offset, track_rows = len(self.swatches), []
        for track in self.tracks:
            codes = {level: offset + i for i, level in enumerate(track.colours)}
            track_rows.append([codes[track.values[c]] for c in cols])
            offset += len(track.colours)
        cells = np.full((len(rows), len(cols)), np.nan)
        r, c = {v: i for i, v in enumerate(rows)}, {v: i for i, v in enumerate(cols)}
        for row, col, level in self.frame.select(self.row, self.column, self.level).rows():
            if level is not None:
                cells[r[str(row)], c[str(col)]] = index[level]
        gap = np.full((GAP if self.tracks else 0, len(cols)), np.nan)
        labels = [t.name for t in self.tracks] + [""] * len(gap) + rows
        return np.vstack([np.array(track_rows).reshape(-1, len(cols)), gap, cells]), labels

    def _rgba(self, theme: Theme) -> list[tuple[float, float, float, float]]:
        """The colour of every swatch index: data swatches, then each track's levels."""
        from matplotlib.colors import to_rgba  # noqa: PLC0415 - heavy

        colours = [to_rgba(s.colour or theme.grid, s.alpha) for s in self.swatches]
        return colours + [to_rgba(c) for t in self.tracks for c in t.colours.values()]

    def _display(self, row: str) -> str:
        """The label a row's axis prints; a line break where it wraps."""
        return (self.row_display or {}).get(row, row)

    def _shown(self) -> list[str]:
        """The row axis's labels top down as the static prints them, a line break where one wraps.

        The tracks, the blank row, then every row while there are ``HEAT["label_rows"]`` or
        fewer.
        """
        labels, rows = self._grid()[1], len(self.rows_in_order())
        shown = labels if rows <= HEAT["label_rows"] else labels[: len(labels) - rows]
        return [self._display(label) for label in shown]

    def _row_lines(self) -> int:
        """Lines in the tallest row label: each row is drawn that many lines high."""
        return max((self._display(r).count("\n") + 1 for r in self.rows_in_order()), default=1)

    def _column_label_length(self, size: float) -> float:
        """The longest column label set at ``size``, or 0 where columns are left unlabelled.

        A label runs upward, so this is how far it reaches below the cells; estimated at the
        average glyph it fell 38 px short of a 19-character run name in capitals (2026-10-01).
        """
        columns = self.columns_in_order()
        if len(columns) > HEAT["label_columns"]:
            return 0.0
        return max((text_width(c, size) for c in columns), default=0.0)

    def _inches(self) -> tuple[float, float]:
        """Width by the column rule; height by the rows drawn and the longest column label.

        The rows drawn are capped; a column label runs upward and, unmeasured, squeezed the rows.
        A ``_square`` heatmap's rows are as tall as its columns are wide at the full width,
        every row drawn: capped at 60 rows of 4.5 mm, Step-1's 105 drawn samples stood 2.9 times
        as tall as wide on the page (2026-10-09). The static then holds the cells square.
        """
        if self._square():
            rows_mm = self._static_edges()[-1] * self._pitch_mm()
        else:
            drawn = min(len(self.rows_in_order()), HEAT["max_rows_drawn"]) + len(self.tracks) + GAP
            rows_mm = drawn * HEAT["row_mm"] * self._row_lines()
        label_mm = self._column_label_length(SIZES["tick"]) / 72 * MM_PER_INCH
        return (WIDTHS_MM[self.width] / MM_PER_INCH,
                (rows_mm + HEAT["margin_mm"] + label_mm) / MM_PER_INCH)  # fmt: skip

    def _square(self) -> bool:
        """Whether the cells are drawn square: one unit on both axes, every row holding its label.

        A row holds its label when it goes unlabelled or a column's width is a label's height.

        A labelled row a square cell cannot hold keeps the height its label needs: 40 rows at
        170 mm are 4.25 mm, under a 9 pt label's 4.5 mm. Step-1's 8-sample laptop heatmap,
        21.25 mm a row, once stopped its page loading held square, until the page script
        stopped on an alternation of heights as well as a repeat (2026-10-14).
        """
        if not self.thin:
            return False
        unlabelled = len(self.rows_in_order()) > HEAT["label_rows"]
        return unlabelled or self._pitch_mm() >= HEAT["row_mm"] * self._row_lines()

    def _pitch_mm(self) -> float:
        """A column's width in mm were the cells to take the figure's whole width."""
        return WIDTHS_MM[self.width] / len(self.columns_in_order())

    def _static_edges(self) -> np.ndarray:
        """The static's row edges top down, a data row one unit tall.

        A ``_square`` heatmap's cells can be shorter than a label: its track and gap rows stand
        ``HEAT["row_mm"]`` tall, so two track names never overlap; every other row is one.
        """
        top = len(self.tracks) + (GAP if self.tracks else 0)
        heights = np.ones(top + len(self.rows_in_order()))
        if self._square():
            heights[:top] = max(1.0, HEAT["row_mm"] / self._pitch_mm())
        return np.concatenate([[0.0], np.cumsum(heights)])

    def _column_title(self, width: float) -> str:
        """The column axis's title, wrapped to the width the cells leave beside the row labels.

        A long title centred under cells pushed right by long row labels ran off the page. The
        labels and the title are measured as the files set them: estimated at an average
        glyph, the Step-1 teaching band map's title broke a word early (2026-10-07).
        """
        lines = [line for label in self._shown() for line in label.split("\n")]
        labels = max((text_width(line, SIZES["tick"]) for line in lines), default=0.0)
        title = self.column_title or self.column_label
        return "\n".join(wrap_to(metrics.words(title), width * 72 - labels, SIZES["axis_label"]))

    def static(self, theme: Theme) -> Any:
        """Render the matplotlib figure for one theme; cells are vector quads.

        Args:
            theme: One of ``style.THEMES``.

        Returns:
            A matplotlib Figure, not attached to pyplot.
        """
        from matplotlib.colors import BoundaryNorm, ListedColormap  # noqa: PLC0415 - heavy
        from matplotlib.figure import Figure  # noqa: PLC0415 - heavy
        from matplotlib.patches import Patch  # noqa: PLC0415 - heavy

        apply_matplotlib(theme)
        width, height = self._inches()
        fig = Figure(figsize=(width, height))
        ax = fig.subplots()
        grid, labels = self._grid()
        colours = self._rgba(theme)
        ax.pcolormesh(np.arange(grid.shape[1] + 1), self._static_edges(),
                      np.ma.masked_invalid(grid), cmap=ListedColormap(colours),
                      norm=BoundaryNorm(np.arange(len(colours) + 1) - 0.5, len(colours)),
                      linewidth=LINE_WIDTHS["band_edge"])  # fmt: skip
        ax.invert_yaxis()
        if self._square():
            # One kind of unit on both axes: a cell as tall as it is wide, the cells at the top.
            ax.set_aspect("equal", anchor="N")
        self._label_axes(ax, labels)
        ax.set_xlabel(self._column_title(width))
        for spine in ax.spines.values():
            spine.set_visible(False)
        handles = [Patch(facecolor=c, label=s.level) for s, c in zip(self.swatches, colours,
                                                                      strict=False)]  # fmt: skip
        handles += [Patch(facecolor=c, label=f"{t.name}: {level}")
                    for t in self.tracks for level, c in t.colours.items()]  # fmt: skip
        place_caption(fig, self.caption, width, height, ax=ax, handles=handles)
        return fig

    def _label_axes(self, ax: Any, labels: list[str]) -> None:
        """Row and track names at left, the statistic at right, column names below.

        Past ``HEAT``'s label limits a side is left unlabelled; hover and the CSV name every
        cell.
        """
        rows, cols = self.rows_in_order(), self.columns_in_order()
        shown = self._shown()
        edges = self._static_edges()
        centres = (edges[:-1] + edges[1:]) / 2
        ax.set_yticks(centres[: len(shown)], shown)
        ax.tick_params(axis="both", length=0)
        if len(cols) <= HEAT["label_columns"]:
            ax.set_xticks(np.arange(len(cols)) + 0.5, cols, rotation=BAR_LABEL_ROTATION)
        else:
            ax.set_xticks([])
        if self.row_stats is not None and len(rows) <= HEAT["label_rows"]:
            right = ax.secondary_yaxis("right")
            start = len(labels) - len(rows)
            right.set_yticks(centres[start:], [self._stat(r) for r in rows])
            right.tick_params(axis="y", length=0)
            right.set_ylabel(self.row_stats_label)
            right.spines["right"].set_visible(False)

    def cell_aspect(self) -> float:
        """A cell's height over its width as the static draws it, read off its laid-out axes.

        The static's height adds a fixed margin to its rows and matplotlib gives the cells what
        the legend, labels and caption leave, so the aspect is the layout's outcome, not a
        constant; the page takes it from here (D13). At a fixed pitch, kit v6's M12 page drew
        its cells 0.34 high per unit of width at 1,440 px, where its files draw 0.44
        (2026-10-07).
        """
        fig = self.static(THEMES["white"])
        box = fig.axes[0].get_position()
        width, height = fig.get_size_inches()
        rows, columns = self._static_edges()[-1], len(self.columns_in_order())
        return float((box.height * height / rows) / (box.width * width / columns))

    def interactive(self) -> Any:
        """Render the plotly figure in the static's layout: tracks, the gap row, then the rows.

        Every row of ``_grid`` is a band of a numeric axis with one tick, so a track stands one
        row above the cells as in the static and no label is dropped; every column the static
        labels is labelled. The bands take the static's cell aspect at a module page's width
        (``_wide``), so the cells read as in the files. The row statistic is a second axis at
        right, which automargin keeps inside the figure. Narrower than its labels need, the
        plot scrolls sideways under the page's hint (``meta.least_width``); a phone draws the
        labels wrapped to its plot, each row as tall as its label (``_narrow``).

        Returns:
            A plotly Figure with the light template; the page switches themes.
        """
        import plotly.graph_objects as go  # noqa: PLC0415 - heavy; imported after the profile

        if (stride := self.stride()) > 1:
            return self._paged(stride).interactive()
        rows, cols = self.rows_in_order(), self.columns_in_order()
        words = self._words(phone=False)
        wide = self._wide(words)
        words = replace(words, title=self._title(words, wide.below if self.column_title
                                                 else wide.width))  # fmt: skip
        edges = wide.edges
        first = len(edges) - 1 - len(rows)
        fig = go.Figure(self._cell_trace(go, rows, cols, edges[first:]))
        for at, track in enumerate(self.tracks):
            fig.add_trace(self._track_trace(go, track, cols, edges[at : at + 2]))
        # One legend group per entry: plotly lays ungrouped entries out at the widest one's
        # width, one a line on an overview's plot (2026-10-01); groups flow by their own widths.
        for i, (name, colour) in enumerate(zip(words.legend, self._colours(), strict=True)):
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode="markers", name=name, legendgroup=str(i),
                marker={"symbol": "square", "color": colour, "size": PLOTLY_LEGEND_SWATCH_PX},
            ))  # fmt: skip
        fig.update_layout(
            template=plotly_templates()["white"], height=self._wide_height(wide, words),
            showlegend=True,
            # At the div's top left, not the plot's, which a phone's labels push right; never
            # scrolled, as plotly scrolls a legend past half the figure's height.
            legend={"xref": "container", "yref": "container", "x": 0, "y": 1,
                    "yanchor": "top", "maxheight": 1, "traceorder": "grouped",
                    "tracegroupgap": 0},
            **self._axes(edges, words, wide),
        )  # fmt: skip
        if self._shows_stats():
            # plotly draws no axis that no trace is on, so an empty one holds it.
            fig.add_trace(go.Scatter(x=[None], y=[None], yaxis="y2", showlegend=False,
                                     hoverinfo="skip"))  # fmt: skip
            fig.update_layout(yaxis2={
                **_BARE, "matches": "y", "overlaying": "y", "side": "right", "tickmode": "array",
                "tickvals": _centres(edges)[first:], "ticktext": [self._stat(r) for r in rows],
                "title": {"text": self.row_stats_label},
            })  # fmt: skip
        if self.column_title:
            fig.add_annotation(self._title_note(words.title))
            fig.update_layout(margin_b=self._bottom(words.title))
        fig.update_layout(meta=self._meta(wide))
        return fig

    def _axes(self, edges: list[float], words: "_Words", wide: "_Wide") -> dict[str, Any]:
        """The column and row axes of the static's layout.

        Columns sit at their positions on a numeric axis, each named by its tick, as the rows
        do: plotly thins a categorical axis's labels to the width it draws at, and M7's page
        labelled every other column of 21; and plotly holds an aspect only between two axes of
        one type. Where ``wide.held``, the row axis holds the bands at ``wide.ratio`` of a
        column's width, its length shrinking toward the legend where the plot is narrower
        than ``wide.width``, and the column axis's toward the row labels where it is wider.
        """
        cols = self.columns_in_order()
        labelled = self._labelled_columns() > 0
        held = {"scaleanchor": "x", "scaleratio": wide.ratio, "constrain": "domain",
                "constraintoward": "top"} if wide.held else {}  # fmt: skip
        return {
            "xaxis": {**_BARE, "range": [-0.5, len(cols) - 0.5], "tickmode": "array",
                      "tickvals": list(range(len(cols))) if labelled else [],
                      "ticktext": cols if labelled else [],
                      "title": {"text": "" if self.column_title else words.title},
                      "tickangle": -BAR_LABEL_ROTATION, "showticklabels": labelled,
                      "constrain": "domain", "constraintoward": "left"},
            "yaxis": {**_BARE, "tickmode": "array", **self._yaxis(edges, words), **held},
        }  # fmt: skip

    def _meta(self, wide: "_Wide") -> dict[str, Any]:
        """The least width the page draws the plot at, and a phone's layout where it can show.

        The least is what a phone's labels need, the narrowest layout's; plotly's own padding
        doubles to 12 px, not 16, opened alone, so twice it holds the labels there too.
        """
        phone = self._words(phone=True)
        least = self._needs(phone)
        meta: dict[str, Any] = {"least_width": [[least, 0]]}
        if wide.below > least:
            meta["narrow"] = self._narrow(phone, wide, least)
        return meta

    def _words(self, phone: bool) -> "_Words":
        """The static's row labels and legend entries, or a phone's; the title is ``_title``'s.

        On a phone each wraps to the room it has there, a row label to half the plot, an entry
        to the plot beside its swatch.
        """
        ticks = {i: "<br>".join(_fit(line, PHONE_PLOT_PX / 2) if phone else line
                                for line in label.split("\n"))
                 for i, label in enumerate(self._shown()) if label}  # fmt: skip
        legend = [s.level for s in self.swatches]
        legend += [f"{t.name}: {level}" for t in self.tracks for level in t.colours]
        if phone:
            room = PHONE_PLOT_PX - PLOTLY_LEGEND_SYMBOL_PX
            legend = [_fit(name, room, SIZES["legend"]) for name in legend]
        return _Words(ticks, legend, "")

    def _title(self, words: "_Words", width: float) -> str:
        """The column title as a layout prints it at its narrowest, ``width`` px.

        ``column_title`` is right-aligned under the column labels, from the figure's left
        margin to the plot's right edge. ``column_label`` is the axis title, centred under the
        plot, which the page script wraps to twice the room from the plot's centre to the
        nearer edge of the figure: a phone's layout was sized for it wrapped to half the plot,
        while the script wrapped it wider (2026-10-07).
        """
        right = math.ceil(self._right_px())
        if self.column_title:
            room = width - PLOTLY_MARGIN_PX["l"] - right
            return _fit(metrics.words(self.column_title), room, SIZES["axis_label"])
        left = math.ceil(self._left_px(words.ticks))
        return _fit(metrics.words(self.column_label), width - abs(left - right),
                    SIZES["axis_label"])  # fmt: skip

    def _wide(self, words: "_Words") -> "_Wide":
        """The static's layout on the page, laid out at the width a module page draws a plot at.

        A row's band is the static's cell height for the page's cell width there, every band
        that tall as in the static, and no band less than keeps two labels apart: the page's
        labels are larger beside its plot than the files' beside theirs, and the files' aspect,
        0.945, would put M7's 22.2 px labels in rows of 21.2 px at its 687 px. Past
        ``HEAT["label_rows"]`` the rows are unlabelled and keep the files' aspect; a track's
        band alone grows to its label. The page draws the layout from ``below`` up, where its
        labels fit and the plot is wider than a phone's screen. plotly holds the aspect at
        every such width where each band still keeps its label at the narrowest (``held``).
        Elsewhere, as on the band maps whose two-line labels would have pushed every desktop
        narrower than 686 px onto a phone's layout, the bands keep their heights; there the
        legend and title can take more lines, and a band that holds a label grows by what that
        takes from it, so it keeps the label at the narrowest width too.
        """
        n = len(self.columns_in_order())
        sides = self._sides(words.ticks)
        needs = self._needs(words)
        width = max(PAGE_PLOT_PX["beside"], needs)
        column = (width - sides) / n
        floors = self._floors(words.ticks)
        rows = len(self.rows_in_order())
        band = max(self.cell_aspect() * column, *floors[len(floors) - rows :])
        below = max(needs, PHONE_PX)
        scale = (below - sides) / n / column
        # A square heatmap's track names keep their floor down to the narrowest width, so plotly
        # holds its square cells at every width: unheld, Step-1's 105 samples kept the bands of
        # 686 px and drew cells wider than tall at 974 px (2026-10-09).
        top = len(floors) - rows if self._square() else 0
        heights = [max(band, floor / scale if i < top else floor) for i, floor in enumerate(floors)]
        held = all(h * scale >= floor for h, floor in zip(heights, floors, strict=True))
        if not held:
            squeeze = self._chrome(words, below, below) - self._chrome(words, width, below)
            grow = max([1.0, *(f / h + squeeze / sum(heights)
                               for h, f in zip(heights, floors, strict=True) if f)])  # fmt: skip
            band, heights = band * grow, [h * grow for h in heights]
        edges = [0.0, *itertools.accumulate(h / band for h in heights)]
        return _Wide(width, band, band / column, below, held, edges)

    def _chrome(self, words: "_Words", width: float, below: float) -> float:
        """The figure's height outside its bands in px, at a plot ``width`` px wide.

        The legend as plotly flows it over the cells there, then the column labels and the
        title: an axis title as the page script wraps it there, a ``column_title`` as wrapped
        for the narrowest plot, ``below``.
        """
        title = self._title(words, below if self.column_title else width)
        legend = self._legend_px(words.legend, width - self._sides(words.ticks))
        return self._top(legend) + math.ceil(self._bottom(title))

    def _wide_height(self, wide: "_Wide", words: "_Words") -> float:
        """The figure's height in px at ``wide.width``: its legend, bands, labels and title."""
        return self._chrome(words, wide.width, wide.below) + wide.edges[-1] * wide.band

    def _edges(self, words: "_Words") -> list[float]:
        """A phone's band edges top down, in row pitches: each as many as its label's lines.

        A ``_square`` heatmap's rows are as tall as its columns are wide on a phone's plot: at
        a label's pitch, Step-1's 105 drawn samples stood 8 times as tall as wide at 375 px
        (2026-10-09).
        """
        heights = [float(_lines(words.ticks.get(i, ""))) for i in range(len(self._grid()[1]))]
        if self._square():
            rows, columns = len(self.rows_in_order()), len(self.columns_in_order())
            plot = PHONE_PLOT_PX - self._left_px(words.ticks) - self._right_px()
            heights[len(heights) - rows :] = [plot / columns / PLOTLY_HEAT_ROW_PX] * rows
        return [0.0, *itertools.accumulate(heights)]

    def _yaxis(self, edges: list[float], words: "_Words") -> dict[str, Any]:
        """The row axis's range and ticks: one tick at the middle of each labelled band."""
        centres = _centres(edges)
        return {"range": [edges[-1], edges[0]], "tickvals": [centres[i] for i in words.ticks],
                "ticktext": list(words.ticks.values())}  # fmt: skip

    def _narrow(self, words: "_Words", wide: "_Wide", least: float) -> dict[str, Any]:
        """A phone's layout, which the page script swaps in narrower than ``wide.below`` px.

        At 375 px the band map's two-line labels started up to 132 px left of the plot and lost
        their unit (2026-10-01); wrapped for a phone on every screen, a desktop's rows were as
        tall as a phone's labels (work order item 2, 2026-10-01). Each label wraps to a phone's
        plot and each row grows to its label's lines, free of the static's aspect: new band
        edges for the heatmap traces, names for the legend's (``restyle``), ticks, title, the
        aspect's anchor and the height (``relayout``).
        """
        words = replace(words, title=self._title(words, max(PHONE_PLOT_PX, least)))
        edges = self._edges(words)
        first = len(edges) - 1 - len(self.rows_in_order())
        relayout = {f"yaxis.{key}": value for key, value in self._yaxis(edges, words).items()}
        if wide.held:
            relayout["yaxis.scaleanchor"] = None
        if self._shows_stats():
            relayout["yaxis2.tickvals"] = _centres(edges)[first:]
        if self.column_title:
            relayout |= {"annotations[0].text": words.title, "margin.b": self._bottom(words.title)}
        else:
            relayout["xaxis.title.text"] = words.title
        relayout["height"] = self._page_height(edges, words)
        heat = list(range(len(self.tracks) + 1))
        y = [edges[first:], *(edges[at : at + 2] for at in range(len(self.tracks)))]
        keys = list(range(len(heat), len(heat) + len(words.legend)))
        return {"below": wide.below, "relayout": relayout,
                "restyle": [[{"y": y}, heat], [{"name": words.legend}, keys]]}  # fmt: skip

    def _labelled_columns(self) -> int:
        """How many columns the static labels: all of them, up to ``HEAT["label_columns"]``."""
        n = len(self.columns_in_order())
        return n if n <= HEAT["label_columns"] else 0

    def _needs(self, words: "_Words") -> int:
        """The narrowest plot that prints ``words`` whole, in px.

        The row labels' margin, a label's box across for each labelled column, upright, and the
        statistic's margin; or the widest legend entry beside its swatch. A column under its
        label's box overlapped its neighbour's.
        """
        columns = self._labelled_columns() * PLOTLY_LABEL_BOX_EM * on_page(SIZES["tick"])
        size = on_page(SIZES["legend"])
        entries = max((text_width(line, size) for name in words.legend
                       for line in name.split("<br>")), default=0.0)  # fmt: skip
        return math.ceil(max(self._sides(words.ticks) + columns, entries + PLOTLY_LEGEND_SYMBOL_PX))

    def _floors(self, ticks: Mapping[int, str]) -> list[float]:
        """Each band's least height in px: what keeps its label off a neighbour's.

        For a labelled band with a labelled neighbour, the box of the tallest label's lines, as
        the static gives every row the tallest label's; none for any other band.
        """
        lines = max((_lines(t) for t in ticks.values()), default=1)
        box = ((lines - 1) * LINE_SPACING + PLOTLY_LABEL_BOX_EM) * on_page(SIZES["tick"])
        return [box if i in ticks and (i - 1 in ticks or i + 1 in ticks) else 0.0
                for i in range(len(self._grid()[1]))]  # fmt: skip

    def _shows_stats(self) -> bool:
        """Whether the page prints the row statistic, as the static does, on a right axis."""
        return self.row_stats is not None and len(self.rows_in_order()) <= HEAT["label_rows"]

    def _left_px(self, ticks: Mapping[int, str]) -> float:
        """The left margin plotly draws, in px; the template's where no row is labelled.

        With row labels: its padding, half the axis line, the labels' standoff and the widest
        label line (plotly.js autoMargin).
        """
        tick = on_page(SIZES["tick"])
        widest = max((text_width(line, tick) for t in ticks.values()
                      for line in t.split("<br>")), default=0.0)  # fmt: skip
        if not widest:
            return PLOTLY_MARGIN_PX["l"]
        labels = LINE_WIDTHS["axes"] / 2 + PLOTLY_TICK_PX + widest
        return max(PLOTLY_MARGIN_PX["l"], PLOTLY_PUSH_PAD_PX + labels)

    def _right_px(self) -> float:
        """The right margin plotly draws, in px; the template's where no statistic is printed.

        With the statistic's axis: its padding, half the axis line, the labels' standoff, the
        widest statistic, then the title's cap height and standoff (plotly.js autoMargin and
        axes.drawTitle). Counted as a whole line of title on top of the template's margin, it
        came to 93.8 px on M7's page where plotly drew 82 (2026-10-07).
        """
        if not self._shows_stats():
            return PLOTLY_MARGIN_PX["r"]
        tick, title = on_page(SIZES["tick"]), on_page(SIZES["axis_label"])
        widest = max((text_width(self._stat(r), tick) for r in self.rows_in_order()), default=0.0)
        labels = LINE_WIDTHS["axes"] / 2 + PLOTLY_TICK_PX + widest
        return PLOTLY_PUSH_PAD_PX + labels + PLOTLY_CAP_SHIFT * title + PLOTLY_TITLE_STANDOFF_PX

    def _sides(self, ticks: Mapping[int, str]) -> int:
        """The left and right margins together, as plotly rounds each, in whole px."""
        return math.ceil(self._left_px(ticks)) + math.ceil(self._right_px())

    def _title_drop(self) -> float:
        """How far under the cells a right-aligned title starts, in px.

        The longest column label as set, plus the tick, plotly's own gap before a label and the
        title's standoff. Estimated at the average glyph it fell 38 px short of 19-character run
        names in capitals (2026-10-01).
        """
        labels = self._column_label_length(on_page(SIZES["tick"]))
        return labels + PLOTLY_TICK_PX + 2 * PLOTLY_TITLE_STANDOFF_PX

    def _title_note(self, text: str) -> dict[str, Any]:
        """``column_title`` as an annotation right-aligned under the column labels.

        Centred under cells that long row labels push to the right, it ran off the plot and was
        clipped (Rytis 2026-10-01). Placed on the axes' domains, it follows the cells where
        plotly shortens the row axis to hold their aspect.
        """
        return {"text": text, "xref": "x domain", "x": 1, "xanchor": "right",
                "yref": "y domain", "y": 0, "yanchor": "top", "yshift": -self._title_drop(),
                "align": "right", "showarrow": False,
                "font": {"size": on_page(SIZES["axis_label"])}}  # fmt: skip

    def _labels_depth(self) -> float:
        """How far the upright column labels reach under the cells, in px; 0 when unlabelled.

        plotly.js turns each about a point half its size above its anchor, which sits its cap
        shift under the axis line past the labels' standoff (axes.makeLabelFns).
        """
        longest = self._column_label_length(on_page(SIZES["tick"]))
        if not longest:
            return 0.0
        offset = (PLOTLY_CAP_SHIFT - 0.5) * on_page(SIZES["tick"])
        return LINE_WIDTHS["axes"] / 2 + PLOTLY_TICK_PX + offset + longest

    def _bottom(self, title: str) -> float:
        """The page plot's height under its cells, in px: column labels, title and margins.

        A right-aligned ``column_title`` sits ``_title_drop`` under the cells, in a bottom
        margin set to hold it. An axis title takes the margin plotly draws: its padding, the
        labels' depth, the title's cap height and its standoff (plotly.js autoMargin).
        """
        size = on_page(SIZES["axis_label"])
        if self.column_title:
            title_px = _lines(title) * size * LINE_SPACING
            return PLOTLY_MARGIN_PX["b"] + self._title_drop() + title_px
        caption = size * (PLOTLY_CAP_SHIFT + (_lines(title) - 1) * LINE_SPACING)
        pushed = PLOTLY_PUSH_PAD_PX + self._labels_depth() + caption + PLOTLY_TITLE_STANDOFF_PX
        return max(PLOTLY_MARGIN_PX["b"], pushed)

    def _top(self, legend: float) -> int:
        """The top margin plotly draws over a legend ``legend`` px tall, in whole px."""
        return max(PLOTLY_MARGIN_PX["t"], math.ceil(legend + PLOTLY_PUSH_PAD_PX))

    def _legend_px(self, names: list[str], room: float) -> float:
        """The legend's height as plotly lays it out in a plot ``room`` px wide.

        Entries flow left to right, each its swatch, its text and a gap, onto a new line where
        the next would cross ``room``; each line is as tall as its tallest entry
        (legend/draw.js, grouped and horizontal).
        """
        size = on_page(SIZES["legend"])
        heights, used, tallest = [], 0.0, 0.0
        for name in names:
            lines = name.split("<br>")
            entry = (PLOTLY_LEGEND_SYMBOL_PX + max(text_width(line, size) for line in lines)
                     + PLOTLY_LEGEND_GAP_PX)  # fmt: skip
            if used and used + entry > room:
                heights.append(tallest)
                used, tallest = 0.0, 0.0
            used += entry
            tallest = max(tallest, PLOTLY_LEGEND_LINE_PX + (len(lines) - 1) * size * LINE_SPACING)
        return sum([*heights, tallest]) + PLOTLY_LEGEND_PAD_PX

    def _page_height(self, edges: list[float], words: "_Words") -> float:
        """A phone's figure height in px: its rows at the label pitch, then the rest.

        As ``_inches`` does, the rows drawn are capped and the column labels measured. The
        legend has a line per entry, the most a phone stacks (M7's nine on nine lines at
        375 px, 2026-10-01); a wider plot gives the rows what it leaves.
        """
        rows = len(self.rows_in_order())
        above = edges[len(edges) - 1 - rows]  # the tracks and the blank row
        whole = rows <= HEAT["max_rows_drawn"] or self._square()
        cells = edges[-1] - above if whole else HEAT["max_rows_drawn"]
        lines = sum(_lines(name) for name in words.legend)
        top = self._top(lines * PLOTLY_LEGEND_LINE_PX + PLOTLY_LEGEND_PAD_PX)
        return top + (above + cells) * PLOTLY_HEAT_ROW_PX + math.ceil(self._bottom(words.title))

    def _cell_trace(self, go: Any, rows: list[str], cols: list[str], edges: list[float]) -> Any:
        """The cells as one heatmap trace between the row ``edges``; hover names each cell."""
        index = {s.level: i for i, s in enumerate(self.swatches)}
        drawn = self.frame.select(self.row, self.column, self.level, self.value).rows()
        cells = {(str(r), str(c)): (level, value) for r, c, level, value in drawn}
        z, text = [], []
        for r in rows:
            found = [cells.get((r, c), (None, None)) for c in cols]
            z.append([None if level is None else index[level] for level, _ in found])
            text.append([self._hover(r, c, *cell) for c, cell in zip(cols, found, strict=True)])
        return go.Heatmap(
            z=z, x=list(range(len(cols))), y=edges, text=text, showscale=False, xgap=0, ygap=0,
            hoverongaps=False, zmin=-0.5, zmax=len(self.swatches) - 0.5,
            colorscale=plotly_steps([_page_colour(s) for s in self.swatches]),
            hovertemplate="%{text}<extra></extra>",
        )  # fmt: skip

    def _track_trace(self, go: Any, track: Track, cols: list[str], edges: list[float]) -> Any:
        """One track as a heatmap strip between its two ``edges``; hover names its level."""
        codes = {level: i for i, level in enumerate(track.colours)}
        hover = [_hover_lines([_named(track.name, track.values[c]), _named(self.column_label, c)])
                 for c in cols]  # fmt: skip
        return go.Heatmap(
            z=[[codes[track.values[c]] for c in cols]], x=list(range(len(cols))), y=edges,
            text=[hover],
            showscale=False, colorscale=plotly_steps(list(track.colours.values())),
            zmin=-0.5, zmax=len(codes) - 0.5, hovertemplate="%{text}<extra></extra>",
        )  # fmt: skip

    def _hover(self, row: str, column: str, level: str | None, value: float | None) -> str:
        """A cell's hover: its row and column, its number named with its unit, then its level.

        A row is named on one line. A band map's cells carry a band and no value; their hover
        read "None" (2026-10-01).
        """
        lines = [_named(self.row_label, row), _named(self.column_label, column)]
        if value is not None and not np.isnan(value):
            unit = self.columns[self.value][0]
            lines.append(_named(axis_title(self.value, unit), f"{value:.3g}"))
        return _hover_lines([*lines, *([level] if level is not None else [])])

    def _colours(self) -> list[str]:
        """The static legend's colours on the page, in its order."""
        return [*(_page_colour(s) for s in self.swatches),
                *(c for t in self.tracks for c in t.colours.values())]  # fmt: skip

    def _stat(self, row: str) -> str:
        """A row's statistic as the figure prints it; empty where it has none."""
        value = (self.row_stats or {}).get(row)
        return "" if value is None else f"{value:{self.row_stats_format}}"


# An axis drawn bare: no line, ticks or grid; labels clear of the cells; margins grown to hold
# them.
_BARE = {"showline": False, "ticks": "", "ticklabelstandoff": PLOTLY_TICK_PX, "zeroline": False,
         "showgrid": False, "automargin": True}  # fmt: skip


@dataclass(frozen=True)
class _Words:
    """What one page layout prints, each text with ``<br>`` where it breaks.

    Row labels by grid row, legend entries and the column title.
    """

    ticks: Mapping[int, str]
    legend: list[str]
    title: str


@dataclass(frozen=True)
class _Wide:
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


def _page_colour(swatch: Swatch) -> str:
    """A swatch's colour on the page; "nothing to show" is a grey that reads in both themes."""
    colour, alpha = (swatch.colour, 1.0) if swatch.colour else PLOTLY_NO_LEVEL
    return tint(colour, alpha * swatch.alpha)


def _centres(edges: list[float]) -> list[float]:
    """The middle of each band between consecutive ``edges``."""
    return [(top + bottom) / 2 for top, bottom in itertools.pairwise(edges)]


def _named(label: str, value: str) -> str:
    """``label: value``; a label written as a sentence loses its full stop first."""
    return f"{label.removesuffix('.')}: {value}"


def _hover_lines(lines: list[str]) -> str:
    """Hover lines, each wrapped to a phone's plot, which a wider hover box ran off."""
    return "<br>".join(_fit(line, PLOTLY_HOVER_PX) for line in lines)


def _fit(text: str | list[str], width_px: float, size: float = SIZES["tick"]) -> str:
    """Text wrapped at spaces to lines of at most ``width_px`` on the page, in plotly's ``<br>``.

    Measured as ``ctx.viz.metrics`` measures every wrapped text. A word longer than a line
    stays whole: split, it would read as two words. A list is the text's words already, as
    ``metrics.words`` holds a parenthesised unit whole.
    """
    lines = wrap_to(text, width_px, on_page(size))
    return "<br>".join(lines) if lines else (text if isinstance(text, str) else "")


def _lines(text: str) -> int:
    """How many lines a plotly text runs to."""
    return text.count("<br>") + 1
