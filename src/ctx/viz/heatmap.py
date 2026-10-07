"""Heatmap: every cell holds a level, every level has a swatch; vector cells, a swatch legend.

Two uses share the contract. A QC band map colours samples × metrics good, usable, bad or no
band. An expression map bins z-scores into a stepped diverging scale, two data colours at
three opacities (``diverging_levels``), orders rows and columns by average-linkage clustering
on 1 − Pearson (``cluster_order``), draws category tracks above the columns and a statistic
beside each row. ``static()`` and ``interactive()`` read the same ``HeatmapData``.
"""

import itertools
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import polars as pl

from ctx.viz.annotate import place_caption, wrap
from ctx.viz.axes import axis_title
from ctx.viz.metrics import text_width, wrap_to
from ctx.viz.style import (
    AVERAGE_CHAR_EM,
    BAR_LABEL_ROTATION,
    HEAT,
    HEAT_ALPHAS,
    LINE_SPACING,
    LINE_WIDTHS,
    MM_PER_INCH,
    PHONE_PLOT_PX,
    PLOTLY_HEAT_BASE_PX,
    PLOTLY_HEAT_ROW_PX,
    PLOTLY_HOVER_PX,
    PLOTLY_LEGEND_LINE_PX,
    PLOTLY_LEGEND_PAD_PX,
    PLOTLY_LEGEND_SWATCH_PX,
    PLOTLY_LEGEND_SYMBOL_PX,
    PLOTLY_MARGIN_PX,
    PLOTLY_NO_LEVEL,
    PLOTLY_TICK_PX,
    PLOTLY_TITLE_STANDOFF_PX,
    ROLES,
    SIZES,
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
        """Every cell: the page draws them all."""
        return self.table()

    def caption_text(self) -> str:
        """The caption, unwrapped."""
        return self.caption

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
        """
        drawn = min(len(self.rows_in_order()), HEAT["max_rows_drawn"]) + len(self.tracks) + GAP
        label_mm = self._column_label_length(SIZES["tick"]) / 72 * MM_PER_INCH
        rows_mm = drawn * HEAT["row_mm"] * self._row_lines()
        return (WIDTHS_MM[self.width] / MM_PER_INCH,
                (rows_mm + HEAT["margin_mm"] + label_mm) / MM_PER_INCH)  # fmt: skip

    def _column_title(self, width: float) -> str:
        """The column axis's title, wrapped to the width the cells leave beside the row labels.

        A long title centred under cells pushed right by long row labels ran off the page.
        """
        lines = [line for r in self.rows_in_order() for line in self._display(r).split("\n")]
        labels_pt = max(map(len, lines), default=0) * SIZES["tick"] * AVERAGE_CHAR_EM
        room = width * 72 - labels_pt
        return wrap(self.column_title or self.column_label,
                    max(int(room / (SIZES["axis_label"] * AVERAGE_CHAR_EM)), 12))  # fmt: skip

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
        ax.pcolormesh(np.ma.masked_invalid(grid), cmap=ListedColormap(colours),
                      norm=BoundaryNorm(np.arange(len(colours) + 1) - 0.5, len(colours)),
                      linewidth=LINE_WIDTHS["band_edge"])  # fmt: skip
        ax.invert_yaxis()
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
        labels = [self._display(label) for label in labels]
        shown = labels if len(rows) <= HEAT["label_rows"] else labels[: len(labels) - len(rows)]
        ax.set_yticks(np.arange(len(shown)) + 0.5, shown)
        ax.tick_params(axis="both", length=0)
        if len(cols) <= HEAT["label_columns"]:
            ax.set_xticks(np.arange(len(cols)) + 0.5, cols, rotation=BAR_LABEL_ROTATION)
        else:
            ax.set_xticks([])
        if self.row_stats is not None and len(rows) <= HEAT["label_rows"]:
            right = ax.secondary_yaxis("right")
            start = len(labels) - len(rows)
            right.set_yticks(np.arange(len(rows)) + start + 0.5,
                             [self._stat(r) for r in rows])  # fmt: skip
            right.tick_params(axis="y", length=0)
            right.set_ylabel(self.row_stats_label)
            right.spines["right"].set_visible(False)

    def interactive(self) -> Any:
        """Render the plotly figure in the static's layout: tracks, the gap row, then the rows.

        Every row of ``_grid`` is a band of a numeric axis with one tick, so a track stands one
        row above the cells as in the static and no label is dropped. A band is as many
        pitches tall as its own label has lines, so one long label raises its own row alone.
        The row statistic is a second axis at right, which automargin keeps inside the figure.
        Rows, legend and column title read as the static's; narrower than they need, the page
        swaps in ``layout.meta.narrow`` (``_narrow``), each wrapped to a phone's plot. Hover
        lines fit a phone's plot at every width; plotly thins the column labels to the width
        the plot is drawn at.

        Returns:
            A plotly Figure with the light template; the page switches themes.
        """
        import plotly.graph_objects as go  # noqa: PLC0415 - heavy; imported after the profile

        rows, cols = self.rows_in_order(), self.columns_in_order()
        words = self._words(phone=False)
        edges = self._edges(words)
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
        bare = {"showline": False, "ticks": "", "ticklabelstandoff": PLOTLY_TICK_PX,
                "zeroline": False, "showgrid": False, "automargin": True}  # fmt: skip
        fig.update_layout(
            template=plotly_templates()["white"], height=self._page_height(edges, words),
            showlegend=True,
            # At the div's top left, not the plot's, which a phone's labels push right; never
            # scrolled, as plotly scrolls a legend past half the figure's height.
            legend={"xref": "container", "yref": "container", "x": 0, "y": 1,
                    "yanchor": "top", "maxheight": 1, "traceorder": "grouped",
                    "tracegroupgap": 0},
            xaxis={**bare, "type": "category", "categoryorder": "array", "categoryarray": cols,
                   "title": {"text": "" if self.column_title else words.title},
                   "tickangle": -BAR_LABEL_ROTATION,
                   "showticklabels": len(cols) <= HEAT["label_columns"]},
            yaxis={**bare, "tickmode": "array", **self._yaxis(edges, words)},
        )  # fmt: skip
        if self._shows_stats():
            # plotly draws no axis that no trace is on, so an empty one holds it.
            fig.add_trace(go.Scatter(x=[None], y=[None], yaxis="y2", showlegend=False,
                                     hoverinfo="skip"))  # fmt: skip
            fig.update_layout(yaxis2={
                **bare, "matches": "y", "overlaying": "y", "side": "right", "tickmode": "array",
                "tickvals": _centres(edges)[first:], "ticktext": [self._stat(r) for r in rows],
                "title": {"text": self.row_stats_label},
            })  # fmt: skip
        if self.column_title:
            fig.add_annotation(self._title_note(words.title))
            fig.update_layout(margin_b=self._bottom(words.title))
        if (narrow := self._narrow(words)) is not None:
            fig.update_layout(meta={"narrow": narrow})
        return fig

    def _words(self, phone: bool) -> "_Words":
        """The static's row labels, legend entries and column title, or a phone's.

        On a phone each wraps to the room it has there, a row label to half the plot.
        """
        labels, rows = self._grid()[1], len(self.rows_in_order())
        shown = labels if rows <= HEAT["label_rows"] else labels[: len(labels) - rows]
        ticks = {i: "<br>".join(_fit(line, PHONE_PLOT_PX / 2) if phone else line
                                for line in self._display(label).split("\n"))
                 for i, label in enumerate(shown) if label}  # fmt: skip
        legend = [s.level for s in self.swatches]
        legend += [f"{t.name}: {level}" for t in self.tracks for level in t.colours]
        if phone:
            room = PHONE_PLOT_PX - PLOTLY_LEGEND_SYMBOL_PX
            legend = [_fit(name, room, SIZES["legend"]) for name in legend]
        return _Words(ticks, legend, self._title(ticks, legend, phone))

    def _title(self, ticks: Mapping[int, str], legend: list[str], phone: bool) -> str:
        """The column title as a layout prints it.

        ``column_title`` is right-aligned under the column labels, wrapped to the plot at the
        narrowest width its layout is drawn at. ``column_label`` is the axis title, which the
        page script wraps to the room under the plot; a phone's layout makes room for it
        wrapped to half the plot. Centred under a phone's cells, a wider one ran off the plot.
        """
        if not self.column_title:
            if not phone:
                return self.column_label
            return _fit(self.column_label, PHONE_PLOT_PX / 2, SIZES["axis_label"])
        width = PHONE_PLOT_PX if phone else self._below(ticks, legend)
        room = width - PLOTLY_MARGIN_PX["l"] - PLOTLY_MARGIN_PX["r"] - self._right_px()
        return _fit(self.column_title, room, SIZES["axis_label"])

    def _edges(self, words: "_Words") -> list[float]:
        """Band edges top down, in row pitches: each band as many pitches as its label's lines."""
        heights = [_lines(words.ticks.get(i, "")) for i in range(len(self._grid()[1]))]
        return [0.0, *itertools.accumulate(map(float, heights))]

    def _yaxis(self, edges: list[float], words: "_Words") -> dict[str, Any]:
        """The row axis's range and ticks: one tick at the middle of each labelled band."""
        centres = _centres(edges)
        return {"range": [edges[-1], edges[0]], "tickvals": [centres[i] for i in words.ticks],
                "ticktext": list(words.ticks.values())}  # fmt: skip

    def _narrow(self, wide: "_Words") -> dict[str, Any] | None:
        """The layout the page script swaps in while the plot is narrower than ``below`` px.

        At 375 px the band map's two-line labels started up to 132 px left of the plot and lost
        their unit (2026-10-01); wrapped for a phone on every screen, a desktop's rows were as
        tall as a phone's labels (work order item 2, 2026-10-01). Narrower than the static's
        labels need, each wraps to a phone's plot and each row grows to its label's lines: new
        band edges for the heatmap traces, names for the legend's (``restyle``), ticks, title
        and height (``relayout``). None where nothing wraps.
        """
        words = self._words(phone=True)
        if words == wide:
            return None
        edges = self._edges(words)
        first = len(edges) - 1 - len(self.rows_in_order())
        relayout = {f"yaxis.{key}": value for key, value in self._yaxis(edges, words).items()}
        if self._shows_stats():
            relayout["yaxis2.tickvals"] = _centres(edges)[first:]
        if self.column_title:
            relayout |= {"annotations[0].text": words.title, "margin.b": self._bottom(words.title)}
        relayout["height"] = self._page_height(edges, words)
        heat = list(range(len(self.tracks) + 1))
        y = [edges[first:], *(edges[at : at + 2] for at in range(len(self.tracks)))]
        keys = list(range(len(heat), len(heat) + len(words.legend)))
        return {"below": self._below(wide.ticks, wide.legend), "relayout": relayout,
                "restyle": [[{"y": y}, heat], [{"name": words.legend}, keys]]}  # fmt: skip

    def _below(self, ticks: Mapping[int, str], legend: list[str]) -> float:
        """The narrowest plot the static's labels fit, in px.

        The widest row label beside a line's width for each labelled column, upright, and the
        statistic's axis; or the widest legend entry beside its swatch.
        """
        tick = on_page(SIZES["tick"])
        widest = max((text_width(line, tick) for t in ticks.values() for line in t.split("<br>")),
                     default=0.0)  # fmt: skip
        columns = min(len(self.columns_in_order()), HEAT["label_columns"])
        plot = (PLOTLY_MARGIN_PX["l"] + widest + PLOTLY_TICK_PX + columns * tick * LINE_SPACING
                + self._right_px() + PLOTLY_MARGIN_PX["r"])  # fmt: skip
        size = on_page(SIZES["legend"])
        entries = max((text_width(name, size) for name in legend), default=0.0)
        return max(plot, entries + PLOTLY_LEGEND_SYMBOL_PX)

    def _shows_stats(self) -> bool:
        """Whether the page prints the row statistic, as the static does, on a right axis."""
        return self.row_stats is not None and len(self.rows_in_order()) <= HEAT["label_rows"]

    def _right_px(self) -> float:
        """What the statistic's axis takes at right, in px; nothing where none is printed.

        Its widest label, the label's and the title's standoffs, and the title's line.
        """
        if not self._shows_stats():
            return 0.0
        tick, title = on_page(SIZES["tick"]), on_page(SIZES["axis_label"])
        widest = max((text_width(self._stat(r), tick) for r in self.rows_in_order()), default=0.0)
        return widest + PLOTLY_TICK_PX + PLOTLY_TITLE_STANDOFF_PX + title * LINE_SPACING

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
        clipped (Rytis 2026-10-01).
        """
        return {"text": text, "xref": "paper", "x": 1, "xanchor": "right", "yref": "paper",
                "y": 0, "yanchor": "top", "yshift": -self._title_drop(), "align": "right",
                "showarrow": False, "font": {"size": on_page(SIZES["axis_label"])}}  # fmt: skip

    def _bottom(self, title: str) -> float:
        """The page plot's height under its cells: column labels, the column title, margins.

        An axis title takes plotly's own room under the labels, which automargin grows to;
        a right-aligned title sits ``_title_drop`` under the cells, in a bottom margin set to
        hold it.
        """
        title_px = _lines(title) * on_page(SIZES["axis_label"]) * LINE_SPACING
        if self.column_title:
            return PLOTLY_MARGIN_PX["b"] + self._title_drop() + title_px
        return PLOTLY_HEAT_BASE_PX + self._column_label_length(on_page(SIZES["tick"])) + title_px

    def _page_height(self, edges: list[float], words: "_Words") -> float:
        """The page figure's height in px: the static's rows at the page pitch, then the rest.

        As ``_inches`` does, the rows drawn are capped and the column labels measured. The
        legend has a line per entry, the most a phone stacks (M7's nine on nine lines at
        375 px, 2026-10-01); a wider page gives the rows what it leaves.
        """
        rows = len(self.rows_in_order())
        above = edges[len(edges) - 1 - rows]  # the tracks and the blank row
        cells = edges[-1] - above if rows <= HEAT["max_rows_drawn"] else HEAT["max_rows_drawn"]
        lines = sum(_lines(name) for name in words.legend)
        top = max(PLOTLY_MARGIN_PX["t"], lines * PLOTLY_LEGEND_LINE_PX + PLOTLY_LEGEND_PAD_PX)
        return top + (above + cells) * PLOTLY_HEAT_ROW_PX + self._bottom(words.title)

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
            z=z, x=cols, y=edges, text=text, showscale=False, xgap=0, ygap=0,
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
            z=[[codes[track.values[c]] for c in cols]], x=cols, y=edges, text=[hover],
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


@dataclass(frozen=True)
class _Words:
    """What one page layout prints, each text with ``<br>`` where it breaks.

    Row labels by grid row, legend entries and the column title.
    """

    ticks: Mapping[int, str]
    legend: list[str]
    title: str


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


def _fit(text: str, width_px: float, size: float = SIZES["tick"]) -> str:
    """Text wrapped at spaces to lines of at most ``width_px`` on the page, in plotly's ``<br>``.

    Measured as ``ctx.viz.metrics`` measures every wrapped text. A word longer than a line
    stays whole: split, it would read as two words.
    """
    return "<br>".join(wrap_to(text, width_px, on_page(size))) or text


def _lines(text: str) -> int:
    """How many lines a plotly text runs to."""
    return text.count("<br>") + 1
