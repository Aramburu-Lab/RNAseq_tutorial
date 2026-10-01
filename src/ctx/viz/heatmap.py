"""Heatmap: every cell holds a level, every level has a swatch; vector cells, a swatch legend.

Two uses share the contract. A QC band map colours samples × metrics good, usable, bad or no
band. An expression map bins z-scores into a stepped diverging scale, two data colours at
three opacities (``diverging_levels``), orders rows and columns by average-linkage clustering
on 1 − Pearson (``cluster_order``), draws category tracks above the columns and a statistic
beside each row. ``static()`` and ``interactive()`` read the same ``HeatmapData``.
"""

import itertools
import textwrap
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import polars as pl

from ctx.viz.annotate import place_caption
from ctx.viz.axes import axis_title
from ctx.viz.style import (
    AVERAGE_CHAR_EM,
    BAR_LABEL_ROTATION,
    HEAT,
    HEAT_ALPHAS,
    LINE_SPACING,
    LINE_WIDTHS,
    MM_PER_INCH,
    PLOTLY_HEAT_BASE_PX,
    PLOTLY_HEAT_ROW_PX,
    PLOTLY_HOVER_PX,
    PLOTLY_LEGEND_LINE_PX,
    PLOTLY_LEGEND_PAD_PX,
    PLOTLY_LEGEND_SWATCH_PX,
    PLOTLY_LEGEND_SYMBOL_PX,
    PLOTLY_MARGIN_PX,
    PLOTLY_NARROW_PX,
    PLOTLY_NO_LEVEL,
    PLOTLY_TICK_PX,
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
        column_label: What a column is.
        swatches: Every level's colour, in legend order.
        caption: The caption, unwrapped.
        row_order: Rows top to bottom; empty keeps first appearance.
        column_order: Columns left to right; empty keeps first appearance.
        tracks: Categories drawn above the columns.
        row_stats: Row to the number printed at its right, e.g. an adjusted p.
        row_stats_label: What that number is.
        row_stats_format: Its format spec on the figure; the CSV keeps the number.
        width: ``"single"`` or ``"double"`` column.
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

    def _inches(self) -> tuple[float, float]:
        """Width by the column rule; height by the rows drawn and the longest column label.

        The rows drawn are capped; a column label runs upward and, unmeasured, squeezed the rows.
        """
        drawn = min(len(self.rows_in_order()), HEAT["max_rows_drawn"]) + len(self.tracks) + GAP
        columns = self.columns_in_order()
        labels = max(map(len, columns)) if len(columns) <= HEAT["label_columns"] else 0
        label_mm = labels * SIZES["tick"] * AVERAGE_CHAR_EM / 72 * MM_PER_INCH
        return (WIDTHS_MM[self.width] / MM_PER_INCH,
                (drawn * HEAT["row_mm"] + HEAT["margin_mm"] + label_mm) / MM_PER_INCH)  # fmt: skip

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
        shown = labels if len(rows) <= HEAT["label_rows"] else labels[: len(labels) - len(rows)]
        ax.set_yticks(np.arange(len(shown)) + 0.5, shown)
        ax.tick_params(axis="both", length=0)
        if len(cols) <= HEAT["label_columns"]:
            ax.set_xticks(np.arange(len(cols)) + 0.5, cols, rotation=BAR_LABEL_ROTATION)
        else:
            ax.set_xticks([])
        ax.set_xlabel(self.column_label)
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
        Labels, the column title, legend entries and hover lines wrap to fit a phone's plot
        (``_fit``); plotly thins the column labels to the width the plot is drawn at.

        Returns:
            A plotly Figure with the light template; the page switches themes.
        """
        import plotly.graph_objects as go  # noqa: PLC0415 - heavy; imported after the profile

        rows, cols = self.rows_in_order(), self.columns_in_order()
        labels = self._grid()[1]
        ticks = self._ticks(labels)
        edges = _edges([_lines(ticks.get(i, "")) for i in range(len(labels))])
        centres = [(top + bottom) / 2 for top, bottom in itertools.pairwise(edges)]
        first = len(labels) - len(rows)
        fig = go.Figure(self._cell_trace(go, rows, cols, edges[first:]))
        for at, track in enumerate(self.tracks):
            fig.add_trace(self._track_trace(go, track, cols, edges[at : at + 2]))
        # One legend group per entry: plotly lays ungrouped entries out at the widest one's
        # width, one a line on an overview's plot (2026-10-01); groups flow by their own widths.
        for i, (name, colour) in enumerate(self._legend()):
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode="markers", name=name, legendgroup=str(i),
                marker={"symbol": "square", "color": colour, "size": PLOTLY_LEGEND_SWATCH_PX},
            ))  # fmt: skip
        bare = {"showline": False, "ticks": "", "ticklabelstandoff": PLOTLY_TICK_PX,
                "zeroline": False, "showgrid": False, "automargin": True}  # fmt: skip
        fig.update_layout(
            template=plotly_templates()["white"], height=self._page_height(edges),
            showlegend=True,
            # At the div's top left, not the plot's, which a phone's labels push right; never
            # scrolled, as plotly scrolls a legend past half the figure's height.
            legend={"xref": "container", "yref": "container", "x": 0, "y": 1,
                    "yanchor": "top", "maxheight": 1, "traceorder": "grouped",
                    "tracegroupgap": 0},
            xaxis={**bare, "type": "category", "categoryorder": "array", "categoryarray": cols,
                   "title": {"text": self._title()}, "tickangle": -BAR_LABEL_ROTATION,
                   "showticklabels": len(cols) <= HEAT["label_columns"]},
            yaxis={**bare, "range": [edges[-1], edges[0]], "tickmode": "array",
                   "tickvals": [centres[i] for i in ticks], "ticktext": list(ticks.values())},
        )  # fmt: skip
        if self.row_stats is not None and len(rows) <= HEAT["label_rows"]:
            # plotly draws no axis that no trace is on, so an empty one holds it.
            fig.add_trace(go.Scatter(x=[None], y=[None], yaxis="y2", showlegend=False,
                                     hoverinfo="skip"))  # fmt: skip
            fig.update_layout(yaxis2={
                **bare, "matches": "y", "overlaying": "y", "side": "right", "tickmode": "array",
                "tickvals": centres[first:], "ticktext": [self._stat(r) for r in rows],
                "title": {"text": self.row_stats_label},
            })  # fmt: skip
        return fig

    def _ticks(self, labels: list[str]) -> dict[int, str]:
        """Grid row to label, for every row the static labels, wrapped to a phone's half."""
        rows = len(self.rows_in_order())
        shown = labels if rows <= HEAT["label_rows"] else labels[: len(labels) - rows]
        return {i: _fit(label, PLOTLY_NARROW_PX / 2) for i, label in enumerate(shown) if label}

    def _title(self) -> str:
        """The column axis title, wrapped as wide as a row label.

        Centred under a phone's cells, a wider title runs off the plot's right edge.
        """
        return _fit(self.column_label, PLOTLY_NARROW_PX / 2, SIZES["axis_label"])

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
        """A cell's hover: its row and column, its number named with its unit, then its level."""
        lines = [_named(self.row_label, row), _named(self.column_label, column)]
        if value is not None and not np.isnan(value):
            unit = self.columns[self.value][0]
            lines.append(_named(axis_title(self.value, unit), f"{value:.3g}"))
        return _hover_lines([*lines, *([level] if level is not None else [])])

    def _legend(self) -> list[tuple[str, str]]:
        """The static legend's entries, wrapped to a phone's width, with their page colours."""
        swatches = [(s.level, _page_colour(s)) for s in self.swatches]
        tracks = [(f"{t.name}: {level}", colour)
                  for t in self.tracks for level, colour in t.colours.items()]  # fmt: skip
        width = PLOTLY_NARROW_PX - PLOTLY_LEGEND_SYMBOL_PX
        return [(_fit(name, width, SIZES["legend"]), c) for name, c in [*swatches, *tracks]]

    def _page_height(self, edges: list[float]) -> float:
        """The page figure's height in px: the static's rows at the page pitch, then the rest.

        As ``_inches`` does, the rows drawn are capped and the column labels measured. The
        legend has a line per entry, the most a phone stacks (M7's nine on nine lines at
        375 px, 2026-10-01); a wider page gives the rows what it leaves.
        """
        rows = len(self.rows_in_order())
        above = edges[len(edges) - 1 - rows]  # the tracks and the blank row
        cells = edges[-1] - above if rows <= HEAT["max_rows_drawn"] else HEAT["max_rows_drawn"]
        columns = self.columns_in_order()
        shown = len(columns) <= HEAT["label_columns"]
        label_px = max(map(len, columns)) * on_page(SIZES["tick"]) * AVERAGE_CHAR_EM if shown else 0
        title_px = _lines(self._title()) * on_page(SIZES["axis_label"]) * LINE_SPACING
        lines = sum(_lines(n) for n, _ in self._legend())
        top = max(PLOTLY_MARGIN_PX["t"], lines * PLOTLY_LEGEND_LINE_PX + PLOTLY_LEGEND_PAD_PX)
        rows_px = (above + cells) * PLOTLY_HEAT_ROW_PX
        return PLOTLY_HEAT_BASE_PX + top + label_px + title_px + rows_px

    def _stat(self, row: str) -> str:
        """A row's statistic as the figure prints it; empty where it has none."""
        value = (self.row_stats or {}).get(row)
        return "" if value is None else f"{value:{self.row_stats_format}}"


def _page_colour(swatch: Swatch) -> str:
    """A swatch's colour on the page; "nothing to show" is a grey that reads in both themes."""
    colour, alpha = (swatch.colour, 1.0) if swatch.colour else PLOTLY_NO_LEVEL
    return tint(colour, alpha * swatch.alpha)


def _edges(heights: list[int]) -> list[float]:
    """Band edges top down, in row pitches, for bands this many pitches tall."""
    return [0.0, *itertools.accumulate(map(float, heights))]


def _named(label: str, value: str) -> str:
    """``label: value``; a label written as a sentence loses its full stop first."""
    return f"{label.removesuffix('.')}: {value}"


def _hover_lines(lines: list[str]) -> str:
    """Hover lines, each wrapped to a phone's plot, which a wider hover box ran off."""
    return "<br>".join(_fit(line, PLOTLY_HOVER_PX) for line in lines)


def _fit(text: str, width_px: float, size: float = SIZES["tick"]) -> str:
    """Text wrapped at spaces to lines of at most ``width_px`` on the page, in plotly's ``<br>``.

    A word longer than a line stays whole: split, it would read as two words.
    """
    chars = int(width_px / (on_page(size) * AVERAGE_CHAR_EM))
    lines = textwrap.wrap(text, chars, break_on_hyphens=False, break_long_words=False)
    return "<br>".join(lines) or text


def _lines(text: str) -> int:
    """How many lines a plotly text runs to."""
    return text.count("<br>") + 1
