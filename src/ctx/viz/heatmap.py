"""Heatmap: every cell holds a level, every level has a swatch; vector cells, a swatch legend.

Two uses share the contract. A QC band map colours samples × metrics good, usable, bad or no
band. An expression map bins z-scores into a stepped diverging scale, two data colours at
three opacities (``diverging_levels``), orders rows and columns by average-linkage clustering
on 1 − Pearson (``cluster_order``), draws category tracks above the columns and a statistic
beside each row. ``static()`` and ``interactive()`` read the same ``HeatmapData``.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import polars as pl

from ctx.viz.annotate import place_caption
from ctx.viz.style import (
    BAR_LABEL_ROTATION,
    HEAT,
    HEAT_ALPHAS,
    LINE_WIDTHS,
    MM_PER_INCH,
    PLOTLY_HEAT_BASE_PX,
    PLOTLY_HEAT_ROW_PX,
    ROLES,
    THEMES,
    WIDTHS_MM,
    Theme,
    apply_matplotlib,
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
        """Width by the column rule; height grows with the rows drawn, capped."""
        drawn = min(len(self.rows_in_order()), HEAT["max_rows_drawn"]) + len(self.tracks) + GAP
        return (WIDTHS_MM[self.width] / MM_PER_INCH,
                (drawn * HEAT["row_mm"] + HEAT["margin_mm"]) / MM_PER_INCH)  # fmt: skip

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
        """Render the plotly figure: one heatmap trace per block, swatches as legend entries.

        Returns:
            A plotly Figure with the light template; the page switches themes.
        """
        import plotly.graph_objects as go  # noqa: PLC0415 - heavy; imported after the profile

        light = THEMES["white"]
        rows, cols = self.rows_in_order(), self.columns_in_order()
        fig = go.Figure()
        index = {s.level: i for i, s in enumerate(self.swatches)}
        cells = self.table()
        lookup = {(r, c): (lv, v) for r, c, lv, v in
                  cells.select(self.row, self.column, self.level, self.value).rows()}  # fmt: skip
        z = [[index.get(lookup.get((r, c), (None, None))[0]) for c in cols] for r in rows]
        text = [[f"{lookup.get((r, c), (None, None))[1]}" for c in cols] for r in rows]
        fig.add_trace(go.Heatmap(
            z=z, x=cols, y=rows, text=text, showscale=False, xgap=0, ygap=0,
            colorscale=plotly_steps([tint(s.colour or light.grid, s.alpha) for s in self.swatches]),
            zmin=-0.5, zmax=len(self.swatches) - 0.5,
            hovertemplate=f"{self.row_label}: %{{y}}<br>{self.column_label}: %{{x}}<br>"
            "%{text}<extra></extra>",
        ))  # fmt: skip
        for track in self.tracks:
            codes = {level: i for i, level in enumerate(track.colours)}
            fig.add_trace(go.Heatmap(
                z=[[codes[track.values[c]] for c in cols]], x=cols, y=[track.name],
                text=[[track.values[c] for c in cols]], showscale=False,
                colorscale=plotly_steps(list(track.colours.values())), zmin=-0.5,
                zmax=len(codes) - 0.5, hovertemplate=f"{track.name}: %{{text}}<extra></extra>",
            ))  # fmt: skip
        for s in self.swatches:
            fig.add_trace(
                go.Scatter(
                    x=[None],
                    y=[None],
                    mode="markers",
                    name=s.level,
                    marker={"symbol": "square", "color": tint(s.colour or light.grid, s.alpha)},
                )
            )
        order = [*reversed(rows), *reversed([t.name for t in self.tracks])]
        fig.update_layout(
            template=plotly_templates()["white"],
            height=PLOTLY_HEAT_BASE_PX + PLOTLY_HEAT_ROW_PX * (len(rows) + len(self.tracks)),
            xaxis={"title": {"text": self.column_label}, "showline": False,
                   "showticklabels": len(cols) <= HEAT["label_columns"]},
            yaxis={"categoryorder": "array", "categoryarray": order, "showline": False,
                   "showticklabels": len(rows) <= HEAT["label_rows"]},
            annotations=self._stat_notes(rows),
        )  # fmt: skip
        return fig

    def _stat(self, row: str) -> str:
        """A row's statistic as the figure prints it; empty where it has none."""
        value = (self.row_stats or {}).get(row)
        return "" if value is None else f"{value:{self.row_stats_format}}"

    def _stat_notes(self, rows: list[str]) -> list[dict[str, Any]]:
        """The row statistics beside the plot, one annotation per row."""
        if self.row_stats is None or len(rows) > HEAT["label_rows"]:
            return []
        return [{"x": 1.0, "xref": "paper", "xanchor": "left", "y": r, "yref": "y",
                 "text": self._stat(r), "showarrow": False} for r in rows]  # fmt: skip
