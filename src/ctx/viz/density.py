"""Density contours: where most samples of each group lie on two metrics, and who stands out.

Each group of at least ``DENSITY_MIN_POINTS`` samples is drawn as its 50 % (solid) and 90 %
(dashed) highest-density regions, from one smoothed histogram on a grid every group shares
(``ctx.viz.contour``); smaller groups are drawn as points, and so are the highlighted samples.
``static()`` and ``interactive()`` read the same ``DensityData`` and draw the same levels.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from functools import cached_property
from typing import Any

import numpy as np
import polars as pl

from ctx.viz.annotate import files_room, fit_title, place_caption, refit_y_title
from ctx.viz.axes import axis_title, hover_field, hover_rows, page_y_title
from ctx.viz.contour import (
    DENSITY_BINS,
    DENSITY_MARGIN,
    DENSITY_MIN_POINTS,
    DENSITY_SIGMA,
    HDR_MASSES,
    density_grid,
    hdr_levels,
)
from ctx.viz.scatter import log_ticks
from ctx.viz.style import (
    DENSITY_DASH,
    FIGURE_WIDTH,
    LINE_WIDTHS,
    MARKER_SIZE,
    MARKERS,
    MAX_CATEGORIES,
    PLOTLY_HEIGHT_PX,
    PLOTLY_MARKER_PX,
    POINT_ALPHA,
    Theme,
    apply_matplotlib,
    figure_inches,
    plotly_templates,
)

REGIONS = ("inner", "outer", "outside")


@dataclass(frozen=True)
class Density:
    """One group's density on the shared grid and its two highest-density levels."""

    group: str
    n: int
    density: np.ndarray
    levels: tuple[float, float]


@dataclass(frozen=True)
class DensityData:
    """The contract both renderers read.

    Attributes:
        frame: One row per sample.
        columns: (unit, meaning) for every column of ``frame``.
        x: Column on the x axis.
        y: Column on the y axis.
        group: Column whose levels each get a density, in their colour.
        x_label: x-axis title, without the unit.
        y_label: y-axis title, without the unit.
        colours: Group level to hex colour, in legend order.
        caption: The caption, unwrapped; it states n per group and the bandwidth.
        highlight: Boolean column of samples drawn as points over the contours.
        highlight_label: Legend text of those points.
        region_labels: Legend text of the inner and outer regions, e.g. "50 %" and "90 %".
        hover: Columns shown on hover of a drawn point.
        log_x: Log-scaled x axis; the density is computed on log10 x.
        log_y: Log-scaled y axis; likewise.
        width: The column; every figure takes ``FIGURE_WIDTH``, the double column.
    """

    frame: pl.DataFrame
    columns: Mapping[str, tuple[str, str]]
    x: str
    y: str
    group: str
    x_label: str
    y_label: str
    colours: Mapping[str, str]
    caption: str
    highlight: str | None = None
    highlight_label: str = ""
    region_labels: tuple[str, str] = ("50 %", "90 %")
    hover: tuple[str, ...] = ()
    log_x: bool = False
    log_y: bool = False
    width: str = FIGURE_WIDTH

    def __post_init__(self) -> None:
        """Refuse a contract the renderers could only draw wrongly."""
        named = [self.x, self.y, self.group, self.highlight, *self.hover]
        if missing := [c for c in named if c is not None and c not in self.frame.columns]:
            raise ValueError(f"missing columns {missing}")
        if len(self.colours) > MAX_CATEGORIES:
            raise ValueError(f"{len(self.colours)} groups; at most {MAX_CATEGORIES}")
        if stray := set(self.frame[self.group].unique()) - set(self.colours):
            raise ValueError(f"{self.group}: levels {sorted(stray)} have no colour")
        for column, log in ((self.x, self.log_x), (self.y, self.log_y)):
            if log and not (self.frame[column] > 0).all():
                raise ValueError(f"{column}: a log axis needs positive values")

    def _scaled(self, column: str, log: bool) -> np.ndarray:
        """A column on the scale its axis draws: log10 for a log axis."""
        values = self.frame[column].to_numpy().astype(float)
        return np.log10(values) if log else values

    @cached_property
    def extent(self) -> tuple[float, float, float, float]:
        """The shared grid's bounds, on the axes' scales, a margin beyond every point."""
        bounds = []
        for column, log in ((self.x, self.log_x), (self.y, self.log_y)):
            v = self._scaled(column, log)
            pad = DENSITY_MARGIN * (v.max() - v.min() or 1.0)
            bounds += [float(v.min() - pad), float(v.max() + pad)]
        return bounds[0], bounds[1], bounds[2], bounds[3]

    @cached_property
    def densities(self) -> dict[str, Density]:
        """The groups with enough samples, each with its density and levels."""
        x, y = self._scaled(self.x, self.log_x), self._scaled(self.y, self.log_y)
        groups = self.frame[self.group].to_numpy()
        out = {}
        for level in self.colours:
            keep = groups == level
            if keep.sum() >= DENSITY_MIN_POINTS:
                _, _, d = density_grid(x[keep], y[keep], self.extent)
                inner, outer = hdr_levels(d)
                out[level] = Density(level, int(keep.sum()), d, (inner, outer))
        return out

    def centres(self) -> tuple[np.ndarray, np.ndarray]:
        """The grid's bin centres in data units, as both renderers place the contours."""
        x0, x1, y0, y1 = self.extent
        xs = x0 + (np.arange(DENSITY_BINS) + 0.5) * (x1 - x0) / DENSITY_BINS
        ys = y0 + (np.arange(DENSITY_BINS) + 0.5) * (y1 - y0) / DENSITY_BINS
        return (10**xs if self.log_x else xs), (10**ys if self.log_y else ys)

    def bandwidth(self) -> str:
        """The smoothing, for the caption: sigma in bins and in axis units."""
        x0, x1, y0, y1 = self.extent
        scale = ("log10 " if self.log_x else "", "log10 " if self.log_y else "")
        return (f"Gaussian sigma {DENSITY_SIGMA:g} of {DENSITY_BINS} bins: "
                f"{DENSITY_SIGMA * (x1 - x0) / DENSITY_BINS:.3g} {scale[0]}x units, "
                f"{DENSITY_SIGMA * (y1 - y0) / DENSITY_BINS:.3g} {scale[1]}y units")  # fmt: skip

    def _regions(self) -> list[str | None]:
        """Each sample's region in its own group's density; None for an undrawn group."""
        x0, x1, y0, y1 = self.extent
        x, y = self._scaled(self.x, self.log_x), self._scaled(self.y, self.log_y)
        col = np.clip(((x - x0) / (x1 - x0) * DENSITY_BINS).astype(int), 0, DENSITY_BINS - 1)
        row = np.clip(((y - y0) / (y1 - y0) * DENSITY_BINS).astype(int), 0, DENSITY_BINS - 1)
        out: list[str | None] = []
        for i, level in enumerate(self.frame[self.group].to_list()):
            d = self.densities.get(level)
            if d is None:
                out.append(None)
                continue
            value = float(d.density[row[i], col[i]])
            inner, outer = d.levels
            out.append(REGIONS[0] if value >= inner else REGIONS[1] if value >= outer
                       else REGIONS[2])  # fmt: skip
        return out

    def table(self) -> pl.DataFrame:
        """Every sample, its group's n, and the region of its group's density it falls in."""
        n = {g: d.n for g, d in self.densities.items()}
        return self.frame.with_columns(
            group_n=pl.col(self.group).replace_strict(n, default=None, return_dtype=pl.Int32),
            region=pl.Series(self._regions(), dtype=pl.String),
        )

    def table_columns(self) -> Mapping[str, tuple[str, str]]:
        """(unit, meaning) for every column of ``table()``."""
        inner, outer = (f"{m:.0%}" for m in HDR_MASSES)
        return {**self.columns,
                "group_n": ("samples", "Samples in the group; empty when too few for a density."),
                "region": ("none", f"inner: in the group's {inner} highest-density region; "
                           f"outer: in its {outer} region only; outside: in neither.")}  # fmt: skip

    def page_table(self) -> pl.DataFrame:
        """Every sample, with its region: the page draws some as points, all as contours."""
        return self.table()

    def caption_text(self) -> str:
        """The caption, unwrapped."""
        return self.caption

    def _drawn(self) -> pl.Expr:
        """Samples drawn as points: highlighted, or in a group too small for a density."""
        small = ~pl.col(self.group).is_in(list(self.densities))
        return small | pl.col(self.highlight) if self.highlight else small

    def _titles(self) -> tuple[str, str]:
        """The axis titles with their units."""
        return (axis_title(self.x_label, self.columns[self.x][0]),
                axis_title(self.y_label, self.columns[self.y][0]))  # fmt: skip

    def _inches(self) -> tuple[float, float]:
        """The figure files' size, whose height also sets where the page breaks its y title."""
        return figure_inches(self.width, "contour")

    def static(self, theme: Theme) -> Any:
        """Render the matplotlib figure for one theme.

        Args:
            theme: One of ``style.THEMES``.

        Returns:
            A matplotlib Figure, not attached to pyplot.
        """
        import seaborn as sns  # noqa: PLC0415 - heavy
        from matplotlib.figure import Figure  # noqa: PLC0415 - heavy

        apply_matplotlib(theme)
        width, height = self._inches()
        fig = Figure(figsize=(width, height))
        ax = fig.subplots()
        xs, ys = self.centres()
        for group, d in self.densities.items():
            ax.contour(xs, ys, d.density, levels=sorted(d.levels), colors=self.colours[group],
                       linewidths=LINE_WIDTHS["line"],
                       linestyles=[DENSITY_DASH["outer"][0], DENSITY_DASH["inner"][0]])  # fmt: skip
        points = self.frame.filter(self._drawn())
        for group, colour in self.colours.items():
            rows = points.filter(pl.col(self.group) == group)
            flagged = rows[self.highlight].to_list() if self.highlight else [False] * rows.height
            for flag in (False, True):
                keep = [i for i, f in enumerate(flagged) if f == flag]
                ax.scatter(rows[self.x].to_numpy()[keep], rows[self.y].to_numpy()[keep],
                           s=MARKER_SIZE**2, c=colour, marker=MARKERS[int(flag)][0],
                           alpha=POINT_ALPHA, linewidths=LINE_WIDTHS["marker_edge"])  # fmt: skip
        self._axes(ax, height)
        sns.despine(fig=fig)
        place_caption(fig, self.caption, width, height, ax=ax, handles=self._handles(theme))
        refit_y_title(fig, ax, self._titles()[1])
        return fig

    def _axes(self, ax: Any, height: float) -> None:
        """Scales, plain log labels and titles."""
        ax.set_xscale("log" if self.log_x else "linear")
        ax.set_yscale("log" if self.log_y else "linear")
        for axis, log, column in ((ax.xaxis, self.log_x, self.x), (ax.yaxis, self.log_y, self.y)):
            if log:
                log_ticks(axis, self.frame[column])
        x_title, y_title = self._titles()
        ax.set_xlabel(x_title)
        ax.set_ylabel(fit_title(y_title, height))

    def _handles(self, theme: Theme) -> list[Any]:
        """Legend entries: each group, the two regions, the highlight triangle."""
        from matplotlib.lines import Line2D  # noqa: PLC0415 - heavy

        present = set(self.frame[self.group].unique())
        handles = [Line2D([], [], color=c, label=str(g))
                   for g, c in self.colours.items() if g in present]  # fmt: skip
        regions = zip(("inner", "outer"), self.region_labels, strict=True)
        handles += [Line2D([], [], color=theme.ink, linestyle=DENSITY_DASH[key][0], label=text)
                    for key, text in regions]  # fmt: skip
        if self.highlight and self.frame[self.highlight].any():
            handles.append(Line2D([], [], marker=MARKERS[1][0], linestyle="", color=theme.ink,
                                  label=self.highlight_label))  # fmt: skip
        return handles

    def interactive(self) -> Any:
        """Render the plotly figure: the same levels on the same grid, points with hover.

        Returns:
            A plotly Figure with the light template; the page switches themes.
        """
        import plotly.graph_objects as go  # noqa: PLC0415 - heavy; imported after the profile

        xs, ys = self.centres()
        fig = go.Figure()
        for group, d in self.densities.items():
            for key, level in zip(("inner", "outer"), d.levels, strict=True):
                fig.add_trace(go.Contour(
                    x=xs, y=ys, z=d.density, autocontour=False, showscale=False,
                    contours={"start": level, "end": level, "size": 1, "coloring": "none"},
                    line={"color": self.colours[group], "width": LINE_WIDTHS["line"],
                          "dash": DENSITY_DASH[key][1]},
                    name=f"{group} {self.region_labels[key == 'outer']}", hoverinfo="skip",
                ))  # fmt: skip
        x_title, y_title = self._titles()
        schema = self.frame.schema
        template = "<br>".join([
            f"{x_title}: {hover_field('x', schema[self.x])}",
            f"{y_title}: {hover_field('y', schema[self.y])}",
            *(f"{c}: {hover_field(f'customdata[{i}]', schema[c])}"
              for i, c in enumerate(self.hover)),
        ])  # fmt: skip
        points = self.frame.filter(self._drawn())
        for group, colour in self.colours.items():
            rows = points.filter(pl.col(self.group) == group)
            if rows.height:
                flags = rows[self.highlight].to_list() if self.highlight else [False] * rows.height
                fig.add_trace(go.Scatter(
                    x=rows[self.x].to_list(), y=rows[self.y].to_list(), mode="markers",
                    name=str(group),
                    customdata=hover_rows(rows, self.hover) if self.hover else None,
                    hovertemplate=template + "<extra></extra>",
                    marker={"color": colour, "size": PLOTLY_MARKER_PX, "opacity": POINT_ALPHA,
                            "symbol": [MARKERS[int(f)][1] for f in flags]},
                ))  # fmt: skip
        fig.update_layout(
            template=plotly_templates()["white"], height=PLOTLY_HEIGHT_PX["contour"],
            xaxis={"title": {"text": x_title}, "type": "log" if self.log_x else "linear"},
            yaxis={"title": {"text": page_y_title(y_title, files_room(self.static))},
                   "type": "log" if self.log_y else "linear"},
        )  # fmt: skip
        return fig
