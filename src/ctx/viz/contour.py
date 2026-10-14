"""Contours of a continuous value over two continuous axes: computed surfaces and densities.

``SurfaceData`` fills a value computed on a grid (a two-knob map, a power map) in eight
cividis bands, with labelled isolines and the current setting marked. Cells at an optional
floor, such as settings that keep no read, are drawn in the reference grey and the bands span
the cells above it: spanning a fastp grid's 0 % put its 30 other cells in the top two bands.
The density helpers, ``density_grid`` and ``hdr_levels``, compute once the smoothed histogram
and the 50 % and 90 % highest-density levels that ``ctx.viz.density`` draws in both renderers.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import polars as pl

from ctx.viz.annotate import fit_title, place_caption, refit_y_title
from ctx.viz.axes import axis_title
from ctx.viz.style import (
    CIVIDIS,
    FIGURE_WIDTH,
    INK_COLOUR,
    LABEL_BOX_PAD,
    LINE_WIDTHS,
    MARKER_SIZE,
    PLOTLY_HEIGHT_PX,
    PLOTLY_LEGEND_SWATCH_PX,
    PLOTLY_MARKER_PX,
    SIZES,
    THEMES,
    Theme,
    apply_matplotlib,
    figure_inches,
    plotly_steps,
    plotly_templates,
)

# The smoothed 2D histogram both renderers draw a density from: bins per axis, the Gaussian
# kernel's sigma in bins, the margin around the points as a share of their range, and the
# fewest points a group needs for its density to be drawn (CONVENTIONS, Figures).
DENSITY_BINS, DENSITY_SIGMA, DENSITY_MARGIN, DENSITY_MIN_POINTS = 64, 2.0, 0.1, 30
HDR_MASSES = (0.5, 0.9)


def density_grid(
    x: np.ndarray, y: np.ndarray, extent: tuple[float, float, float, float]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Bin centres and the smoothed density of points on a shared grid, summing to 1.

    Args:
        x: Point x values, on the scale the axis draws (log10 for a log axis).
        y: Point y values, likewise.
        extent: (x low, x high, y low, y high) of the grid every group shares.

    Returns:
        x centres, y centres, and the density with rows along y.
    """
    from scipy.ndimage import gaussian_filter  # noqa: PLC0415 - heavy

    counts, xe, ye = np.histogram2d(x, y, bins=DENSITY_BINS, range=[extent[:2], extent[2:]])
    smooth = gaussian_filter(counts.T, DENSITY_SIGMA)
    return (xe[:-1] + xe[1:]) / 2, (ye[:-1] + ye[1:]) / 2, smooth / smooth.sum()


def hdr_levels(density: np.ndarray, masses: Sequence[float] = HDR_MASSES) -> list[float]:
    """The density above which each mass lies: the highest-density region's threshold.

    The region {density ≥ level} is the smallest set of cells holding at least ``mass``.
    """
    ordered = np.sort(density.ravel())[::-1]
    held = np.cumsum(ordered)
    return [float(ordered[min(int(np.searchsorted(held, m)), ordered.size - 1)]) for m in masses]


def above_floor(values: np.ndarray, floor: float | None) -> np.ndarray:
    """Whether each value lies above the floor and so takes a band: every value without one."""
    return np.ones(values.size, dtype=bool) if floor is None else values > floor


def box_labels(labels: Sequence[Any], theme: Theme) -> None:
    """Sit each isoline label on a box of the grid grey, readable on any cividis band."""
    for label in labels:
        label.set_bbox({"facecolor": theme.grid, "edgecolor": "none",
                        "pad": LABEL_BOX_PAD * SIZES["annotation"]})  # fmt: skip


@dataclass(frozen=True)
class Mark:
    """A point on a surface worth naming, e.g. the current setting."""

    x: float
    y: float
    label: str


@dataclass(frozen=True)
class SurfaceData:
    """The contract both renderers read.

    Attributes:
        frame: One row per grid point; every x with every y, once.
        columns: (unit, meaning) for every column of ``frame``.
        x: Column on the x axis.
        y: Column on the y axis.
        z: Column of the value filled in bands.
        x_label: x-axis title, without the unit.
        y_label: y-axis title, without the unit.
        z_label: Colour bar title, without the unit.
        caption: The caption, unwrapped.
        isolines: Values drawn as labelled lines, e.g. 80 % of reads kept.
        isoline_format: Format spec of their labels.
        mark: The current setting, or None.
        width: The column; every figure takes ``FIGURE_WIDTH``, the double column.
        floor: A value at or below which a cell takes no band and is drawn in the reference
            grey, e.g. a setting that keeps no read; None puts every cell in the bands.
        floor_label: The floor's legend key.
        live: Settings of the page control that moves the mark, with its ``kind``; None for a
            surface the page does not change. Set, the PDF and SVG buttons export the view as
            set (D-033).
    """

    frame: pl.DataFrame
    columns: Mapping[str, tuple[str, str]]
    x: str
    y: str
    z: str
    x_label: str
    y_label: str
    z_label: str
    caption: str
    isolines: tuple[float, ...] = ()
    isoline_format: str = ".3g"
    mark: Mark | None = None
    width: str = FIGURE_WIDTH
    floor: float | None = None
    floor_label: str = ""
    live: Mapping[str, Any] | None = None

    def __post_init__(self) -> None:
        """Refuse anything but a full grid of finite values that vary above the floor."""
        if missing := [c for c in (self.x, self.y, self.z) if c not in self.frame.columns]:
            raise ValueError(f"missing columns {missing}")
        nx, ny = self.frame[self.x].n_unique(), self.frame[self.y].n_unique()
        if self.frame.height != nx * ny or self.frame.select(self.x, self.y).is_duplicated().any():
            raise ValueError(f"not a full grid: {self.frame.height} rows for {nx} x {ny}")
        z = self.frame[self.z].cast(pl.Float64)
        if z.is_null().any() or not z.is_finite().all():
            raise ValueError(f"{self.z} must be finite everywhere")
        banded = self._values()[self.above()]
        if not banded.size or banded.min() == banded.max():
            raise ValueError(f"{self.z} must vary above the floor")
        if self.floor is not None and not self.floor_label:
            raise ValueError("a floor needs its legend key")

    def _values(self) -> np.ndarray:
        """The z column, one value per row of ``frame``."""
        return self.frame[self.z].cast(pl.Float64).to_numpy()

    def above(self) -> np.ndarray:
        """Whether each row of ``frame`` lies above the floor and so takes a band."""
        return above_floor(self._values(), self.floor)

    def grid(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """The x values, the y values, and z with rows along y."""
        frame = self.frame.sort(self.y, self.x)
        xs = np.sort(self.frame[self.x].unique().to_numpy())
        ys = np.sort(self.frame[self.y].unique().to_numpy())
        return xs, ys, frame[self.z].to_numpy().astype(float).reshape(len(ys), len(xs))

    def edges(self) -> np.ndarray:
        """The value at each band edge: eight equal bands over the cells above the floor."""
        z = self._values()[self.above()]
        return np.linspace(float(z.min()), float(z.max()), len(CIVIDIS) + 1)

    def table(self) -> pl.DataFrame:
        """Every grid point, with the band it is filled in (0 lowest); none on the floor."""
        edges = self.edges()
        band = np.clip(np.searchsorted(edges, self._values(), "right") - 1, 0, len(CIVIDIS) - 1)
        return self.frame.with_columns(band=pl.when(pl.Series(self.above()))
                                       .then(pl.Series(band, dtype=pl.Int32)))  # fmt: skip

    def table_columns(self) -> Mapping[str, tuple[str, str]]:
        """(unit, meaning) for every column of ``table()``."""
        edges = ", ".join(f"{e:.4g}" for e in self.edges())
        meaning = f"Colour band, 0 lowest, between the edges {edges}"
        if self.floor is not None:
            meaning += f"; empty at or below the floor, {self.floor:g}, drawn grey"
        return {**self.columns, "band": ("none", f"{meaning}.")}

    def page_table(self) -> pl.DataFrame:
        """Every grid point: the page draws the whole surface."""
        return self.table()

    def caption_text(self) -> str:
        """The caption, unwrapped."""
        return self.caption

    def controls(self) -> dict[str, Any] | None:
        """The page control that moves the mark, or None."""
        return dict(self.live) if self.live is not None else None

    def inches(self) -> tuple[float, float]:
        """The figure files' size, which a page's own PDF and SVG of a changed view take."""
        return figure_inches(self.width, "contour")

    def export_layers(self) -> dict[str, Any]:
        """Nothing: the page draws every band, isoline and mark the files draw."""
        return {}

    def _titles(self) -> tuple[str, str, str]:
        """Axis and colour bar titles with their units."""
        return tuple(
            axis_title(label, self.columns[c][0])
            for label, c in ((self.x_label, self.x), (self.y_label, self.y), (self.z_label, self.z))  # type: ignore[return-value]
        )

    def _levels(self) -> list[float]:
        """The band edges matplotlib fills between.

        matplotlib closes a band at its top edge, and at its bottom only where that edge is the
        grid's minimum, which a floor lies below; the lowest edge moves one float down, so the
        lowest cell above the floor is filled in band 0, as ``table()`` puts it, not as floor.
        """
        edges = self.edges()
        return [float(np.nextafter(edges[0], -np.inf)), *(float(e) for e in edges[1:])]

    def static(self, theme: Theme) -> Any:
        """Render the matplotlib figure for one theme: filled bands, isolines, colour bar.

        Args:
            theme: One of ``style.THEMES``.

        Returns:
            A matplotlib Figure, not attached to pyplot.
        """
        from matplotlib.figure import Figure  # noqa: PLC0415 - heavy
        from matplotlib.lines import Line2D  # noqa: PLC0415 - heavy
        from matplotlib.patches import Patch  # noqa: PLC0415 - heavy
        from matplotlib.ticker import FuncFormatter  # noqa: PLC0415 - heavy

        apply_matplotlib(theme)
        width, height = figure_inches(self.width, "contour")
        fig = Figure(figsize=(width, height))
        ax = fig.subplots()
        xs, ys, z = self.grid()
        levels = self._levels()
        filled = ax.contourf(xs, ys, z, levels=levels, colors=list(CIVIDIS))
        if not self.above().all():
            ax.contourf(xs, ys, z, levels=[float(z.min()), levels[0]], colors=[theme.reference])
        if self.isolines:
            lines = ax.contour(xs, ys, z, levels=sorted(self.isolines), colors=theme.ink,
                               linewidths=LINE_WIDTHS["line"])  # fmt: skip
            box_labels(ax.clabel(lines, fmt=lambda v: f"{v:{self.isoline_format}}",
                                 fontsize=SIZES["annotation"]), theme)  # fmt: skip
        handles: list[Any] = []
        if self.mark:
            ax.plot([self.mark.x], [self.mark.y], marker="o", color=theme.ink, linestyle="",
                    markersize=MARKER_SIZE)  # fmt: skip
            handles.append(Line2D([], [], marker="o", linestyle="", color=theme.ink,
                                  label=self.mark.label))  # fmt: skip
        if not self.above().all():
            handles.append(Patch(facecolor=theme.reference, linewidth=LINE_WIDTHS["band_edge"],
                                 label=self.floor_label))  # fmt: skip
        x_title, y_title, z_title = self._titles()
        ax.set_xlabel(x_title)
        ax.set_ylabel(fit_title(y_title, height))
        bar = fig.colorbar(filled, ax=ax)
        bar.set_label(fit_title(z_title, height))
        bar.outline.set_linewidth(LINE_WIDTHS["axes"])
        bar.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.3g}"))
        place_caption(fig, self.caption, width, height, ax=ax, handles=handles or None)
        refit_y_title(fig, ax, y_title)
        refit_y_title(fig, bar.ax, z_title, label=bar.set_label)
        return fig

    def interactive(self) -> Any:
        """Render the plotly figure: the same bands, floor and isolines, hover on every point.

        Isolines are black in both themes: they sit on the cividis fill, which the theme does
        not change. The floor keeps the reference grey, one grey for both themes.

        Returns:
            A plotly Figure with the light template; the page switches themes.
        """
        import plotly.graph_objects as go  # noqa: PLC0415 - heavy; imported after the profile

        xs, ys, z = self.grid()
        edges = self.edges()
        x_title, y_title, z_title = self._titles()
        # plotly colours the bands by floor((zmax - zmin) / size) + 1 levels. The difference of
        # two rounded edges can exceed the band width by an ulp, which dropped a level and
        # shifted bands 3 to 6 up a colour on the fastp grid; a width divided by eight is exact.
        size = (edges[-1] - edges[0]) / (len(edges) - 1)
        fig = go.Figure(go.Contour(
            x=xs, y=ys, z=z, zmin=edges[0], zmax=edges[-1], autocontour=False,
            contours={"start": edges[0], "end": edges[-1], "size": size, "coloring": "fill"},
            colorscale=plotly_steps(list(CIVIDIS)), line={"width": LINE_WIDTHS["band_edge"]},
            colorbar={"title": {"text": z_title}},
            hovertemplate=f"{x_title}: %{{x}}<br>{y_title}: %{{y}}<br>{z_title}: %{{z:.4g}}"
            "<extra></extra>",
        ))  # fmt: skip
        reference = THEMES["white"].reference
        if not self.above().all():
            # A constraint contour shades the side its operation rules out: ">" keeps the
            # values above the lowest edge and shades the floor's side in grey.
            fig.add_trace(go.Contour(
                x=xs, y=ys, z=z, showscale=False, hoverinfo="skip", showlegend=False,
                contours={"type": "constraint", "operation": ">", "value": edges[0]},
                fillcolor=reference, line={"width": LINE_WIDTHS["band_edge"]},
            ))  # fmt: skip
        for value in self.isolines:
            fig.add_trace(go.Contour(
                x=xs, y=ys, z=z, autocontour=False, showscale=False, hoverinfo="skip",
                showlegend=False,
                contours={"start": value, "end": value, "size": 1, "coloring": "none",
                          "showlabels": True, "labelformat": self.isoline_format,
                          "labelfont": {"color": INK_COLOUR}},
                line={"color": INK_COLOUR, "width": LINE_WIDTHS["line"]},
            ))  # fmt: skip
        if self.mark:
            fig.add_trace(go.Scatter(
                x=[self.mark.x], y=[self.mark.y], mode="markers", name=self.mark.label,
                marker={"color": INK_COLOUR, "size": PLOTLY_MARKER_PX},
            ))  # fmt: skip
        if not self.above().all():
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode="markers", name=self.floor_label,
                marker={"symbol": "square", "color": reference, "size": PLOTLY_LEGEND_SWATCH_PX},
            ))  # fmt: skip
        fig.update_layout(
            template=plotly_templates()["white"], height=PLOTLY_HEIGHT_PX["contour"],
            xaxis={"title": {"text": x_title}}, yaxis={"title": {"text": y_title}},
        )  # fmt: skip
        return fig
