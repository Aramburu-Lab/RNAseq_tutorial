"""Contours of a continuous value over two continuous axes: computed surfaces and densities.

``SurfaceData`` fills a value computed on a grid (a two-knob map, a power map) in eight
cividis bands, with labelled isolines and the current setting marked. The density helpers,
``density_grid`` and ``hdr_levels``, compute once the smoothed histogram and the 50 % and 90 %
highest-density levels that ``ctx.viz.density`` draws in both renderers.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import polars as pl

from ctx.viz.annotate import fit_title, place_caption
from ctx.viz.axes import axis_title
from ctx.viz.style import (
    CIVIDIS,
    INK_COLOUR,
    LABEL_BOX_PAD,
    LINE_WIDTHS,
    MARKER_SIZE,
    PLOTLY_HEIGHT_PX,
    PLOTLY_MARKER_PX,
    SIZES,
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
        width: ``"single"`` or ``"double"`` column.
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
    width: str = "single"

    def __post_init__(self) -> None:
        """Refuse anything but a full grid of finite values that vary."""
        if missing := [c for c in (self.x, self.y, self.z) if c not in self.frame.columns]:
            raise ValueError(f"missing columns {missing}")
        nx, ny = self.frame[self.x].n_unique(), self.frame[self.y].n_unique()
        if self.frame.height != nx * ny or self.frame.select(self.x, self.y).is_duplicated().any():
            raise ValueError(f"not a full grid: {self.frame.height} rows for {nx} x {ny}")
        z = self.frame[self.z].cast(pl.Float64)
        if z.is_null().any() or not z.is_finite().all() or z.min() == z.max():
            raise ValueError(f"{self.z} must be finite everywhere and vary")

    def grid(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """The x values, the y values, and z with rows along y."""
        frame = self.frame.sort(self.y, self.x)
        xs = np.sort(self.frame[self.x].unique().to_numpy())
        ys = np.sort(self.frame[self.y].unique().to_numpy())
        return xs, ys, frame[self.z].to_numpy().astype(float).reshape(len(ys), len(xs))

    def edges(self) -> np.ndarray:
        """The value at each band edge: eight equal bands from the lowest to the highest."""
        z = self.frame[self.z].cast(pl.Float64)
        return np.linspace(float(z.min()), float(z.max()), len(CIVIDIS) + 1)  # type: ignore[arg-type]

    def table(self) -> pl.DataFrame:
        """Every grid point, with the band it is filled in (0 lowest)."""
        edges = self.edges()
        band = np.clip(np.searchsorted(edges, self.frame[self.z].to_numpy(), "right") - 1,
                       0, len(CIVIDIS) - 1)  # fmt: skip
        return self.frame.with_columns(band=pl.Series(band, dtype=pl.Int32))

    def table_columns(self) -> Mapping[str, tuple[str, str]]:
        """(unit, meaning) for every column of ``table()``."""
        edges = ", ".join(f"{e:.4g}" for e in self.edges())
        return {**self.columns,
                "band": ("none", f"Colour band, 0 lowest, between the edges {edges}.")}  # fmt: skip

    def page_table(self) -> pl.DataFrame:
        """Every grid point: the page draws the whole surface."""
        return self.table()

    def caption_text(self) -> str:
        """The caption, unwrapped."""
        return self.caption

    def _titles(self) -> tuple[str, str, str]:
        """Axis and colour bar titles with their units."""
        return tuple(
            axis_title(label, self.columns[c][0])
            for label, c in ((self.x_label, self.x), (self.y_label, self.y), (self.z_label, self.z))  # type: ignore[return-value]
        )

    def static(self, theme: Theme) -> Any:
        """Render the matplotlib figure for one theme: filled bands, isolines, colour bar.

        Args:
            theme: One of ``style.THEMES``.

        Returns:
            A matplotlib Figure, not attached to pyplot.
        """
        from matplotlib.figure import Figure  # noqa: PLC0415 - heavy
        from matplotlib.lines import Line2D  # noqa: PLC0415 - heavy
        from matplotlib.ticker import FuncFormatter  # noqa: PLC0415 - heavy

        apply_matplotlib(theme)
        width, height = figure_inches(self.width, "contour")
        fig = Figure(figsize=(width, height))
        ax = fig.subplots()
        xs, ys, z = self.grid()
        filled = ax.contourf(xs, ys, z, levels=self.edges(), colors=list(CIVIDIS))
        if self.isolines:
            lines = ax.contour(xs, ys, z, levels=sorted(self.isolines), colors=theme.ink,
                               linewidths=LINE_WIDTHS["line"])  # fmt: skip
            box_labels(ax.clabel(lines, fmt=lambda v: f"{v:{self.isoline_format}}",
                                 fontsize=SIZES["annotation"]), theme)  # fmt: skip
        handles = []
        if self.mark:
            ax.plot([self.mark.x], [self.mark.y], marker="o", color=theme.ink, linestyle="",
                    markersize=MARKER_SIZE)  # fmt: skip
            handles.append(Line2D([], [], marker="o", linestyle="", color=theme.ink,
                                  label=self.mark.label))  # fmt: skip
        x_title, y_title, z_title = self._titles()
        ax.set_xlabel(x_title)
        ax.set_ylabel(fit_title(y_title, height))
        bar = fig.colorbar(filled, ax=ax)
        bar.set_label(fit_title(z_title, height))
        bar.outline.set_linewidth(LINE_WIDTHS["axes"])
        bar.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.3g}"))
        place_caption(fig, self.caption, width, height, ax=ax, handles=handles or None)
        return fig

    def interactive(self) -> Any:
        """Render the plotly figure: the same bands and isolines, hover on every grid point.

        Isolines are black in both themes: they sit on the cividis fill, which the theme does
        not change.

        Returns:
            A plotly Figure with the light template; the page switches themes.
        """
        import plotly.graph_objects as go  # noqa: PLC0415 - heavy; imported after the profile

        xs, ys, z = self.grid()
        edges = self.edges()
        x_title, y_title, z_title = self._titles()
        fig = go.Figure(go.Contour(
            x=xs, y=ys, z=z, zmin=edges[0], zmax=edges[-1], autocontour=False,
            contours={"start": edges[0], "end": edges[-1], "size": edges[1] - edges[0],
                      "coloring": "fill"},
            colorscale=plotly_steps(list(CIVIDIS)), line={"width": LINE_WIDTHS["band_edge"]},
            colorbar={"title": {"text": z_title}},
            hovertemplate=f"{x_title}: %{{x}}<br>{y_title}: %{{y}}<br>{z_title}: %{{z:.4g}}"
            "<extra></extra>",
        ))  # fmt: skip
        for value in self.isolines:
            fig.add_trace(go.Contour(
                x=xs, y=ys, z=z, autocontour=False, showscale=False, hoverinfo="skip",
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
        fig.update_layout(
            template=plotly_templates()["white"], height=PLOTLY_HEIGHT_PX["contour"],
            xaxis={"title": {"text": x_title}}, yaxis={"title": {"text": y_title}},
        )  # fmt: skip
        return fig
