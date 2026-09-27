"""Scatter: points coloured by one category, shaped by another, with labels and lines.

One contract draws PCA, volcano, MA, dispersion, LFC-against-LFC and the GSEA running score.
``static()`` and ``interactive()`` read the same ``ScatterData``.
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import polars as pl

from ctx.viz.annotate import fit_title, place_caption, plotly_caption
from ctx.viz.style import (
    LINE_WIDTHS,
    MARKER_SIZE,
    MARKERS,
    MAX_CATEGORIES,
    PLOTLY_MARKER_PX,
    SIZES,
    THEMES,
    Theme,
    apply_matplotlib,
    figure_inches,
    plotly_templates,
)

# Points not in ``keep`` above this many are thinned on the HTML page; the CSV keeps all.
THIN_MAX = 2000
MAX_LABELS = 10
# Room above the plot for the threshold slider, in pixels.
SLIDER_MARGIN_PX = 130
# Offset of a point label from its point, in points.
LABEL_OFFSET = 3.0


@dataclass(frozen=True)
class Line:
    """A line drawn over the points, e.g. a fitted trend or the identity line."""

    name: str
    x: tuple[float, ...]
    y: tuple[float, ...]


@dataclass(frozen=True)
class ScatterData:
    """The contract both renderers read.

    Attributes:
        frame: One row per point.
        columns: (unit, meaning) for every column of ``frame``; the axes print the units.
        x: Column on the x axis.
        y: Column on the y axis.
        x_label: x-axis title, without the unit.
        y_label: y-axis title, without the unit.
        colour: Column whose levels colour the points.
        colours: Level to hex colour, in legend order; None is the theme's reference grey.
        caption: The statistics sentence, unwrapped.
        symbol: Column whose levels set the point shape, or None.
        label: Column with text printed beside a point (empty for none), or None.
        hover: Columns shown on hover, besides x and y.
        lines: Lines drawn in ink over the points.
        log_x: Log-scaled x axis.
        log_y: Log-scaled y axis.
        width: ``"single"`` or ``"double"`` column.
        steps: (slider label, column of levels) per step of the page's slider, which recolours
            the points; the step whose column is ``colour`` is the one drawn statically.
        keep: Boolean column of points always drawn on the page; the others are thinned.
    """

    frame: pl.DataFrame
    columns: Mapping[str, tuple[str, str]]
    x: str
    y: str
    x_label: str
    y_label: str
    colour: str
    colours: Mapping[str, str | None]
    caption: str
    symbol: str | None = None
    label: str | None = None
    hover: tuple[str, ...] = ()
    lines: tuple[Line, ...] = ()
    log_x: bool = False
    log_y: bool = False
    width: str = "single"
    steps: tuple[tuple[str, str], ...] = ()
    keep: str | None = None

    def __post_init__(self) -> None:
        """Refuse a contract the renderers could only draw wrongly."""
        named = [self.x, self.y, self.colour, self.symbol, self.label, self.keep, *self.hover]
        named += [column for _, column in self.steps]
        if missing := [c for c in named if c is not None and c not in self.frame.columns]:
            raise ValueError(f"missing columns {missing}")
        if len(self.colours) > MAX_CATEGORIES:
            raise ValueError(f"{len(self.colours)} colour levels; at most {MAX_CATEGORIES}")
        for column in {self.colour, *(c for _, c in self.steps)}:
            if stray := set(self.frame[column].unique()) - set(self.colours):
                raise ValueError(f"{column}: levels {sorted(stray)} have no colour")
        if self.symbol and self.frame[self.symbol].n_unique() > len(MARKERS):
            raise ValueError(f"{self.symbol}: more levels than the {len(MARKERS)} marker shapes")

    def table(self) -> pl.DataFrame:
        """Every point, drawn or thinned."""
        return self.frame

    def table_columns(self) -> Mapping[str, tuple[str, str]]:
        """(unit, meaning) for every column of ``table()``."""
        return self.columns

    def page_table(self) -> pl.DataFrame:
        """The points the page draws: thinned like the plot, so the page stays under its cap."""
        return self.on_page()

    def on_page(self) -> pl.DataFrame:
        """The points the HTML page draws: every kept point, the rest thinned by a stride."""
        if self.keep is None:
            return self.frame
        other = self.frame.filter(~pl.col(self.keep)).sort(self.y)
        if other.height <= THIN_MAX:
            return self.frame
        stride = math.ceil(other.height / THIN_MAX)
        return pl.concat([self.frame.filter(pl.col(self.keep)), other.gather_every(stride)])

    def caption_text(self) -> str:
        """The caption, plus the thinning rule when the page draws fewer points than the CSV."""
        drawn = self.on_page().height
        if drawn == self.frame.height or self.keep is None:
            return self.caption
        other = self.frame.filter(~pl.col(self.keep)).height
        stride = math.ceil(other / THIN_MAX)
        return (
            f"{self.caption} This page draws {drawn - (self.frame.height - other):,} of the "
            f"{other:,} points outside the called set, every {stride}th in order of "
            f"{self.y_label}; the CSV holds all {self.frame.height:,}."
        )

    def _levels(self, column: str | None) -> list[Any]:
        """Levels of a column in first-appearance order, or [None] for no column."""
        return self.frame[column].unique(maintain_order=True).to_list() if column else [None]

    def _axis_title(self, label: str, column: str) -> str:
        """Axis title with the unit from the column dictionary."""
        return f"{label} ({self.columns[column][0]})"

    def static(self, theme: Theme) -> Any:
        """Render the matplotlib figure for one theme.

        Args:
            theme: One of ``style.THEMES``.

        Returns:
            A matplotlib Figure, not attached to pyplot.
        """
        import seaborn as sns  # noqa: PLC0415 - heavy; imported after the runtime profile
        from matplotlib.backends.backend_agg import FigureCanvasAgg  # noqa: PLC0415 - heavy
        from matplotlib.figure import Figure  # noqa: PLC0415 - heavy
        from matplotlib.lines import Line2D  # noqa: PLC0415 - heavy

        apply_matplotlib(theme)
        width, height = figure_inches(self.width, "scatter")
        fig = Figure(figsize=(width, height))
        ax = fig.subplots()
        colours = {level: colour or theme.reference for level, colour in self.colours.items()}
        shapes = self._levels(self.symbol)
        for level, colour in colours.items():
            for j, shape in enumerate(shapes):
                rows = self.frame.filter(pl.col(self.colour) == level)
                if self.symbol:
                    rows = rows.filter(pl.col(self.symbol) == shape)
                ax.scatter(
                    rows[self.x].to_numpy(), rows[self.y].to_numpy(), s=MARKER_SIZE**2,
                    c=colour, marker=MARKERS[j][0], linewidths=LINE_WIDTHS["marker_edge"],
                )  # fmt: skip
        for line in self.lines:
            ax.plot(line.x, line.y, color=theme.ink, linewidth=LINE_WIDTHS["line"])
        ax.set_xscale("log" if self.log_x else "linear")
        ax.set_yscale("log" if self.log_y else "linear")
        ax.set_xlabel(self._axis_title(self.x_label, self.x))
        ax.set_ylabel(fit_title(self._axis_title(self.y_label, self.y), height))
        present = set(self.frame[self.colour].unique())
        handles = [
            Line2D([], [], marker="o", linestyle="", color=c, label=str(level))
            for level, c in colours.items()
            if level in present
        ]
        if self.symbol:
            handles += [
                Line2D([], [], marker=MARKERS[j][0], linestyle="", color=theme.ink, label=str(s))
                for j, s in enumerate(shapes)
            ]
        handles += [Line2D([], [], color=theme.ink, label=line.name) for line in self.lines]
        sns.despine(fig=fig)
        place_caption(fig, self.caption, width, height, ax=ax, handles=handles)
        if self.label:
            renderer = FigureCanvasAgg(fig).get_renderer()  # type: ignore[no-untyped-call]
            self._place_labels(fig, ax, renderer, theme)
        return fig

    def _place_labels(self, fig: Any, ax: Any, renderer: Any, theme: Theme) -> None:
        """Label up to ``MAX_LABELS`` points inside the axes, clear of every text on the page.

        A label that finds no free place is left out; the point keeps its name on hover and
        in the CSV. A leader line joins each label to its point.
        """
        from matplotlib.text import Text  # noqa: PLC0415 - heavy; imported after the profile

        assert self.label is not None  # noqa: S101 - called only when a label column is set
        rows = self.frame.filter(pl.col(self.label) != "").sort(self.y, descending=True)
        fig.draw(renderer)  # the legend lays its texts out only when drawn
        placed = [
            text.get_window_extent(renderer)
            for text in fig.findobj(Text)
            if text.get_visible() and text.get_text().strip()
        ]
        frame = ax.get_window_extent(renderer)
        centre = ax.get_window_extent(renderer).x0 + ax.get_window_extent(renderer).width / 2
        for text, x, y in rows.select(self.label, self.x, self.y).head(MAX_LABELS).rows():
            # Labels face the middle of the axes, so a point near an edge keeps its label on
            # the page.
            right = ax.transData.transform((x, y))[0] > centre
            dx = -LABEL_OFFSET if right else LABEL_OFFSET
            note = ax.annotate(
                text, (x, y), xytext=(dx, LABEL_OFFSET), textcoords="offset points",
                ha="right" if right else "left", fontsize=SIZES["annotation"], color=theme.ink,
                arrowprops={"arrowstyle": "-", "color": theme.ink,
                            "linewidth": LINE_WIDTHS["line"], "shrinkA": 0, "shrinkB": 0},
            )  # fmt: skip
            for attempt in range(1, 4 * MAX_LABELS):
                box = note.get_window_extent(renderer)
                inside = frame.x0 <= box.x0 and box.x1 <= frame.x1 and frame.y0 <= box.y0
                if inside and box.y1 <= frame.y1 and not any(box.overlaps(o) for o in placed):
                    placed.append(box)
                    break
                # Alternate above and below, one label height further each time.
                step = box.height * 72 / fig.dpi * ((attempt + 1) // 2)
                note.xyann = (dx, LABEL_OFFSET + (step if attempt % 2 else -step))
            else:
                note.remove()

    def _trace_rows(self, frame: pl.DataFrame, column: str) -> list[pl.DataFrame]:
        """The rows of each (colour level, shape) trace when ``column`` sets the colour."""
        out = []
        for level in self.colours:
            for shape in self._levels(self.symbol):
                rows = frame.filter(pl.col(column) == level)
                out.append(rows.filter(pl.col(self.symbol) == shape) if self.symbol else rows)
        return out

    def interactive(self) -> Any:
        """Render the plotly figure; hover shows x, y and every ``hover`` column.

        Returns:
            A plotly Figure with the light template; the page switches themes.
        """
        import plotly.graph_objects as go  # noqa: PLC0415 - heavy; imported after the profile

        light = THEMES["white"]
        frame = self.on_page()
        hover = [f"{self.x_label}: %{{x}}", f"{self.y_label}: %{{y}}"]
        hover += [f"{c}: %{{customdata[{i}]}}" for i, c in enumerate(self.hover)]
        template = "<br>".join(hover) + "<extra></extra>"
        fig = go.Figure()
        shapes = self._levels(self.symbol)
        pairs = [(level, j, s) for level in self.colours for j, s in enumerate(shapes)]
        for (level, j, shape), rows in zip(
            pairs, self._trace_rows(frame, self.colour), strict=True
        ):
            fig.add_trace(
                go.Scatter(
                    x=rows[self.x].to_list(), y=rows[self.y].to_list(), mode="markers",
                    name=str(level) if shape is None else f"{level}, {shape}",
                    marker={"color": self.colours[level] or light.reference,
                            "symbol": MARKERS[j][1], "size": PLOTLY_MARKER_PX},
                    customdata=rows.select(self.hover).rows(), hovertemplate=template,
                )
            )  # fmt: skip
        # No text labels here: plotly cannot keep them apart; hover names every point.
        for line in self.lines:
            fig.add_trace(
                go.Scatter(
                    x=list(line.x), y=list(line.y), mode="lines", name=line.name,
                    line={"color": light.reference, "width": LINE_WIDTHS["line"]}, hoverinfo="skip",
                )
            )  # fmt: skip
        caption, bottom = plotly_caption(self.caption_text(), figure_inches("double", "scatter")[0])
        fig.update_layout(
            template=plotly_templates()["white"],
            xaxis={"title": self._axis_title(self.x_label, self.x),
                   "type": "log" if self.log_x else "linear"},
            yaxis={"title": self._axis_title(self.y_label, self.y),
                   "type": "log" if self.log_y else "linear"},
            annotations=[caption],
            margin={"b": bottom},
        )  # fmt: skip
        if self.steps:
            fig.update_layout(sliders=[self._slider(frame)], margin={"t": SLIDER_MARGIN_PX})
        return fig

    def _slider(self, frame: pl.DataFrame) -> dict[str, Any]:
        """A slider whose steps restyle each colour trace to that step's calls; no script."""
        indices = list(range(len(self.colours) * len(self._levels(self.symbol))))
        steps = []
        for label, column in self.steps:
            rows = self._trace_rows(frame, column)
            restyle = {
                "x": [r[self.x].to_list() for r in rows],
                "y": [r[self.y].to_list() for r in rows],
                "customdata": [r.select(self.hover).rows() for r in rows],
            }
            steps.append({"label": label, "method": "restyle", "args": [restyle, indices]})
        active = next(i for i, (_, column) in enumerate(self.steps) if column == self.colour)
        # Above the plot: below it, the slider covers the caption the figure carries.
        return {
            "active": active,
            "steps": steps,
            "currentvalue": {"prefix": "threshold: "},
            "x": 0,
            "y": 1.02,
            "yanchor": "bottom",
            "len": 0.6,
            "pad": {"t": 0, "b": 10},
        }
