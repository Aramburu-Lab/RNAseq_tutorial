"""Strip plot: one value per sample, per category, over the bands it is judged against.

The QC atlas draws one per metric: every sample a point, the metric's absolute bands shaded
across the axis where a source gives one, and each stratum's median ± k MAD as a grey bar.
``static()`` and ``interactive()`` read the same ``StripData``.
"""

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import polars as pl

from ctx.viz.annotate import fit_title, place_caption
from ctx.viz.axes import axis_title
from ctx.viz.scatter import FLOAT_HOVER, log_ticks
from ctx.viz.style import (
    BAND_ALPHA,
    BAND_COLOURS,
    BAR_LABEL_ROTATION,
    FIGURE_WIDTH,
    LINE_WIDTHS,
    MARKER_SIZE,
    MARKERS,
    MAX_CATEGORIES,
    PLOTLY_HEIGHT_PX,
    PLOTLY_MARKER_PX,
    POINT_ALPHA,
    STRIP,
    THEMES,
    Theme,
    apply_matplotlib,
    figure_inches,
    plotly_templates,
    tint,
)

# The golden ratio's fraction: successive points of a category land far apart in its slot.
SPREAD = (math.sqrt(5) - 1) / 2
RANGE_COLUMNS = {
    "range_median": "Median of the category's values, the centre of its grey bar.",
    "range_low": "Lower end of the category's grey bar.",
    "range_high": "Upper end of the category's grey bar.",
}


@dataclass(frozen=True)
class Span:
    """An absolute band shaded across the whole axis; None runs to the axis's end."""

    band: str
    label: str
    low: float | None
    high: float | None


def shade(ax: Any, spans: Sequence[Span]) -> None:
    """Shade absolute bands across the axes, behind everything, within the current y range."""
    bottom, top = ax.get_ylim()
    for span in spans:
        low = bottom if span.low is None else max(span.low, bottom)
        high = top if span.high is None else min(span.high, top)
        if low < high:
            ax.axhspan(low, high, color=BAND_COLOURS[span.band], alpha=BAND_ALPHA,
                       linewidth=LINE_WIDTHS["band_edge"], zorder=0)  # fmt: skip
    ax.set_ylim(bottom, top)


def span_handles(spans: Sequence[Span]) -> list[Any]:
    """A legend patch per band, tinted as it is shaded."""
    from matplotlib.patches import Patch  # noqa: PLC0415 - heavy

    return [Patch(facecolor=BAND_COLOURS[s.band], alpha=BAND_ALPHA, label=s.label) for s in spans]


def span_traces(spans: Sequence[Span], x: tuple[float, float], values: pl.Series) -> list[Any]:
    """The bands as filled plotly traces across ``x``, an open end reaching past the values.

    Filled traces, not layout shapes, so a band sits in data coordinates on a log axis too.
    """
    import plotly.graph_objects as go  # noqa: PLC0415 - heavy

    present = values.drop_nulls()
    bottom, top = float(present.min()), float(present.max())  # type: ignore[arg-type]
    traces = []
    for span in spans:
        low = min(bottom, span.high or bottom) if span.low is None else span.low
        high = max(top, span.low or top) if span.high is None else span.high
        traces.append(go.Scatter(
            x=[x[0], x[1], x[1], x[0], x[0]], y=[low, low, high, high, low], fill="toself",
            mode="lines", line={"width": LINE_WIDTHS["band_edge"]}, name=span.label,
            fillcolor=tint(BAND_COLOURS[span.band], BAND_ALPHA), hoverinfo="skip",
        ))  # fmt: skip
    return traces


@dataclass(frozen=True)
class Range:
    """One category's own band: its median and the bar's two ends."""

    level: str
    median: float
    low: float
    high: float


@dataclass(frozen=True)
class StripData:
    """The contract both renderers read.

    Attributes:
        frame: One row per point, in the order categories are drawn.
        columns: (unit, meaning) for every column of ``frame``.
        category: Column whose levels are the slots along x.
        value: Column on the y axis.
        category_label: x-axis title.
        value_label: y-axis title, without the unit.
        colour: Column whose levels colour the points.
        colours: Level to hex colour, in legend order.
        caption: The caption, unwrapped.
        flagged: Boolean column; a flagged point is drawn as a triangle.
        flag_label: Legend text of the triangle.
        hover: Columns shown on hover, besides category and value.
        spans: Absolute bands, shaded behind everything.
        ranges: Per-category bars, e.g. median ± 3 MAD.
        range_label: Legend text of the grey bar.
        log: Log-scaled value axis.
        width: The column; every figure takes ``FIGURE_WIDTH``, the double column.
    """

    frame: pl.DataFrame
    columns: Mapping[str, tuple[str, str]]
    category: str
    value: str
    category_label: str
    value_label: str
    colour: str
    colours: Mapping[str, str]
    caption: str
    flagged: str | None = None
    flag_label: str = ""
    hover: tuple[str, ...] = ()
    spans: tuple[Span, ...] = ()
    ranges: tuple[Range, ...] = ()
    range_label: str = ""
    log: bool = False
    width: str = FIGURE_WIDTH

    def __post_init__(self) -> None:
        """Refuse a contract the renderers could only draw wrongly."""
        named = [self.category, self.value, self.colour, self.flagged, *self.hover]
        if missing := [c for c in named if c is not None and c not in self.frame.columns]:
            raise ValueError(f"missing columns {missing}")
        if len(self.colours) > MAX_CATEGORIES:
            raise ValueError(f"{len(self.colours)} colour levels; at most {MAX_CATEGORIES}")
        if stray := set(self.frame[self.colour].unique()) - set(self.colours):
            raise ValueError(f"{self.colour}: levels {sorted(stray)} have no colour")
        if bands := {s.band for s in self.spans} - set(BAND_COLOURS):
            raise ValueError(f"spans name no band colour: {sorted(bands)}")
        if levels := {r.level for r in self.ranges} - set(self.levels()):
            raise ValueError(f"ranges for absent categories: {sorted(levels)}")

    def levels(self) -> list[str]:
        """Category levels, in first-appearance order: the slots from left to right."""
        return [str(v) for v in self.frame[self.category].unique(maintain_order=True)]

    def positions(self) -> list[float]:
        """Each point's x: its slot, spread by the golden ratio so points seldom cover."""
        slot = {level: i for i, level in enumerate(self.levels())}
        seen: dict[str, int] = {}
        out = []
        for level in self.frame[self.category].cast(pl.String):
            k = seen[level] = seen.get(level, -1) + 1
            out.append(slot[level] + ((k * SPREAD) % 1 - 0.5) * 2 * STRIP["jitter"])
        return out

    def table(self) -> pl.DataFrame:
        """Every point, with its category's bar where one is drawn."""
        bars = pl.DataFrame(
            [(r.level, r.median, r.low, r.high) for r in self.ranges],
            schema={self.category: pl.String, **dict.fromkeys(RANGE_COLUMNS, pl.Float64)},
            orient="row",
        )
        frame = self.frame.with_columns(pl.col(self.category).cast(pl.String))
        return frame.join(bars, on=self.category, how="left", maintain_order="left")

    def table_columns(self) -> Mapping[str, tuple[str, str]]:
        """(unit, meaning) for every column of ``table()``."""
        unit = self.columns[self.value][0]
        return {**self.columns, **{c: (unit, m) for c, m in RANGE_COLUMNS.items()}}

    def page_table(self) -> pl.DataFrame:
        """Every point: a strip plot draws them all."""
        return self.table()

    def caption_text(self) -> str:
        """The caption, unwrapped."""
        return self.caption

    def _value_title(self) -> str:
        """The value axis title with its unit, as both renderers print it."""
        return axis_title(self.value_label, self.columns[self.value][0])

    def _flags(self) -> list[bool]:
        """Whether each point is flagged; none are without a ``flagged`` column."""
        return self.frame[self.flagged].to_list() if self.flagged else [False] * self.frame.height

    def static(self, theme: Theme) -> Any:
        """Render the matplotlib figure for one theme.

        Args:
            theme: One of ``style.THEMES``.

        Returns:
            A matplotlib Figure, not attached to pyplot.
        """
        import seaborn as sns  # noqa: PLC0415 - heavy; imported after the runtime profile
        from matplotlib.figure import Figure  # noqa: PLC0415 - heavy

        apply_matplotlib(theme)
        width, height = figure_inches(self.width, "strip")
        fig = Figure(figsize=(width, height))
        ax = fig.subplots()
        x, flags = self.positions(), self._flags()
        levels = self.frame[self.colour].to_list()
        values = self.frame[self.value].to_list()
        for level, colour in self.colours.items():
            for flag, (marker, _) in ((False, MARKERS[0]), (True, MARKERS[1])):
                keep = [i for i, v in enumerate(levels) if v == level and flags[i] == flag]
                ax.scatter([x[i] for i in keep], [values[i] for i in keep], s=MARKER_SIZE**2,
                           c=colour, marker=marker, alpha=POINT_ALPHA,
                           linewidths=LINE_WIDTHS["marker_edge"])  # fmt: skip
        ax.set_yscale("log" if self.log else "linear")
        if self.log:
            log_ticks(ax.yaxis, self.frame[self.value])
        self._draw_bands(ax, theme)
        slots = self.levels()
        ax.set_xticks(range(len(slots)), slots, rotation=BAR_LABEL_ROTATION)
        ax.set_xlim(-0.5, len(slots) - 0.5)
        ax.set_xlabel(self.category_label)
        ax.set_ylabel(fit_title(self._value_title(), height))
        sns.despine(fig=fig)
        place_caption(fig, self.caption, width, height, ax=ax, handles=self._handles(theme))
        return fig

    def _draw_bands(self, ax: Any, theme: Theme) -> None:
        """Absolute bands across the axis, then each category's bar and median, behind points."""
        from matplotlib.patches import Rectangle  # noqa: PLC0415 - heavy

        bottom, top = ax.get_ylim()
        shade(ax, self.spans)
        slot = {level: i for i, level in enumerate(self.levels())}
        half = STRIP["range"]
        for r in self.ranges:
            i = slot[r.level]
            ax.add_patch(
                Rectangle(
                    (i - half, r.low),
                    2 * half,
                    r.high - r.low,
                    zorder=0,
                    facecolor=theme.reference,
                    alpha=BAND_ALPHA,
                    linewidth=LINE_WIDTHS["band_edge"],
                )
            )
            ax.plot([i - half, i + half], [r.median, r.median], color=theme.reference,
                    linewidth=LINE_WIDTHS["line"], zorder=1)  # fmt: skip
        ax.set_ylim(bottom, top)

    def _handles(self, theme: Theme) -> list[Any]:
        """Legend entries: colour levels, the flag triangle, the bar, the shaded bands."""
        from matplotlib.lines import Line2D  # noqa: PLC0415 - heavy
        from matplotlib.patches import Patch  # noqa: PLC0415 - heavy

        present = set(self.frame[self.colour].unique())
        handles: list[Any] = [Line2D([], [], marker=MARKERS[0][0], linestyle="", color=c,
                                     label=str(level))
                              for level, c in self.colours.items() if level in present]  # fmt: skip
        if any(self._flags()):
            handles.append(Line2D([], [], marker=MARKERS[1][0], linestyle="", color=theme.ink,
                                  label=self.flag_label))  # fmt: skip
        if self.ranges:
            handles.append(Patch(facecolor=theme.reference, alpha=BAND_ALPHA,
                                 label=self.range_label))  # fmt: skip
        return handles + span_handles(self.spans)

    def interactive(self) -> Any:
        """Render the plotly figure; hover shows the category, the value and ``hover``.

        Bands are filled traces, not layout shapes, so they sit in data coordinates on a log
        axis as on a linear one.

        Returns:
            A plotly Figure with the light template; the page switches themes.
        """
        import plotly.graph_objects as go  # noqa: PLC0415 - heavy; imported after the profile

        light, slots = THEMES["white"], self.levels()
        fig = go.Figure()
        for trace in self._band_traces(light, len(slots)):
            fig.add_trace(trace)
        x, flags = self.positions(), self._flags()
        rows = self.frame.with_columns(_x=pl.Series(x), _flag=pl.Series(flags, dtype=pl.Boolean))
        title = self._value_title()
        template = "<br>".join([
            f"{self.category_label}: %{{customdata[0]}}", f"{title}: %{{y{FLOAT_HOVER}}}",
            *(f"{c}: %{{customdata[{i + 1}]}}" for i, c in enumerate(self.hover)),
        ]) + "<extra></extra>"  # fmt: skip
        for level, colour in self.colours.items():
            part = rows.filter(pl.col(self.colour) == level)
            fig.add_trace(go.Scatter(
                x=part["_x"].to_list(), y=part[self.value].to_list(), mode="markers",
                name=str(level), hovertemplate=template,
                customdata=part.select(self.category, *self.hover).rows(),
                marker={"color": colour, "size": PLOTLY_MARKER_PX, "opacity": POINT_ALPHA,
                        "symbol": [MARKERS[int(f)][1] for f in part["_flag"]]},
            ))  # fmt: skip
        fig.update_layout(
            template=plotly_templates()["white"], height=PLOTLY_HEIGHT_PX["strip"],
            xaxis={"title": {"text": self.category_label}, "tickvals": list(range(len(slots))),
                   "ticktext": slots, "range": [-0.5, len(slots) - 0.5]},
            yaxis={"title": {"text": title}, "type": "log" if self.log else "linear"},
        )  # fmt: skip
        return fig

    def _band_traces(self, light: Theme, n: int) -> list[Any]:
        """The shaded bands and bars as filled traces, legend entries once each."""
        import plotly.graph_objects as go  # noqa: PLC0415 - heavy

        traces = span_traces(self.spans, (-0.5, n - 0.5), self.frame[self.value])
        slot = {level: i for i, level in enumerate(self.levels())}
        half = STRIP["range"]
        for k, r in enumerate(self.ranges):
            i = slot[r.level]
            traces.append(go.Scatter(
                x=[i - half, i + half, i + half, i - half, i - half],
                y=[r.low, r.low, r.high, r.high, r.low], fill="toself", mode="lines",
                name=self.range_label, legendgroup="range", showlegend=k == 0,
                fillcolor=tint(light.reference, BAND_ALPHA), line={"width": 0},
                hoverinfo="skip",
            ))  # fmt: skip
            traces.append(go.Scatter(
                x=[i - half, i + half], y=[r.median, r.median], mode="lines",
                legendgroup="range", showlegend=False, hoverinfo="skip",
                line={"color": light.reference, "width": LINE_WIDTHS["line"]},
            ))  # fmt: skip
        return traces
