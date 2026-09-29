"""Scatter: points coloured by one category, shaped by another, with labels and lines.

One contract draws PCA, volcano, MA, dispersion, LFC-against-LFC and the GSEA running score.
``static()`` and ``interactive()`` read the same ``ScatterData``.
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import polars as pl

from ctx.viz.annotate import fit_title, place_caption
from ctx.viz.axes import axis_title
from ctx.viz.style import (
    BUBBLE_KEYS,
    BUBBLE_PT,
    LINE_WIDTHS,
    MARKER_SIZE,
    MARKERS,
    MAX_CATEGORIES,
    PLOTLY_HEIGHT_PX,
    PLOTLY_MARKER_PX,
    POINT_ALPHA,
    PX_PER_PT,
    SIZES,
    THEMES,
    THRESHOLD_DASH,
    Theme,
    apply_matplotlib,
    figure_inches,
    plotly_templates,
)

# Points not in ``keep`` above this many are thinned on the HTML page; the CSV keeps all.
THIN_MAX = 2000
MAX_LABELS = 10
# Offset of a point label from its point, in points.
LABEL_OFFSET = 3.0
# Hover number formats: significant digits for measurements, thousands for counts.
FLOAT_HOVER, INT_HOVER = ":.3~g", ":,d"
# Cutoffs the page offers as presets; the loosest bounds what may be thinned.
CUTOFF_PRESETS = (0.1, 0.05, 0.01, 0.001)
# The strictest cutoff the page accepts; a smaller typed value is raised to it.
MIN_CUTOFF = 0.0001
CUTOFF_NAME = "cutoff"
# A log axis spanning fewer decades than this labels 1, 2 and 5 of each, not the decades alone.
LOG_SUBS, LOG_SUBS_BELOW = (1.0, 2.0, 5.0), 2.0


@dataclass(frozen=True)
class Line:
    """A line drawn over the points, e.g. a fitted trend, the identity line or a cutoff."""

    name: str
    x: tuple[float, ...]
    y: tuple[float, ...]
    dash: bool = False


def _line_style(line: Line, theme: Theme) -> dict[str, Any]:
    """A cutoff is dashed in the reference grey; a fitted or identity line is ink."""
    if line.dash:
        return {"color": theme.reference, "linestyle": THRESHOLD_DASH[0]}
    return {"color": theme.ink}


@dataclass(frozen=True)
class Calls:
    """What the page needs to call points again, live, under the reader's own cutoffs.

    A point is called when it was tested (it has a value in the first adjusted column) and
    its chosen p-value is below the cutoff; up when its ``lfc`` is at least the fold-change
    threshold, down when at most minus it. ``ctx.de.stats.classify`` is that rule's Python
    definition and colours the default view; the page script repeats it on the embedded rows.

    Attributes:
        lfc: Column of raw log2 fold changes.
        p: Column of raw p-values.
        adjusted: Adjusted-p column to its words (e.g. "BH-adjusted p"); the first is the
            default and marks the tested points.
        raw_words: Words for the raw p-value, uncorrected.
        tested: Name of the tested group, which "up" refers to.
        reference: Name of the reference group.
        alpha: Default cutoff.
        lfc_min: Default fold-change threshold, in log2.
        max_cutoff: The loosest cutoff the controls allow; points no cutoff can call may be
            thinned.
        n_tested: Points the adjustment counted over, drawn or thinned.
        p_floor: p-values below this are drawn at it; -log10 of 0 has no place on an axis.
        p_line: The raw p the default cutoff sits at (``ctx.de.stats.p_line``), or None
            when no point was tested.
        reading: "How to read" lines with ``{placeholders}``, filled for any setting.
        summary: The one-line result with ``{placeholders}``.
        line_rule: How the horizontal line reads, with ``{placeholders}``.
        no_line_rule: The same when no point was tested, so no line is drawn.
        levels: The colour column's level for up, down and not called.
    """

    lfc: str
    p: str
    adjusted: Mapping[str, str]
    raw_words: str
    tested: str
    reference: str
    alpha: float
    lfc_min: float
    max_cutoff: float
    n_tested: int
    p_floor: float
    p_line: float | None
    reading: tuple[str, ...]
    summary: str
    line_rule: str
    no_line_rule: str
    levels: tuple[str, str, str] = ("up", "down", "ns")

    def page_columns(self, hover: tuple[str, ...]) -> list[str]:
        """The columns the page embeds: what it needs to call and to name each point."""
        return list(dict.fromkeys([*hover, self.lfc, self.p, *self.adjusted]))

    def y_line(self) -> float | None:
        """Height of the default horizontal cutoff on the -log10 p axis."""
        return None if self.p_line is None else -math.log10(max(self.p_line, self.p_floor))

    def values(self, n_up: int, n_down: int) -> dict[str, str]:
        """The placeholders at the default setting, formatted as the page script formats them."""
        t, line = self.lfc_min, self.y_line()
        values = {
            "n_tested": f"{self.n_tested:,}", "tested": self.tested, "reference": self.reference,
            "fc": f"{t:.2f}", "fold_up": f"{2**t:.2f}", "fold_down": f"{2**-t:.2f}",
            "pct_up": f"{100 * (2**t - 1):.1f}", "pct_down": f"{100 * (1 - 2**-t):.1f}",
            "basis": next(iter(self.adjusted.values())), "cutoff": _js_number(self.alpha),
            "y_line": "" if line is None else f"{line:.2f}", "n_up": f"{n_up:,}",
            "n_down": f"{n_down:,}",
        }  # fmt: skip
        rule = self.no_line_rule if line is None else self.line_rule
        return values | {"line_rule": rule.format_map(values)}


def bubble_pt(values: list[float], top: float) -> list[float]:
    """Point diameters in pt: ``BUBBLE_PT["min"]`` at 0, ``["max"]`` at ``top``, area linear."""
    low, high = BUBBLE_PT["min"], BUBBLE_PT["max"]
    return [math.sqrt(low**2 + (high**2 - low**2) * max(v, 0.0) / top) if top > 0 else low
            for v in values]  # fmt: skip


def _js_number(value: float) -> str:
    """A number as JavaScript's String() writes it, for the cutoffs the controls offer."""
    return f"{value:.10g}"


def plain_log_ticks(axis: Any, subs: tuple[float, ...] | None = None) -> None:
    """Label a log axis 100, 1,000, 10,000.

    Matplotlib's own log labels are mathtext: the exponent is drawn at 0.7 of the tick size,
    6.3 pt at 9 pt, and the SVG names one font for it instead of the chain.

    Args:
        axis: A matplotlib axis already on a log scale.
        subs: Labelled multiples within each decade, e.g. (1, 2, 5) for an axis spanning less
            than one; None keeps matplotlib's decades.
    """
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter  # noqa: PLC0415

    if subs:
        axis.set_major_locator(LogLocator(base=10, subs=subs))
    axis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.10g}"))
    axis.set_minor_formatter(NullFormatter())


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
        keep: Boolean column of points always drawn on the page; the others are thinned.
        names: Legend text per colour level, where the level itself is terse.
        calls: How the page calls points again under the reader's cutoffs, or None.
        size: Non-negative numeric column that sets each point's area, or None.
        size_label: What the size shows, for its legend keys.
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
    keep: str | None = None
    names: Mapping[str, str] | None = None
    calls: Calls | None = None
    size: str | None = None
    size_label: str = ""

    def __post_init__(self) -> None:
        """Refuse a contract the renderers could only draw wrongly."""
        named = [
            self.x,
            self.y,
            self.colour,
            self.symbol,
            self.label,
            self.keep,
            self.size,
            *self.hover,
        ]
        if self.calls:
            named += [self.calls.lfc, self.calls.p, *self.calls.adjusted]
        if missing := [c for c in named if c is not None and c not in self.frame.columns]:
            raise ValueError(f"missing columns {missing}")
        if len(self.colours) > MAX_CATEGORIES:
            raise ValueError(f"{len(self.colours)} colour levels; at most {MAX_CATEGORIES}")
        if stray := set(self.frame[self.colour].unique()) - set(self.colours):
            raise ValueError(f"{self.colour}: levels {sorted(stray)} have no colour")
        if self.calls and set(self.calls.levels) != set(self.colours):
            raise ValueError(f"call levels {self.calls.levels} are not the colour levels")
        if self.calls and self.symbol:
            # The page script finds one trace per colour level; shapes would split them.
            raise ValueError("live calls cannot be combined with a symbol column")
        if self.size and not (self.frame[self.size].fill_null(0) >= 0).all():
            raise ValueError(f"{self.size}: sizes must be non-negative")
        if self.symbol and self.frame[self.symbol].n_unique() > len(MARKERS):
            raise ValueError(f"{self.symbol}: more levels than the {len(MARKERS)} marker shapes")

    def table(self) -> pl.DataFrame:
        """Every point, drawn or thinned."""
        return self.frame

    def table_columns(self) -> Mapping[str, tuple[str, str]]:
        """(unit, meaning) for every column of ``table()``."""
        return self.columns

    def page_table(self) -> pl.DataFrame:
        """The points the page draws, thinned like the plot so the page stays under its cap.

        With ``calls`` the page embeds only the columns it calls and names points with, so the
        rows stay small; the CSV beside it holds every column.
        """
        rows = self.on_page()
        return rows.select(self.calls.page_columns(self.hover)) if self.calls else rows

    def controls(self) -> dict[str, Any] | None:
        """The live-calling settings the page script reads, or None without ``calls``."""
        if self.calls is None:
            return None
        c = self.calls
        first_line = len(self.colours) + len(self.lines)
        return {
            "lfc": c.lfc, "p": c.p, "gate": next(iter(c.adjusted)),
            "bases": {**c.adjusted, c.p: c.raw_words}, "tested": c.tested,
            "reference": c.reference, "alpha": c.alpha, "lfc_min": c.lfc_min,
            "max_cutoff": c.max_cutoff, "min_cutoff": MIN_CUTOFF, "n_tested": c.n_tested,
            "p_floor": c.p_floor,
            "reading": list(c.reading), "summary": c.summary, "line_rule": c.line_rule,
            "no_line_rule": c.no_line_rule,
            "levels": dict(zip(("up", "down", "ns"), c.levels, strict=True)),
            "traces": {level: i for i, level in enumerate(self.colours)},
            "hover": list(self.hover), "cutoff_lines": [first_line + i for i in range(3)],
            "presets": [_js_number(c) for c in CUTOFF_PRESETS if c <= self.calls.max_cutoff],
            "default_summary": self.summary(),
        }  # fmt: skip

    def reading(self) -> tuple[str, ...]:
        """The "How to read" lines at the default setting; empty without ``calls``."""
        if self.calls is None:
            return ()
        up, down, _ = self.calls.levels
        counts = self.frame[self.colour].value_counts()
        n = dict(counts.rows())
        values = self.calls.values(n.get(up, 0), n.get(down, 0))
        return tuple(line.format_map(values) for line in self.calls.reading)

    def summary(self) -> str:
        """The one-line result at the default setting; empty without ``calls``."""
        if self.calls is None:
            return ""
        up, down, _ = self.calls.levels
        n = dict(self.frame[self.colour].value_counts().rows())
        return self.calls.summary.format_map(self.calls.values(n.get(up, 0), n.get(down, 0)))

    def all_lines(self) -> tuple[Line, ...]:
        """``lines``, then with ``calls`` the three default cutoffs the page script moves."""
        if self.calls is None:
            return self.lines
        t, line = self.calls.lfc_min, self.calls.y_line()
        top = float(self.frame[self.y].to_numpy().max()) * 1.05
        x = self.frame[self.x].to_numpy()
        low, high = min(float(x.min()), -t), max(float(x.max()), t)
        height = (line, line) if line is not None else (math.nan, math.nan)
        horizontal = ((low, high), height)
        return (
            *self.lines,
            # One name, so one legend entry toggles all three.
            Line(CUTOFF_NAME, (t, t), (0.0, top), dash=True),
            Line(CUTOFF_NAME, (-t, -t), (0.0, top), dash=True),
            Line(CUTOFF_NAME, horizontal[0], horizontal[1], dash=True),
        )

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
            f"{other:,} points no setting can call, every {stride}th in order of "
            f"{self.y_label}; the CSV holds all {self.frame.height:,}."
        )

    def _levels(self, column: str | None) -> list[Any]:
        """Levels of a column in first-appearance order, or [None] for no column."""
        return self.frame[column].unique(maintain_order=True).to_list() if column else [None]

    def _axis_title(self, label: str, column: str) -> str:
        """Axis title with the unit from the column dictionary."""
        return axis_title(label, self.columns[column][0])

    def _name(self, level: Any) -> str:
        """Legend text of a colour level."""
        return (self.names or {}).get(level, str(level))

    def _diameters(self, rows: pl.DataFrame) -> list[float] | float:
        """Each point's diameter in pt: from ``size``, else the one marker size."""
        if self.size is None:
            return MARKER_SIZE
        top = float(self.frame[self.size].max() or 0.0)  # type: ignore[arg-type]
        return bubble_pt(rows[self.size].fill_null(0.0).to_list(), top)

    def _size_keys(self) -> list[tuple[str, float]]:
        """(legend text, diameter in pt) of the size keys; none without ``size``."""
        if self.size is None:
            return []
        top = float(self.frame[self.size].max() or 0.0)  # type: ignore[arg-type]
        values = [float(f"{share * top:.2g}") for share in BUBBLE_KEYS]
        return [(f"{self.size_label} {v:g}", d) for v, d in zip(values, bubble_pt(values, top),
                                                                   strict=True)]  # fmt: skip

    def _px(self, rows: pl.DataFrame) -> list[float] | float:
        """Plotly marker diameters in px: the static diameters, on screen."""
        diameters = self._diameters(rows)
        if isinstance(diameters, list):
            return [d * PX_PER_PT for d in diameters]
        return PLOTLY_MARKER_PX

    def _hover_format(self, column: str) -> str:
        """Plotly number format for a column: none for text, digits for numbers."""
        dtype = self.frame.schema[column]
        return FLOAT_HOVER if dtype.is_float() else INT_HOVER if dtype.is_integer() else ""

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
                diameters = self._diameters(rows)
                ax.scatter(
                    rows[self.x].to_numpy(), rows[self.y].to_numpy(),
                    s=[d**2 for d in diameters] if isinstance(diameters, list) else diameters**2,
                    c=colour, marker=MARKERS[j][0], linewidths=LINE_WIDTHS["marker_edge"],
                    alpha=POINT_ALPHA,
                )  # fmt: skip
        for line in self.all_lines():
            if math.isnan(line.y[0]):
                continue
            ax.plot(line.x, line.y, linewidth=LINE_WIDTHS["line"], **_line_style(line, theme))
        ax.set_xscale("log" if self.log_x else "linear")
        ax.set_yscale("log" if self.log_y else "linear")
        for axis, log, column in ((ax.xaxis, self.log_x, self.x), (ax.yaxis, self.log_y, self.y)):
            if log:
                values = self.frame[column].filter(self.frame[column] > 0)
                decades = math.log10(float(values.max()) / float(values.min()))  # type: ignore[arg-type]
                plain_log_ticks(axis, subs=LOG_SUBS if decades < LOG_SUBS_BELOW else None)
        ax.set_xlabel(self._axis_title(self.x_label, self.x))
        ax.set_ylabel(fit_title(self._axis_title(self.y_label, self.y), height))
        present = set(self.frame[self.colour].unique())
        handles = [
            Line2D([], [], marker="o", linestyle="", color=c, label=self._name(level))
            for level, c in colours.items()
            if level in present
        ]
        if self.symbol:
            handles += [
                Line2D([], [], marker=MARKERS[j][0], linestyle="", color=theme.ink, label=str(s))
                for j, s in enumerate(shapes)
            ]
        handles += [
            Line2D([], [], marker="o", linestyle="", color=theme.reference, markersize=d, label=t)
            for t, d in self._size_keys()
        ]
        handles += [
            Line2D([], [], label=line.name, **_line_style(line, theme))
            for line in {ln.name: ln for ln in self.all_lines() if not math.isnan(ln.y[0])}.values()
        ]
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
        x_title, y_title = (
            self._axis_title(self.x_label, self.x),
            self._axis_title(self.y_label, self.y),
        )
        hover = [
            f"{x_title}: %{{x{self._hover_format(self.x)}}}",
            f"{y_title}: %{{y{self._hover_format(self.y)}}}",
            *(
                f"{c}: %{{customdata[{i}]{self._hover_format(c)}}}"
                for i, c in enumerate(self.hover)
            ),
        ]
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
                    name=self._name(level) if shape is None else f"{self._name(level)}, {shape}",
                    marker={"color": self.colours[level] or light.reference,
                            "symbol": MARKERS[j][1], "size": self._px(rows),
                            "opacity": POINT_ALPHA},
                    customdata=rows.select(self.hover).rows(), hovertemplate=template,
                )
            )  # fmt: skip
        # No text labels here: plotly cannot keep them apart; hover names every point.
        named = set()
        for line in self.all_lines():
            # The reference grey, never ink: a trace keeps its colour when the theme flips,
            # and light ink vanishes on the dark page. A cutoff nothing reaches stays as an
            # empty trace, so the page script finds every cutoff at a fixed index.
            drawn = not math.isnan(line.y[0])
            fig.add_trace(
                go.Scatter(
                    x=list(line.x), y=list(line.y) if drawn else [None, None], mode="lines",
                    name=line.name, visible=drawn,
                    line={"color": light.reference, "width": LINE_WIDTHS["line"],
                          "dash": THRESHOLD_DASH[1] if line.dash else "solid"},
                    hoverinfo="skip", showlegend=line.name not in named, legendgroup=line.name,
                )
            )  # fmt: skip
            named.add(line.name)
        fig.update_layout(
            template=plotly_templates()["white"],
            height=PLOTLY_HEIGHT_PX["scatter"],
            showlegend=len(fig.data) > 1,
            xaxis={"title": {"text": x_title}, "type": "log" if self.log_x else "linear"},
            yaxis={"title": {"text": y_title}, "type": "log" if self.log_y else "linear"},
        )  # fmt: skip
        return fig
