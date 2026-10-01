"""Curves: one line per series over a shared x, over the bands a value is judged against.

The QC atlas draws FastQC's per-position quality and its GC distribution this way: one line
per sample, coloured by a category, the tool's bands shaded and its thresholds dashed.
``static()`` and ``interactive()`` read the same ``CurveData``.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import polars as pl

from ctx.viz.annotate import fit_title, place_caption
from ctx.viz.axes import axis_title
from ctx.viz.scatter import FLOAT_HOVER, Line
from ctx.viz.strip import Span, shade, span_handles, span_traces
from ctx.viz.style import (
    DASHES,
    FIGURE_WIDTH,
    LINE_WIDTHS,
    MAX_CATEGORIES,
    PLOTLY_HEIGHT_PX,
    POINT_ALPHA,
    THEMES,
    THRESHOLD_DASH,
    MplDash,
    Theme,
    apply_matplotlib,
    distinct_styles,
    figure_inches,
    plotly_templates,
)


@dataclass(frozen=True)
class CurveData:
    """The contract both renderers read.

    Attributes:
        frame: Long table, one row per point of a series, in x order within each series.
        columns: (unit, meaning) for every column of ``frame``.
        series: Column naming the line a point belongs to, e.g. the sample.
        x: Column on the x axis.
        y: Column on the y axis.
        x_label: x-axis title, without the unit.
        y_label: y-axis title, without the unit.
        colour: Column whose levels colour the lines; constant within a series.
        colours: Level to hex colour, in legend order.
        caption: The caption, unwrapped.
        hover: Columns shown on hover, besides series, x and y.
        spans: Bands shaded behind the lines.
        thresholds: Dashed reference lines, e.g. a tool's warn and fail values.
        width: The column; every figure takes ``FIGURE_WIDTH``, the double column.
        distinct: One unique colour and dash per ``colour`` level, drawn opaque and thick, so
            each line stays traceable; ``colours`` is then unused (``style.distinct_styles``).
    """

    frame: pl.DataFrame
    columns: Mapping[str, tuple[str, str]]
    series: str
    x: str
    y: str
    x_label: str
    y_label: str
    colour: str
    colours: Mapping[str, str]
    caption: str
    hover: tuple[str, ...] = ()
    spans: tuple[Span, ...] = ()
    thresholds: tuple[Line, ...] = ()
    width: str = FIGURE_WIDTH
    distinct: bool = False

    def __post_init__(self) -> None:
        """Refuse a contract the renderers could only draw wrongly."""
        named = [self.series, self.x, self.y, self.colour, *self.hover]
        if missing := [c for c in named if c not in self.frame.columns]:
            raise ValueError(f"missing columns {missing}")
        if self.distinct:
            self._styles()  # raises when there are more levels than distinct styles
        elif len(self.colours) > MAX_CATEGORIES:
            raise ValueError(f"{len(self.colours)} colour levels; at most {MAX_CATEGORIES}")
        if not self.distinct and (
            stray := set(self.frame[self.colour].unique()) - set(self.colours)
        ):
            raise ValueError(f"{self.colour}: levels {sorted(stray)} have no colour")
        mixed = self.frame.group_by(self.series).agg(pl.col(self.colour).n_unique().alias("n"))
        if bad := sorted(mixed.filter(pl.col("n") > 1)[self.series].to_list()):
            raise ValueError(f"series with more than one {self.colour}: {bad}")

    def table(self) -> pl.DataFrame:
        """Every point of every line."""
        return self.frame

    def table_columns(self) -> Mapping[str, tuple[str, str]]:
        """(unit, meaning) for every column of ``table()``."""
        return self.columns

    def page_table(self) -> pl.DataFrame:
        """Every point: the page draws every line."""
        return self.frame

    def caption_text(self) -> str:
        """The caption, unwrapped."""
        return self.caption

    def _styles(self) -> dict[str, tuple[str, tuple[MplDash, str], float, float]]:
        """(colour, (matplotlib dash, plotly dash), width, opacity) per ``colour`` level."""
        levels = [str(v) for v in self.frame[self.colour].unique(maintain_order=True)]
        if self.distinct:
            return {level: (hexa, dash, LINE_WIDTHS["sample"], 1.0)
                    for level, (hexa, dash) in distinct_styles(levels).items()}  # fmt: skip
        return {level: (self.colours[level], DASHES[0], LINE_WIDTHS["line"], POINT_ALPHA)
                for level in levels}  # fmt: skip

    def _lines(self) -> list[tuple[Any, str, pl.DataFrame]]:
        """(series, colour level, rows) per line, in first-appearance order."""
        return [(key[0], str(rows[self.colour][0]), rows)
                for key, rows in self.frame.group_by(self.series, maintain_order=True)]  # fmt: skip

    def _extra(self, rows: pl.DataFrame) -> list[tuple[Any, ...]]:
        """The hover columns of each point; empty tuples when there are none."""
        return rows.select(self.hover).rows() if self.hover else [()] * rows.height

    def _titles(self) -> tuple[str, str]:
        """The two axis titles with their units, as both renderers print them."""
        return (axis_title(self.x_label, self.columns[self.x][0]),
                axis_title(self.y_label, self.columns[self.y][0]))  # fmt: skip

    def static(self, theme: Theme) -> Any:
        """Render the matplotlib figure for one theme.

        Args:
            theme: One of ``style.THEMES``.

        Returns:
            A matplotlib Figure, not attached to pyplot.
        """
        import seaborn as sns  # noqa: PLC0415 - heavy; imported after the runtime profile
        from matplotlib.figure import Figure  # noqa: PLC0415 - heavy
        from matplotlib.lines import Line2D  # noqa: PLC0415 - heavy

        apply_matplotlib(theme)
        width, height = figure_inches(self.width, "scatter")
        fig = Figure(figsize=(width, height))
        ax = fig.subplots()
        styles = self._styles()
        for _, level, rows in self._lines():
            hexa, dash, line_width, opacity = styles[level]
            ax.plot(rows[self.x].to_numpy(), rows[self.y].to_numpy(), color=hexa,
                    linestyle=dash[0], alpha=opacity, linewidth=line_width)  # fmt: skip
        for line in self.thresholds:
            ax.plot(line.x, line.y, color=theme.reference, linestyle=THRESHOLD_DASH[0],
                    linewidth=LINE_WIDTHS["line"])  # fmt: skip
        shade(ax, self.spans)
        x_title, y_title = self._titles()
        ax.set_xlabel(x_title)
        ax.set_ylabel(fit_title(y_title, height))
        handles: list[Any] = [Line2D([], [], color=hexa, linestyle=dash[0], linewidth=lw,
                                     label=level)
                              for level, (hexa, dash, lw, _) in styles.items()]  # fmt: skip
        names = dict.fromkeys(t.name for t in self.thresholds)
        handles += [Line2D([], [], color=theme.reference, linestyle=THRESHOLD_DASH[0], label=n)
                    for n in names]  # fmt: skip
        sns.despine(fig=fig)
        place_caption(fig, self.caption, width, height, ax=ax,
                      handles=handles + span_handles(self.spans))  # fmt: skip
        return fig

    def interactive(self) -> Any:
        """Render the plotly figure; hover names the series and every ``hover`` column.

        Returns:
            A plotly Figure with the light template; the page switches themes.
        """
        import plotly.graph_objects as go  # noqa: PLC0415 - heavy; imported after the profile

        light = THEMES["white"]
        x_values = self.frame[self.x]
        x_range = (float(x_values.min()), float(x_values.max()))  # type: ignore[arg-type]
        fig = go.Figure(span_traces(self.spans, x_range, self.frame[self.y]))
        x_title, y_title = self._titles()
        template = "<br>".join([
            f"{self.series}: %{{customdata[0]}}", f"{x_title}: %{{x}}",
            f"{y_title}: %{{y{FLOAT_HOVER}}}",
            *(f"{c}: %{{customdata[{i + 1}]}}" for i, c in enumerate(self.hover)),
        ]) + "<extra></extra>"  # fmt: skip
        shown: set[str] = set()
        styles = self._styles()
        for name, level, rows in self._lines():
            hexa, dash, line_width, opacity = styles[level]
            fig.add_trace(go.Scatter(
                x=rows[self.x].to_list(), y=rows[self.y].to_list(), mode="lines", name=level,
                legendgroup=level, showlegend=level not in shown, opacity=opacity,
                line={"color": hexa, "width": line_width, "dash": dash[1]},
                customdata=[(str(name), *r) for r in self._extra(rows)],
                hovertemplate=template,
            ))  # fmt: skip
            shown.add(level)
        for line in self.thresholds:
            fig.add_trace(go.Scatter(
                x=list(line.x), y=list(line.y), mode="lines", name=line.name, hoverinfo="skip",
                line={"color": light.reference, "width": LINE_WIDTHS["line"],
                      "dash": THRESHOLD_DASH[1]},
            ))  # fmt: skip
        fig.update_layout(
            template=plotly_templates()["white"], height=PLOTLY_HEIGHT_PX["scatter"],
            xaxis={"title": {"text": x_title}}, yaxis={"title": {"text": y_title}},
        )  # fmt: skip
        return fig
