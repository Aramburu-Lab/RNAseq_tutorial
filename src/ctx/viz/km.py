"""Kaplan-Meier curves with their bands, censor marks and the number at risk under the plot.

One step line per group from 1 at time 0, its log(-log) Greenwood band shaded behind it, a
vertical tick where a patient is censored, and below the plot one row per group with the
patients still followed at each of ``risk_times``, on the plot's own time axis, so a zoom on
the page moves both. ``static()`` and ``interactive()`` read the same ``KmData``; with
``live`` set, the M10 page's survival panel redraws the traces in place (``_live.js.j2``) and
the PDF and SVG buttons export the view as set (D-033).
"""

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

import polars as pl

from ctx.viz.annotate import fit_title, place_caption
from ctx.viz.axes import axis_title, page_y_title, screen_title
from ctx.viz.style import (
    BAND_ALPHA,
    FIGURE_WIDTH,
    KM,
    LINE_WIDTHS,
    MAX_CATEGORIES,
    PAGE_PX_PER_PT,
    PLOTLY_HEIGHT_PX,
    SIZES,
    Theme,
    apply_matplotlib,
    figure_inches,
    plotly_templates,
)

COLUMNS = ("group", "time", "at_risk", "events", "censored", "survival", "lower", "upper")
# Decimals of the rows a page embeds: months to about an hour, fractions to 1e-5.
PAGE_DIGITS = {"time": 4, "fraction": 5}


@dataclass(frozen=True)
class KmData:
    """The contract both renderers read.

    Attributes:
        frame: One row per group and distinct time: ``COLUMNS``, as
            ``ctx.clinical.survival.Curve`` holds them, in time order within a group.
        columns: (unit, meaning) for every column of ``frame``.
        groups: Group names in legend order; a group without rows draws nothing.
        colours: Group to hex colour.
        x_label: Time axis title, without the unit.
        y_label: Survival axis title, without the unit.
        risk_label: The number-at-risk rows' title.
        risk_times: Times at which the patients still followed are counted.
        caption: The caption, unwrapped.
        width: The column; every figure takes ``FIGURE_WIDTH``, the double column.
        live: The survival panel's settings, which make this a figure with page controls; None
            for a figure the page does not redraw.
    """

    frame: pl.DataFrame
    columns: Mapping[str, tuple[str, str]]
    groups: tuple[str, ...]
    colours: Mapping[str, str]
    x_label: str
    y_label: str
    risk_label: str
    risk_times: tuple[float, ...]
    caption: str
    width: str = FIGURE_WIDTH
    live: Mapping[str, Any] | None = field(default=None)

    def __post_init__(self) -> None:
        """Refuse a contract the renderers could only draw wrongly."""
        if missing := [c for c in COLUMNS if c not in self.frame.columns]:
            raise ValueError(f"missing columns {missing}")
        if len(self.groups) > MAX_CATEGORIES:
            raise ValueError(f"{len(self.groups)} groups; at most {MAX_CATEGORIES}")
        if stray := set(self.frame["group"].unique()) - set(self.groups):
            raise ValueError(f"groups {sorted(stray)} are not named in groups")
        if uncoloured := [g for g in self.groups if g not in self.colours]:
            raise ValueError(f"groups {uncoloured} have no colour")

    def table(self) -> pl.DataFrame:
        """Every step of every curve."""
        return self.frame

    def table_columns(self) -> Mapping[str, tuple[str, str]]:
        """(unit, meaning) for every column of ``table()``."""
        return self.columns

    def page_table(self) -> pl.DataFrame:
        """Every step, rounded to what a reader reads.

        The plot carries the full values and the CSV every digit: unrounded, the rows of
        the breast cohort's four ER figures took 604,978 bytes of its page, rounded 466,356.
        """
        fractions = pl.col("survival", "lower", "upper").round(PAGE_DIGITS["fraction"])
        return self.frame.with_columns(pl.col("time").round(PAGE_DIGITS["time"]), fractions)

    def caption_text(self) -> str:
        """The caption, unwrapped."""
        return self.caption

    def controls(self) -> dict[str, Any] | None:
        """The survival panel's settings, or None for a figure the page does not redraw."""
        return {"kind": "survival", **self.live} if self.live is not None else None

    def inches(self) -> tuple[float, float]:
        """The figure files' size, which a page's own PDF and SVG of a changed view take."""
        return figure_inches(self.width, "km")

    def export_layers(self) -> dict[str, Any]:
        """Nothing: the page draws every line, band, tick and count the files draw."""
        return {}

    def _rows(self, group: str) -> pl.DataFrame:
        """One group's rows."""
        return self.frame.filter(pl.col("group") == group)

    def _titles(self) -> tuple[str, str]:
        """The two axis titles with their units."""
        return (axis_title(self.x_label, self.columns["time"][0]),
                axis_title(self.y_label, self.columns["survival"][0]))  # fmt: skip

    def at_risk(self, group: str) -> list[int]:
        """Patients of a group still followed at each of ``risk_times``.

        The at-risk count of the first distinct time at or after it, since no patient leaves
        between two distinct times; 0 past the last.
        """
        rows = self._rows(group)
        times, risk = rows["time"].to_list(), rows["at_risk"].to_list()
        return [next((n for t, n in zip(times, risk, strict=True) if t >= at), 0)
                for at in self.risk_times]  # fmt: skip

    def steps(self, group: str) -> dict[str, list[float]]:
        """The step line from 1 at time 0, as plotted: x, survival, lower and upper."""
        rows = self._rows(group).with_columns(pl.col("lower", "upper").cast(pl.Float64)
                                              .fill_null(float("nan")))  # fmt: skip
        out = {"x": [0.0, *rows["time"].to_list()], "y": [1.0, *rows["survival"].to_list()],
               "lower": [1.0, *rows["lower"].to_list()],
               "upper": [1.0, *rows["upper"].to_list()]}  # fmt: skip
        return out

    def band(self, group: str) -> tuple[list[float], list[float]]:
        """The band as one closed polygon of steps: the upper bound out, the lower back.

        Each bound holds its value until the next time, as the curve does; an undefined bound
        (survival 0) closes the band where it stops.
        """
        step = self.steps(group)
        x, upper, lower = step["x"], step["upper"], step["lower"]
        keep = [i for i in range(len(x)) if upper[i] == upper[i] and lower[i] == lower[i]]
        out_x, out_y, back_x, back_y = [], [], [], []
        for k, i in enumerate(keep):
            end = x[keep[k + 1]] if k + 1 < len(keep) else x[i]
            out_x += [x[i], end]
            out_y += [upper[i], upper[i]]
            back_x += [end, x[i]]
            back_y += [lower[i], lower[i]]
        return out_x + back_x[::-1], out_y + back_y[::-1]

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
        width, height = self.inches()
        fig = Figure(figsize=(width, height))
        grid = fig.add_gridspec(2, 1, height_ratios=(1 - KM["risk_share"], KM["risk_share"]))
        ax = fig.add_subplot(grid[0])
        risk = fig.add_subplot(grid[1], sharex=ax)
        drawn = [g for g in self.groups if self._rows(g).height]
        handles: list[Any] = []
        for g in drawn:
            step, rows, colour = self.steps(g), self._rows(g), self.colours[g]
            ax.fill_between(step["x"], step["lower"], step["upper"], step="post", color=colour,
                            alpha=BAND_ALPHA, linewidth=LINE_WIDTHS["band_edge"])  # fmt: skip
            ax.step(step["x"], step["y"], where="post", color=colour,
                    linewidth=LINE_WIDTHS["line"])  # fmt: skip
            censored = rows.filter(pl.col("censored") > 0)
            ax.plot(censored["time"].to_numpy(), censored["survival"].to_numpy(), linestyle="none",
                    marker="|", markersize=KM["censor_pt"], color=colour,
                    markeredgewidth=LINE_WIDTHS["line"])  # fmt: skip
            handles.append(Line2D([], [], color=colour, linewidth=LINE_WIDTHS["line"], label=g))
        x_title, y_title = self._titles()
        ax.set_ylim(0, 1.02)
        ax.set_ylabel(fit_title(y_title, height * (1 - KM["risk_share"])))
        ax.tick_params(labelbottom=False)
        for row, g in enumerate(drawn):
            for at, n in zip(self.risk_times, self.at_risk(g), strict=True):
                risk.text(at, row, f"{n:,}", ha="center", va="center", fontsize=SIZES["tick"])
        risk.set_yticks(range(len(drawn)), labels=drawn)
        risk.set_ylim(len(drawn) - 0.5, -0.5)
        risk.set_xticks(list(self.risk_times))
        risk.set_xlabel(x_title)
        risk.set_title(self.risk_label, fontsize=SIZES["tick"], loc="left")
        risk.tick_params(axis="y", length=0)
        sns.despine(fig=fig)
        sns.despine(ax=risk, left=True)
        place_caption(fig, self.caption, width, height, ax=ax, handles=handles)
        return fig

    def _traces(self, g: str, hover: str) -> list[Any]:
        """One group's band, line and censor ticks, in that order, as plotly draws them."""
        import plotly.graph_objects as go  # noqa: PLC0415 - heavy; imported after the profile

        step, rows, colour = self.steps(g), self._rows(g), self.colours[g]
        shown = rows.height > 0
        band_x, band_y = self.band(g)
        risk = [rows["at_risk"][0] if shown else 0, *rows["at_risk"].to_list()]
        custom = [list(v) for v in zip(step["lower"], step["upper"], risk, strict=True)]
        censored = rows.filter(pl.col("censored") > 0)
        line = {"color": colour, "width": LINE_WIDTHS["line"]}
        tick = {"symbol": "line-ns-open", "color": colour, "line": line,
                "size": KM["censor_pt"] * PAGE_PX_PER_PT}  # fmt: skip
        common = {"name": g, "legendgroup": g, "visible": shown}
        return [
            go.Scatter(x=band_x, y=band_y, fill="toself", fillcolor=colour, opacity=BAND_ALPHA,
                       line={"width": 0}, hoverinfo="skip", showlegend=False, **common),
            go.Scatter(x=step["x"], y=step["y"], mode="lines", line={**line, "shape": "hv"},
                       customdata=custom, hovertemplate=hover, **common),
            go.Scatter(x=censored["time"].to_list(), y=censored["survival"].to_list(),
                       mode="markers", marker=tick, showlegend=False, hoverinfo="skip",
                       **common),
        ]  # fmt: skip

    def interactive(self) -> Any:
        """Render the plotly figure.

        Per group a band, a line and its censor ticks, in that order, then one row of at-risk
        counts per group under a shared time axis: the trace order the page script redraws.

        Returns:
            A plotly Figure with the light template; the page switches themes.
        """
        import plotly.graph_objects as go  # noqa: PLC0415 - heavy; imported after the profile
        from plotly.subplots import make_subplots  # noqa: PLC0415 - heavy

        fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.08,
                            row_heights=[1 - KM["risk_share"], KM["risk_share"]])  # fmt: skip
        x_title, y_title = self._titles()
        hover = (f"%{{fullData.name}}<br>{x_title}: %{{x:.3~g}}<br>{y_title}: %{{y:.3f}}"
                 "<br>95 % CI %{customdata[0]:.3f} to %{customdata[1]:.3f}"
                 "<br>at risk %{customdata[2]:,d}<extra></extra>")  # fmt: skip
        for g in self.groups:
            for trace in self._traces(g, hover):
                fig.add_trace(trace, row=1, col=1)
        for g in self.groups:
            counts = [f"{n:,}" for n in self.at_risk(g)]
            row = go.Scatter(x=list(self.risk_times), y=[g] * len(self.risk_times), mode="text",
                             text=counts, showlegend=False, name=g, hoverinfo="skip",
                             visible=self._rows(g).height > 0,
                             textfont={"color": self.colours[g]})  # fmt: skip
            fig.add_trace(row, row=2, col=1)
        drawn = [g for g in self.groups if self._rows(g).height]
        fig.update_layout(
            template=plotly_templates()["white"], height=PLOTLY_HEIGHT_PX["km"],
            yaxis={"title": {"text": page_y_title(y_title, self.inches()[1])}, "range": [0, 1.02]},
            # A row's count is centred on its category: half a row of room at each end. Each
            # row takes its group's colour and the legend names it: written beside the rows,
            # the names left a 375 px phone a plot 20 px wide (2026-10-08).
            yaxis2={"type": "category", "categoryorder": "array",
                    "categoryarray": list(reversed(drawn)), "range": [-0.5, len(drawn) - 0.5],
                    "showgrid": False, "zeroline": False, "showticklabels": False,
                    "title": {"text": self.risk_label}},
            # The top plot shares the time axis and names it below the risk rows alone.
            xaxis={"title": {"text": ""}},
            xaxis2={"title": {"text": screen_title(x_title)}, "tickvals": list(self.risk_times)},
        )  # fmt: skip
        return fig
