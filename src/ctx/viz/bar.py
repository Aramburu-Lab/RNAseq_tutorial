"""Bars: one value per category and series, grouped or stacked, vertical or horizontal.

Draws the M1 composition and library sizes, p-value histograms, ORA and GSEA bars and call
counts. ``static()`` and ``interactive()`` read the same ``BarData``.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
import polars as pl

from ctx.viz.annotate import fit_title, place_caption, plotly_caption
from ctx.viz.style import (
    BAR_LABEL_ROTATION,
    BAR_MARGIN_MM,
    BAR_ROW_MM,
    LINE_WIDTHS,
    MAX_CATEGORIES,
    MM_PER_INCH,
    THEMES,
    Theme,
    apply_matplotlib,
    figure_inches,
    plotly_templates,
)

# Above this many categories the category axis carries no tick labels; hover names each bar.
MAX_TICK_LABELS = 30
BAR_SPAN = 0.8


@dataclass(frozen=True)
class BarData:
    """The contract both renderers read.

    Attributes:
        frame: One row per bar: category, series and value.
        columns: (unit, meaning) for every column of ``frame``.
        category: Column naming the bar's position on the category axis.
        value: Column with the bar's length.
        category_label: Title of the category axis.
        value_label: Title of the value axis, without the unit.
        series: Column whose levels colour the bars.
        colours: Level to hex colour, in legend order; None is the theme's reference grey.
        caption: The statistics sentence, unwrapped.
        horizontal: Bars run left to right, first category on top.
        stacked: Series stack within a category instead of standing side by side.
        hline: (value, name) of a reference line across the value axis, or None.
        hover: Columns shown on hover, besides category, series and value.
        width: ``"single"`` or ``"double"`` column.
    """

    frame: pl.DataFrame
    columns: Mapping[str, tuple[str, str]]
    category: str
    value: str
    category_label: str
    value_label: str
    series: str
    colours: Mapping[str, str | None]
    caption: str
    horizontal: bool = False
    stacked: bool = False
    hline: tuple[float, str] | None = None
    hover: tuple[str, ...] = ()
    width: str = "single"

    def __post_init__(self) -> None:
        """Refuse a contract the renderers could only draw wrongly."""
        if self.frame.is_empty():
            raise ValueError("no bars to draw: the table is empty")
        named = [self.category, self.value, self.series, *self.hover]
        if missing := [c for c in named if c not in self.frame.columns]:
            raise ValueError(f"missing columns {missing}")
        if len(self.colours) > MAX_CATEGORIES:
            raise ValueError(f"{len(self.colours)} series; at most {MAX_CATEGORIES}")
        if stray := set(self.frame[self.series].unique()) - set(self.colours):
            raise ValueError(f"{self.series}: levels {sorted(stray)} have no colour")
        if self.frame.select(self.category, self.series).is_duplicated().any():
            raise ValueError("more than one row for a (category, series) pair")

    def table(self) -> pl.DataFrame:
        """Every bar."""
        return self.frame

    def table_columns(self) -> Mapping[str, tuple[str, str]]:
        """(unit, meaning) for every column of ``table()``."""
        return self.columns

    def page_table(self) -> pl.DataFrame:
        """Every bar."""
        return self.frame

    def caption_text(self) -> str:
        """The caption, unwrapped."""
        return self.caption

    def _categories(self) -> list[Any]:
        """Categories in first-appearance order."""
        return self.frame[self.category].unique(maintain_order=True).to_list()

    def _values(self, level: str) -> np.ndarray:
        """One series' values in category order; 0 where the series has no bar."""
        rows = dict(self.frame.filter(pl.col(self.series) == level).select(
            self.category, self.value).rows())  # fmt: skip
        return np.array([rows.get(c, 0.0) for c in self._categories()], dtype=float)

    def _value_title(self) -> str:
        """Value-axis title with the unit from the column dictionary."""
        return f"{self.value_label} ({self.columns[self.value][0]})"

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
        from matplotlib.patches import Patch  # noqa: PLC0415 - heavy

        apply_matplotlib(theme)
        categories = self._categories()
        width, height = figure_inches(self.width, "bar")
        if self.horizontal:
            rows = (len(categories) * BAR_ROW_MM + BAR_MARGIN_MM) / MM_PER_INCH
            height = max(height, rows)
        fig = Figure(figsize=(width, height))
        ax = fig.subplots()
        positions = np.arange(len(categories), dtype=float)
        present = [level for level in self.colours if level in set(self.frame[self.series])]
        thickness = BAR_SPAN if self.stacked else BAR_SPAN / len(present)
        base = np.zeros(len(categories))
        # Any: bar and barh name their base edge differently, which one call cannot type.
        draw: Any = ax.barh if self.horizontal else ax.bar
        for i, level in enumerate(present):
            values = self._values(level)
            offset = 0.0 if self.stacked else (i - (len(present) - 1) / 2) * thickness
            edge = {"left" if self.horizontal else "bottom": base if self.stacked else None}
            draw(
                positions + offset, values, thickness, **edge,
                color=self.colours[level] or theme.reference,
                linewidth=LINE_WIDTHS["marker_edge"],
            )  # fmt: skip
            base = base + values
        handles: list[Any] = [
            Patch(color=self.colours[level] or theme.reference, label=str(level))
            for level in present
        ]
        if self.hline:
            line = ax.axvline if self.horizontal else ax.axhline
            line(self.hline[0], color=theme.ink, linestyle="--", linewidth=LINE_WIDTHS["line"])
            handles.append(Line2D([], [], color=theme.ink, linestyle="--", label=self.hline[1]))
        ticks = (ax.set_yticks, ax.set_ylabel, ax.set_xlabel) if self.horizontal else (
            ax.set_xticks, ax.set_xlabel, ax.set_ylabel)  # fmt: skip
        shown = len(categories) <= MAX_TICK_LABELS
        rotation = 0.0 if self.horizontal else BAR_LABEL_ROTATION
        labels = [str(c) for c in categories]
        ticks[0](positions if shown else [], labels if shown else [], rotation=rotation)
        ticks[1](self.category_label)
        ticks[2](self._value_title() if self.horizontal else fit_title(self._value_title(), height))
        if self.horizontal:
            ax.invert_yaxis()
        sns.despine(fig=fig)
        place_caption(fig, self.caption, width, height, ax=ax, handles=handles)
        return fig

    def interactive(self) -> Any:
        """Render the plotly figure; hover shows category, series, value and ``hover``.

        Returns:
            A plotly Figure with the light template; the page switches themes.
        """
        import plotly.graph_objects as go  # noqa: PLC0415 - heavy; imported after the profile

        light = THEMES["white"]
        hover = [f"{self.category_label}: %{{{'y' if self.horizontal else 'x'}}}"]
        hover += [f"{self.value_label}: %{{{'x' if self.horizontal else 'y'}}}"]
        hover += [f"{c}: %{{customdata[{i}]}}" for i, c in enumerate(self.hover)]
        fig = go.Figure()
        for level in self.colours:
            rows = self.frame.filter(pl.col(self.series) == level)
            if rows.height == 0:
                continue
            axes = (rows[self.value].to_list(), [str(c) for c in rows[self.category]])
            x, y = axes if self.horizontal else axes[::-1]
            fig.add_trace(
                go.Bar(
                    x=x, y=y, name=str(level), orientation="h" if self.horizontal else "v",
                    marker_color=self.colours[level] or light.reference,
                    customdata=rows.select(self.hover).rows(),
                    hovertemplate="<br>".join([f"{self.series}: {level}", *hover])
                    + "<extra></extra>",
                )
            )  # fmt: skip
        if self.hline:
            add = fig.add_vline if self.horizontal else fig.add_hline
            add(self.hline[0], line_dash="dash", line_color=light.reference,
                annotation_text=self.hline[1])  # fmt: skip
        caption, bottom = plotly_caption(self.caption, figure_inches("double", "bar")[0])
        category_axis = {"title": self.category_label, "categoryorder": "array",
                         "categoryarray": [str(c) for c in self._categories()],
                         "showticklabels": len(self._categories()) <= MAX_TICK_LABELS}  # fmt: skip
        if self.horizontal:
            category_axis["autorange"] = "reversed"
        value_axis = {"title": self._value_title()}
        fig.update_layout(
            template=plotly_templates()["white"],
            barmode="stack" if self.stacked else "group",
            xaxis=value_axis if self.horizontal else category_axis,
            yaxis=category_axis if self.horizontal else value_axis,
            annotations=[caption],
            margin={"b": bottom},
        )
        return fig
