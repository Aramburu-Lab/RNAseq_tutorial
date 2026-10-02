"""Bars: one value per category and series, grouped or stacked, vertical or horizontal.

Draws the M1 composition and library sizes, p-value histograms, ORA and GSEA bars and call
counts. ``static()`` and ``interactive()`` read the same ``BarData``; with ``filters`` the
page script (``templates/_filter.js.j2``) redraws the bars and the median line from the rows
the reader picks.
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
import polars as pl

from ctx.viz.annotate import fit_title, place_caption
from ctx.viz.axes import axis_title
from ctx.viz.style import (
    BAR_LABEL_ROTATION,
    BAR_MARGIN_MM,
    BAR_ROW_MM,
    DASHES,
    LINE_WIDTHS,
    MAX_CATEGORIES,
    MM_PER_INCH,
    PLOTLY_BAR_BASE_PX,
    PLOTLY_BAR_ROW_PX,
    PLOTLY_HEIGHT_PX,
    THEMES,
    THRESHOLD_DASH,
    Theme,
    apply_matplotlib,
    figure_inches,
    plotly_templates,
)

# Above this many categories the category axis carries no tick labels; hover names each bar.
MAX_TICK_LABELS = 30
BAR_SPAN = 0.8
# Reference lines take these dashes in turn, the first the threshold dash every figure uses;
# the median line is solid, so no two lines of one figure share a style.
REFERENCE_DASHES = (THRESHOLD_DASH, *DASHES[2:])
MEDIAN_DASH = DASHES[0]
# The median the page prints, to this many significant digits (ctxFilter.number).
SIGNIFICANT = 3


def significant(value: float) -> str:
    """A number to three significant digits with thousands separators, as ctxFilter.number.

    Rounded half up by the same arithmetic as the page script, so both print the same digits.
    """
    decimals = max(0, SIGNIFICANT - 1 - (math.floor(math.log10(abs(value))) if value else 0))
    scale = 10**decimals
    return f"{math.floor(value * scale + 0.5) / scale:,.{decimals}f}"


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
        names: Legend text per series level, where the level itself is terse.
        lines: (value, name) of further reference lines, drawn after ``hline``.
        median: Name of a line at the median bar length, empty for none. The page moves it
            to the median of the bars the filters show; the PDF and SVG show every bar.
        filters: Columns of ``frame`` the page offers as filters, one select each.
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
    names: Mapping[str, str] | None = None
    lines: tuple[tuple[float, str], ...] = ()
    median: str = ""
    filters: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Refuse a contract the renderers could only draw wrongly."""
        if self.frame.is_empty():
            raise ValueError("no bars to draw: the table is empty")
        if len(self._references()) > len(REFERENCE_DASHES):
            raise ValueError(f"at most {len(REFERENCE_DASHES)} reference lines")
        named = [self.category, self.value, self.series, *self.hover, *self.filters]
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
        return axis_title(self.value_label, self.columns[self.value][0])

    def _name(self, level: str) -> str:
        """Legend text of a series level."""
        return (self.names or {}).get(level, str(level))

    def _references(self) -> list[tuple[float, str]]:
        """``hline``, then ``lines``."""
        return [*([self.hline] if self.hline else []), *self.lines]

    def lengths(self) -> np.ndarray:
        """What the median is taken over: each category's total when stacked, else each bar."""
        if not self.stacked:
            return self.frame[self.value].to_numpy().astype(float)
        totals = self.frame.group_by(self.category, maintain_order=True).agg(
            pl.col(self.value).sum()
        )
        return totals[self.value].to_numpy().astype(float)

    def rules(self) -> list[tuple[float, str, tuple[str, str]]]:
        """Every line across the value axis: (value, name, (matplotlib, plotly) dash).

        The reference lines, then the median; both renderers draw exactly these.
        """
        drawn: list[tuple[float, str, tuple[str, str]]] = [
            (v, name, dash)
            for (v, name), dash in zip(self._references(), REFERENCE_DASHES, strict=False)
        ]
        if self.median:
            drawn.append((float(np.median(self.lengths())), self.median, MEDIAN_DASH))
        return drawn

    def controls(self) -> dict[str, Any] | None:
        """What the page script needs to filter the bars, or None without ``filters``.

        ``start`` fills the summary line before the reader moves a filter, with the numbers
        the script computes for "All" (``tests/test_page_script.py`` compares the two).
        """
        if not self.filters:
            return None
        unit = self.columns[self.value][0]
        unit = "" if unit == "none" else unit
        total = len(self.lengths())
        median = f"{significant(float(np.median(self.lengths())))} {unit}".strip()
        present = set(self.frame[self.series])
        return {
            "kind": "filter", "category": self.category, "value": self.value,
            "series": self.series, "stacked": self.stacked, "horizontal": self.horizontal,
            "traces": [str(level) for level in self.colours if level in present],
            "labels": [str(c) for c in self.frame[self.category]], "hover": list(self.hover),
            "filters": {c: sorted({str(v) for v in self.frame[c]}) for c in self.filters},
            "median_shape": len(self._references()) if self.median else None,
            "max_tick_labels": MAX_TICK_LABELS, "unit": unit, "total": total,
            "start": {"shown": f"{total:,}", "total": f"{total:,}", "median": median},
        }  # fmt: skip

    def inches(self) -> tuple[float, float]:
        """The figure's size, which a page's own PDF and SVG of a filtered view take too.

        Horizontal bars grow by a row per category.
        """
        width, height = figure_inches(self.width, "bar")
        if self.horizontal:
            rows = (len(self._categories()) * BAR_ROW_MM + BAR_MARGIN_MM) / MM_PER_INCH
            height = max(height, rows)
        return width, height

    def export_layers(self) -> dict[str, Any]:
        """Nothing: the page draws every bar and line the figure files draw."""
        return {}

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
        width, height = self.inches()
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
            Patch(color=self.colours[level] or theme.reference, label=self._name(level))
            for level in present
        ]
        line = ax.axvline if self.horizontal else ax.axhline
        for at, name, (dash, _) in self.rules():
            style: dict[str, Any] = {"color": theme.reference, "linestyle": dash,
                                     "linewidth": LINE_WIDTHS["line"]}  # fmt: skip
            line(at, **style)
            handles.append(Line2D([], [], label=name, **style))
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
        # A count prints whole, a measurement to three significant digits.
        digits = ":,d" if self.frame.schema[self.value].is_integer() else ":.3~g"
        value = f"%{{{'x' if self.horizontal else 'y'}{digits}}}"
        hover = [f"{self.category_label}: %{{{'y' if self.horizontal else 'x'}}}"]
        hover += [f"{self._value_title()}: {value}"]
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
                    x=x, y=y, name=self._name(level), orientation="h" if self.horizontal else "v",
                    marker={"color": self.colours[level] or light.reference,
                            "line": {"width": LINE_WIDTHS["marker_edge"]}},
                    customdata=rows.select(self.hover).rows(),
                    hovertemplate="<br>".join([f"{self.series}: {self._name(level)}", *hover])
                    + "<extra></extra>",
                )
            )  # fmt: skip
        # Shapes, in rules() order, named in the legend as the static figure names them; the
        # page script moves the median's shape (controls()["median_shape"]).
        add = fig.add_vline if self.horizontal else fig.add_hline
        for at, name, (_, dash) in self.rules():
            add(at, line_dash=dash, line_color=light.reference, line_width=LINE_WIDTHS["line"],
                name=name, showlegend=True)  # fmt: skip
        categories = [str(c) for c in self._categories()]
        category_axis = {"title": {"text": self.category_label}, "categoryorder": "array",
                         "categoryarray": categories,
                         "showticklabels": len(categories) <= MAX_TICK_LABELS}  # fmt: skip
        if self.horizontal:
            category_axis["autorange"] = "reversed"
        value_axis = {"title": {"text": self._value_title()}}
        rows_px = PLOTLY_BAR_BASE_PX + PLOTLY_BAR_ROW_PX * len(categories)
        fig.update_layout(
            template=plotly_templates()["white"],
            height=max(PLOTLY_HEIGHT_PX["bar"], rows_px) if self.horizontal
            else PLOTLY_HEIGHT_PX["bar"],
            showlegend=len(fig.data) > 1 or bool(self.rules()),
            barmode="stack" if self.stacked else "group",
            xaxis=value_axis if self.horizontal else category_axis,
            yaxis=category_axis if self.horizontal else value_axis,
        )  # fmt: skip
        return fig
