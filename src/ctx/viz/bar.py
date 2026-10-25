"""Bars: one value per category and series, grouped or stacked, vertical or horizontal.

Draws the M1 composition and library sizes, p-value histograms, ORA and GSEA bars and call
counts. ``static()`` and ``interactive()`` read the same ``BarData``; with ``filters`` the
page script (``templates/_filter.js.j2``) redraws the bars and the median line from the rows
the reader picks.
"""

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
import polars as pl

from ctx.viz.annotate import files_room, fit_title, place_caption, refit_y_title
from ctx.viz.axes import (
    axis_title,
    hover_field,
    hover_rows,
    page_y_title,
    phone_legend,
    screen_title,
)
from ctx.viz.metrics import text_width, wrap_to
from ctx.viz.style import (
    BAR_LABEL_ROTATION,
    BAR_MARGIN_MM,
    BAR_ROW_MM,
    DASHES,
    FIGURE_WIDTH,
    LINE_SPACING,
    LINE_WIDTHS,
    MAX_CATEGORIES,
    MM_PER_INCH,
    PHONE_PLOT_PX,
    PLOTLY_BAR_BASE_PX,
    PLOTLY_BAR_ROW_PX,
    PLOTLY_HEIGHT_PX,
    PLOTLY_MARGIN_PX,
    PLOTLY_TICK_PX,
    SIZES,
    THEMES,
    THRESHOLD_DASH,
    Theme,
    apply_matplotlib,
    figure_inches,
    on_page,
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
        width: The column; every figure takes ``FIGURE_WIDTH``, the double column.
        names: Legend text per series level, where the level itself is terse.
        lines: (value, name) of further reference lines, drawn after ``hline``.
        median: Name of a line at the median bar length, empty for none. The page moves it
            to the median of the bars the filters show; the PDF and SVG show every bar.
        filters: Columns of ``frame`` the page offers as filters, one select each.
        labels: Words for a column on hover and on its filter, where its name is terse.
        summary_median: Whether the filter's live line prints the median bar length; off where
            that median is no number a reader asks for: bars stacked to 100 %, or several
            settings of one category pooled.
        live: Settings of a page control that redraws the bars, with its ``kind``, such as M6's
            live over-representation (GS.5); None for none. Not with ``filters``.
        drawable: Every category a ``live`` redraw may draw besides the built ones, such as
            every gene set M6's panel tests. With it, horizontal bars take a phone's layout as
            fixed ones do, at their threshold (D-105), and ``controls()`` carries what the
            redraw rebuilds it from (PK.4).
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
    width: str = FIGURE_WIDTH
    names: Mapping[str, str] | None = None
    lines: tuple[tuple[float, str], ...] = ()
    median: str = ""
    filters: tuple[str, ...] = ()
    labels: Mapping[str, str] | None = None
    summary_median: bool = True
    live: Mapping[str, Any] | None = None
    drawable: tuple[str, ...] = ()

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

    def _label(self, column: str) -> str:
        """A column's words on hover and on its filter: ``labels``' own, else its name."""
        return (self.labels or {}).get(column, column)

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

    def _height(self) -> float:
        """The plotly figure's height as built, in px: horizontal bars grow by the row."""
        if not self.horizontal:
            return PLOTLY_HEIGHT_PX["bar"]
        return max(PLOTLY_HEIGHT_PX["bar"],
                   PLOTLY_BAR_BASE_PX + PLOTLY_BAR_ROW_PX * len(self._categories()))  # fmt: skip

    def _labelled(self) -> bool:
        """Whether the category axis prints its labels."""
        return len(self._categories()) <= MAX_TICK_LABELS

    def _redraws_labels(self) -> bool:
        """Whether a live redraw rebuilds the phone's labels for categories of its own."""
        live = self.live is not None and bool(self.drawable)
        return self.horizontal and live and self._labelled()

    def _phone(self) -> dict[str, Any]:
        """What a live redraw rebuilds its phone's layout from, read by ``ctxOra.narrow``.

        Every category it may draw as ``category_lines`` gives it, the margins the threshold
        adds, and the height as built, which the redraw keeps.
        """
        categories = dict.fromkeys([*map(str, self._categories()), *self.drawable])
        return {"lines": category_lines(categories), "room": PHONE_PLOT_PX / 2,
                "base": PLOTLY_BAR_BASE_PX, "height": self._height(),
                "margins": PLOTLY_MARGIN_PX["l"] + PLOTLY_MARGIN_PX["r"]}  # fmt: skip

    def controls(self) -> dict[str, Any] | None:
        """What the page script needs to filter the bars, or None without ``filters``.

        ``start`` fills the summary line before the reader moves a filter, with the numbers
        the script computes for "All" (``tests/test_page_script.py`` compares the two).
        """
        if not self.filters:
            if self.live is None:
                return None
            return {**self.live, **({"phone": self._phone()} if self._redraws_labels() else {})}
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
            "filter_labels": {c: self._label(c) for c in self.filters},
            "summary_median": self.summary_median,
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
        if not self.horizontal:
            refit_y_title(fig, ax, self._value_title())
        return fig

    def trace_style(self, level: str) -> dict[str, Any]:
        """One series' bar trace as ``interactive`` draws it, points aside, whether drawn or not.

        ``ctxFdr.mark`` redraws M6's bars marked by FDR at another cutoff from these, a series
        the build left empty included (P2S.5).
        """
        schema = self.frame.schema
        value = hover_field("x" if self.horizontal else "y", schema[self.value])
        hover = [f"{self.category_label}: %{{{'y' if self.horizontal else 'x'}}}"]
        hover += [f"{self._value_title()}: {value}"]
        hover += [f"{self._label(c)}: {hover_field(f'customdata[{i}]', schema[c])}"
                  for i, c in enumerate(self.hover)]  # fmt: skip
        first = f"{self._label(self.series)}: {self._name(level)}"
        return {
            "name": self._name(level), "orientation": "h" if self.horizontal else "v",
            "marker": {"color": self.colours[level] or THEMES["white"].reference,
                       "line": {"width": LINE_WIDTHS["marker_edge"]}},
            "hovertemplate": "<br>".join([first, *hover]) + "<extra></extra>",
        }  # fmt: skip

    def interactive(self) -> Any:
        """Render the plotly figure; hover shows category, series, value and ``hover``.

        Returns:
            A plotly Figure with the light template; the page switches themes.
        """
        import plotly.graph_objects as go  # noqa: PLC0415 - heavy; imported after the profile

        light = THEMES["white"]
        fig = go.Figure()
        for level in self.colours:
            rows = self.frame.filter(pl.col(self.series) == level)
            if rows.height == 0:
                continue
            axes = (rows[self.value].to_list(), [str(c) for c in rows[self.category]])
            x, y = axes if self.horizontal else axes[::-1]
            fig.add_trace(go.Bar(x=x, y=y, customdata=hover_rows(rows, self.hover),
                                 **self.trace_style(level)))  # fmt: skip
        # Shapes, in rules() order, named in the legend as the static figure names them; the
        # page script moves the median's shape (controls()["median_shape"]).
        add = fig.add_vline if self.horizontal else fig.add_hline
        for at, name, (_, dash) in self.rules():
            add(at, line_dash=dash, line_color=light.reference, line_width=LINE_WIDTHS["line"],
                name=name, showlegend=True)  # fmt: skip
        categories = [str(c) for c in self._categories()]
        category_axis = {"title": {"text": self.category_label}, "categoryorder": "array",
                         "categoryarray": categories,
                         "showticklabels": self._labelled()}  # fmt: skip
        if self.horizontal:
            category_axis["autorange"] = "reversed"
        # The page script rewraps an x title for the screen; the export copies a y title.
        value_title = (screen_title(self._value_title()) if self.horizontal
                       else page_y_title(self._value_title(), files_room(self.static)))  # fmt: skip
        value_axis = {"title": {"text": value_title}}
        fig.update_layout(
            template=plotly_templates()["white"],
            height=self._height(),
            showlegend=len(fig.data) > 1 or bool(self.rules()),
            barmode="stack" if self.stacked else "group",
            xaxis=value_axis if self.horizontal else category_axis,
            yaxis=category_axis if self.horizontal else value_axis,
        )  # fmt: skip
        if self.horizontal:
            # The page's fit gives each category this much of the drawing area (D-119).
            fig.update_layout(meta={**(fig.layout.meta or {}), "row": PLOTLY_BAR_ROW_PX})
        self._phone_layout(fig, categories)
        return fig

    def _phone_layout(self, fig: Any, categories: list[str]) -> None:
        """Give ``fig`` a phone's layout, its legend's and its labels', where either needs one.

        A live redraw draws categories of its own: M1's recount, with no ``drawable``, keeps
        its labels as built; M6's over-representation is built as fixed bars are and rebuilds
        the phone's layout from ``controls()`` (PK.4, D-105).
        """
        legend = phone_legend(fig)
        fixed = self.horizontal and self.live is None and self._labelled()
        wraps = fixed or self._redraws_labels()
        labels = narrow_categories(categories, fig.layout.height) if wraps else None
        if labels is not None:
            # Named as ticks, so the narrow layout swaps wrapped ones in and these back.
            fig.update_layout(yaxis={"tickmode": "array", "tickvals": categories,
                                     "ticktext": categories})  # fmt: skip
        if legend is not None or labels is not None:
            found = [n for n in (legend, labels) if n is not None]
            narrow = {"below": max(n["below"] for n in found),
                      "relayout": {k: v for n in found for k, v in n["relayout"].items()},
                      "restyle": [pair for n in found for pair in n["restyle"]]}  # fmt: skip
            fig.update_layout(meta={**(fig.layout.meta or {}), "narrow": narrow})


def category_lines(categories: Iterable[str]) -> dict[str, list[Any]]:
    """Each label as a phone draws it: [wrapped to half a phone's plot, width, row height].

    Widths and heights in px; a row is the height its lines need, at least a bar's row.
    """
    size, room = on_page(SIZES["tick"]), PHONE_PLOT_PX / 2
    out: dict[str, list[Any]] = {}
    for c in categories:
        wrapped = "<br>".join(wrap_to(c, room, size))
        lines = wrapped.count("<br>") + 1
        out[c] = [wrapped, text_width(c, size),
                  max(PLOTLY_BAR_ROW_PX, lines * size * LINE_SPACING + PLOTLY_TICK_PX)]  # fmt: skip
    return out


def narrow_categories(
    categories: list[str], height: float, lines: Mapping[str, list[Any]] | None = None
) -> dict[str, Any] | None:
    """Horizontal bars' category labels wrapped to half a phone's plot, where one is wider.

    At 375 px M5's genes-called bars set "Neutrophils, PDAC against steady state, batch in the
    design" (530 px) left of a 317 px plot, and plotly's margin could not hold it: every label
    ran off the plot's left edge (E2b). The page script swaps the wrapped labels in where the
    plot is narrower than twice the widest label, so the bars keep half the plot, and each row
    takes the height its lines need. Every bar figure takes this one threshold (Rytis
    2026-10-21, D-105): M6's live over-representation, which PK.4 held to the widest label
    plus a pad, wraps at 1,280 px as its GSEA bars do. Its twin is ``ctxOra.narrow``.

    Args:
        categories: The labels, in axis order.
        height: The figure's height as built, in px.
        lines: ``category_lines`` of at least ``categories``; measured here when None.

    Returns:
        The narrow layout, or None when every label fits half a phone's plot.
    """
    lines = lines if lines is not None else category_lines(categories)
    widest = max((lines[c][1] for c in categories), default=0.0)
    if widest <= PHONE_PLOT_PX / 2:
        return None
    row = max(lines[c][2] for c in categories)
    return {
        "below": 2 * widest + PLOTLY_MARGIN_PX["l"] + PLOTLY_MARGIN_PX["r"],
        "relayout": {"yaxis.ticktext": [lines[c][0] for c in categories],
                     "height": max(height, PLOTLY_BAR_BASE_PX + row * len(categories))},
        "restyle": [],
    }  # fmt: skip
