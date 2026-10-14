"""Box plot: one value per sample, two groups, one or more features.

``static()`` and ``interactive()`` read the same ``BoxData``, so the two renderings cannot
disagree about the data or the statistics.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import polars as pl

from ctx.viz.annotate import (
    Caption,
    compare_groups,
    files_room,
    fit_title,
    place_caption,
    refit_y_title,
)
from ctx.viz.axes import axis_title, hover_field, hover_rows, hover_value, page_y_title
from ctx.viz.metrics import text_width
from ctx.viz.style import (
    BOX_GAP_PX,
    FIGURE_WIDTH,
    HEADROOM,
    LINE_SPACING,
    LINE_WIDTHS,
    MARKER_SIZE,
    PLOTLY_HEIGHT_PX,
    PLOTLY_MARGIN_PX,
    PLOTLY_MARKER_PX,
    POINT_ALPHA,
    ROLES,
    SIZES,
    Theme,
    apply_matplotlib,
    figure_inches,
    on_page,
    plotly_templates,
)

REQUIRED = ("feature", "group", "sample", "value")
# The numbers the figure prints, shipped beside the values they were computed from.
STAT_COLUMNS = {
    "n_group": ("samples", "Non-missing values of this feature in this row's group."),
    "u": ("pairs", "Mann-Whitney U of the comparison group against the reference group."),
    "p": ("dimensionless", "Two-sided Mann-Whitney U p-value, before adjustment."),
    "p_adj": ("dimensionless", "Benjamini-Hochberg adjusted p over the plotted features."),
    "r": ("dimensionless", "Rank-biserial correlation; > 0 when the comparison group is higher."),
}


@dataclass(frozen=True)
class BoxData:
    """The contract both renderers read.

    Attributes:
        frame: Long table with columns ``feature, group, sample, value`` plus any extras.
        columns: (unit, meaning) for every column of ``frame``; becomes the README dictionary.
        feature_label: What a feature is, e.g. ``"Gene"``; the x-axis title.
        value_label: What the value is, e.g. ``"Expression"``.
        unit: Unit of ``value``, taken from ``columns["value"]``; printed on the y axis.
        groups: (reference, comparison).
        hover: Columns shown on hover, besides feature, group and value.
        caption: The statistics, computed once.
    """

    frame: pl.DataFrame
    columns: Mapping[str, tuple[str, str]]
    feature_label: str
    value_label: str
    unit: str
    groups: tuple[str, str]
    hover: tuple[str, ...]
    caption: Caption

    def table(self) -> pl.DataFrame:
        """The values and the statistics the figure prints; see ``table()``."""
        return table(self)

    def table_columns(self) -> dict[str, tuple[str, str]]:
        """(unit, meaning) for every column of ``table()``."""
        return {**self.columns, **STAT_COLUMNS}

    def page_table(self) -> pl.DataFrame:
        """Every row: a box plot draws every sample."""
        return table(self)

    def static(self, theme: Theme) -> Any:
        """The matplotlib figure for one theme; see ``static()``."""
        return static(self, theme)

    def interactive(self) -> Any:
        """The plotly figure; see ``interactive()``."""
        return interactive(self)

    def caption_text(self) -> str:
        """The caption, unwrapped."""
        return self.caption.text()


def box_data(
    frame: pl.DataFrame,
    *,
    columns: Mapping[str, tuple[str, str]],
    feature_label: str,
    value_label: str,
    groups: tuple[str, str],
    hover: tuple[str, ...] = ("sample",),
) -> BoxData:
    """Check the frame against the contract and compute the statistics once.

    Args:
        frame: Long table; see ``BoxData.frame``.
        columns: (unit, meaning) per column.
        feature_label: x-axis title.
        value_label: y-axis title; the unit is read from ``columns["value"]``, stated once.
        groups: (reference, comparison).
        hover: Extra columns to show on hover.

    Returns:
        The figure's data contract.

    Raises:
        ValueError: If a required or hover column is missing.
    """
    if missing := [c for c in (*REQUIRED, *hover) if c not in frame.columns]:
        raise ValueError(f"missing columns {missing}")
    if "value" not in columns:
        raise ValueError("columns has no (unit, meaning) for value")
    caption = compare_groups(frame, feature="feature", group="group", value="value", groups=groups)
    unit = columns["value"][0]
    return BoxData(frame, columns, feature_label, value_label, unit, groups, hover, caption)


def table(data: BoxData) -> pl.DataFrame:
    """The values plus, on every row, the statistics the figure prints for that feature.

    Args:
        data: The contract.

    Returns:
        ``data.frame`` with the ``STAT_COLUMNS`` appended, rows in their original order.
    """
    stats = pl.DataFrame(
        [
            {"feature": t.feature, "n_ref": t.n_reference, "n_cmp": t.n_comparison,
             "u": t.u, "p": t.p, "p_adj": t.p_adj, "r": t.r}
            for t in data.caption.tests
        ]
    )  # fmt: skip
    n_group = pl.when(pl.col("group") == data.groups[0]).then("n_ref").otherwise("n_cmp")
    joined = data.frame.join(stats, on="feature", how="left", maintain_order="left")
    return joined.with_columns(n_group.alias("n_group")).select(*data.frame.columns, *STAT_COLUMNS)


def static(data: BoxData, theme: Theme) -> Any:
    """Render the matplotlib figure for one theme.

    Args:
        data: The contract.
        theme: One of ``style.THEMES``.

    Returns:
        A matplotlib Figure, not attached to pyplot; the caller saves it.
    """
    import seaborn as sns  # noqa: PLC0415 - heavy; imported after the runtime profile
    from matplotlib.figure import Figure  # noqa: PLC0415 - heavy; imported after the profile

    apply_matplotlib(theme)
    features = data.frame["feature"].unique(maintain_order=True).to_list()
    width, height = _inches()
    # Figure, not pyplot: a GUI backend rounds the canvas to whole pixels (84.84 mm, not 85).
    fig = Figure(figsize=(width, height))
    ax = fig.subplots()
    # matplotlib's own boxplot, not seaborn's: seaborn 0.13.2 passes bxp the vert argument
    # matplotlib 3.11 deprecates. Geometry as seaborn dodges two groups: boxes 0.4 wide at
    # +-0.2 around each feature, caps half a box.
    ink = {"color": theme.ink, "linewidth": LINE_WIDTHS["box"]}
    for offset, role, group in zip((-0.2, 0.2), ("reference", "comparison"), data.groups,
                                   strict=True):  # fmt: skip
        values = [_values(data.frame, f, group) for f in features]
        positions = [i + offset for i in range(len(features))]
        ax.boxplot(
            values, positions=positions, widths=0.4, capwidths=0.2, orientation="vertical",
            patch_artist=True, showfliers=False, manage_ticks=False, label=group,
            boxprops={"facecolor": ROLES[role], "edgecolor": theme.ink,
                      "linewidth": LINE_WIDTHS["box"]},
            medianprops={**ink, "solid_capstyle": "butt"},
            whiskerprops={**ink, "solid_capstyle": "butt"}, capprops=ink,
        )  # fmt: skip
        # No jitter: every run draws the same points, so the same data give the same file.
        for x, v in zip(positions, values, strict=True):
            ax.scatter([x] * len(v), v, s=MARKER_SIZE**2, color=theme.ink, alpha=POINT_ALPHA,
                       linewidths=LINE_WIDTHS["marker_edge"], zorder=3)  # fmt: skip
    ax.set_xticks(range(len(features)), features)
    ax.set_xlim(-0.5, len(features) - 0.5)
    top = float(data.frame["value"].to_numpy().max())
    for i, feature in enumerate(features):
        label = data.caption.label(feature)
        ax.text(i, top, label, ha="center", va="bottom", fontsize=SIZES["annotation"])
    low, high = ax.get_ylim()
    ax.set_ylim(top=high + (high - low) * HEADROOM)
    ax.set_xlabel(data.feature_label)
    ax.set_ylabel(fit_title(axis_title(data.value_label, data.unit), height))
    # Above the axes: inside, it collides with the p/r labels over the boxes.
    ax.legend(title=None, loc="lower left", bbox_to_anchor=(0, 1), ncols=2, borderaxespad=0)
    sns.despine(fig=fig)
    place_caption(fig, data.caption.text(), width, height)
    refit_y_title(fig, ax, axis_title(data.value_label, data.unit))
    return fig


def _inches() -> tuple[float, float]:
    """The figure files' size, whose height also sets where the page breaks its y title."""
    return figure_inches(FIGURE_WIDTH, "box")


def _values(frame: pl.DataFrame, feature: str, group: str) -> list[float]:
    """One feature's measured values in one group; a missing value is not drawn."""
    rows = frame.filter((pl.col("feature") == feature) & (pl.col("group") == group))
    return rows["value"].drop_nulls().drop_nans().to_list()


def interactive(data: BoxData) -> Any:
    """Render the plotly figure; hover shows every variable a reader would ask about.

    Args:
        data: The contract.

    Returns:
        A plotly Figure with the light template; the report page switches themes.
    """
    import plotly.graph_objects as go  # noqa: PLC0415 - heavy; imported after the runtime profile

    stats = {t.feature: t for t in data.caption.tests}
    schema, title = data.frame.schema, axis_title(data.value_label, data.unit)
    fig = go.Figure()
    for role, group in zip(("reference", "comparison"), data.groups, strict=True):
        rows = data.frame.filter(pl.col("group") == group)
        extra = [
            [
                stats[f].n_reference if role == "reference" else stats[f].n_comparison,
                stats[f].p_adj,
                stats[f].r,
            ]
            for f in rows["feature"]
        ]
        custom = [[*h, *map(hover_value, e)]
                  for h, e in zip(hover_rows(rows, data.hover), extra, strict=True)]  # fmt: skip
        k = len(data.hover)
        hover_lines = [f"{c}: {hover_field(f'customdata[{i}]', schema[c])}"
                       for i, c in enumerate(data.hover)]  # fmt: skip
        fig.add_trace(
            go.Box(
                x=rows["feature"].to_list(),
                y=rows["value"].to_list(),
                name=group,
                marker={"color": ROLES[role], "opacity": POINT_ALPHA, "size": PLOTLY_MARKER_PX},
                boxpoints="all",
                jitter=0,
                pointpos=0,
                customdata=custom,
                hovertemplate="<br>".join(
                    [
                        f"{data.feature_label}: %{{x}}",
                        f"group: {group}",
                        *hover_lines,
                        f"{title}: {hover_field('y', schema['value'])}",
                        f"n: {hover_field(f'customdata[{k}]', pl.Int64())}",
                        f"p_adj: %{{customdata[{k + 1}]:.2g}}",
                        f"r: %{{customdata[{k + 2}]:.2f}}",
                    ]
                )
                + "<extra></extra>",
            )
        )
    top = float(data.frame["value"].to_numpy().max())
    labels = [
        {
            "x": f,
            "y": top,
            "text": data.caption.label(f).replace("\n", "<br>"),
            "showarrow": False,
            "yanchor": "bottom",
        }
        for f in stats
    ]
    low, high = float(data.frame["value"].to_numpy().min()), top
    span = (high - low) or 1.0
    fig.update_layout(
        template=plotly_templates()["white"],
        height=PLOTLY_HEIGHT_PX["box"],
        boxmode="group",
        xaxis_title=data.feature_label,
        # Room above the tallest box for the p and r labels, as the static figure leaves.
        yaxis={"title": {"text": page_y_title(title, files_room(lambda t: static(data, t)))},
               "range": [low - 0.05 * span, high + HEADROOM * span]},
        annotations=labels,
    )  # fmt: skip
    y_range = (low - 0.05 * span, high + HEADROOM * span)
    wide, narrow = least_width(data, y_range), narrow_width(data, y_range)
    fig.update_layout(meta={"least_width": narrow,
                            "narrow": _narrow(data, sum(wide[0]), len(labels))})  # fmt: skip
    return fig


def _narrow(data: BoxData, below: float, n_labels: int) -> dict[str, Any]:
    """A phone's layout, swapped in below ``below`` px: each feature's p and r under its name.

    At 375 px the kidney page's four genes needed 480 px of upright labels and the box
    scrolled inside its card (C3.4). A feature's name and its p and r line go into one tick
    label turned a quarter, whose slot needs only two lines' height, and the labels over the
    boxes hide; plotly's automargin gives the turned labels their room below the axis.
    """
    features = data.frame["feature"].unique(maintain_order=True).to_list()
    ticks = [f"{f}<br>{data.caption.label(f, 'label_line')}" for f in features]
    relayout: dict[str, Any] = {"xaxis.tickangle": -90, "xaxis.tickmode": "array",
                                "xaxis.tickvals": features, "xaxis.ticktext": ticks}  # fmt: skip
    relayout |= {f"annotations[{i}].visible": False for i in range(n_labels)}
    return {"below": below, "relayout": relayout, "restyle": []}


def narrow_width(data: BoxData, y_range: tuple[float, float]) -> list[list[float]]:
    """The narrowest plot the phone's layout (``_narrow``) keeps apart.

    Each feature's turned label is as many lines tall as it has, with ``BOX_GAP_PX`` between,
    then the y axis's ticks and title. The page draws the box no narrower and scrolls it
    sideways only below that.

    Returns:
        One [scaled, fixed] px pair, as ``least_width``.
    """
    tick = on_page(SIZES["tick"])
    features = data.frame["feature"].unique(maintain_order=True).to_list()
    lines = max(2 + f.count("\n") for f in features)
    y_ticks = max(text_width(f"{v:.0f}", tick) for v in y_range)
    slot = lines * tick * LINE_SPACING + BOX_GAP_PX
    scaled = len(features) * slot + y_ticks + on_page(SIZES["axis_label"])
    return [[scaled, float(PLOTLY_MARGIN_PX["l"] + PLOTLY_MARGIN_PX["r"])]]


def least_width(data: BoxData, y_range: tuple[float, float]) -> list[list[float]]:
    """The narrowest plot that keeps each feature's tick label and p and r label apart.

    plotly gives every feature an equal slot, so each takes the widest label any needs plus
    ``BOX_GAP_PX``; the left margin holds the y tick labels and the axis title. At 375 px a
    phone squeezed four TCGA kidney genes into 63 px slots: NDUFA4L2, 101 px at the page's
    tick size, ran over VEGFA, and the p labels over each other (2026-10-09). The page draws
    the box's upright layout no narrower than this, and a phone's below it (``_narrow``).

    Returns:
        One [scaled, fixed] px pair: at zoom s the plot needs s * scaled + fixed.
    """
    tick, note = on_page(SIZES["tick"]), on_page(SIZES["annotation"])
    features = data.frame["feature"].unique(maintain_order=True).to_list()
    widest = max(
        max(
            text_width(f, tick),
            *(text_width(line, note) for line in data.caption.label(f).split("\n")),
        )
        for f in features
    )
    y_ticks = max(text_width(f"{v:.0f}", tick) for v in y_range)
    scaled = len(features) * (widest + BOX_GAP_PX) + y_ticks + on_page(SIZES["axis_label"])
    return [[scaled, float(PLOTLY_MARGIN_PX["l"] + PLOTLY_MARGIN_PX["r"])]]
