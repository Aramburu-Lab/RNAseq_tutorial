"""Box plot: one value per sample, two groups, one or more features.

``static()`` and ``interactive()`` read the same ``BoxData``, so the two renderings cannot
disagree about the data or the statistics.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import polars as pl

from ctx.viz.annotate import Caption, compare_groups, place_caption, plotly_caption
from ctx.viz.style import (
    HEADROOM,
    LINE_WIDTHS,
    MARKER_SIZE,
    ROLES,
    SINGLE_COLUMN_MAX_FEATURES,
    SIZES,
    Theme,
    apply_matplotlib,
    figure_inches,
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
    column = "single" if len(features) <= SINGLE_COLUMN_MAX_FEATURES else "double"
    width, height = figure_inches(column, "box")
    # Figure, not pyplot: a GUI backend rounds the canvas to whole pixels (84.84 mm, not 85).
    fig = Figure(figsize=(width, height))
    ax = fig.subplots()
    reference, comparison = data.groups
    order = {"x": "feature", "y": "value", "hue": "group", "hue_order": list(data.groups)}
    palette = {reference: ROLES["reference"], comparison: ROLES["comparison"]}
    sns.boxplot(
        data.frame, **order, palette=palette, linecolor=theme.ink, saturation=1,
        linewidth=LINE_WIDTHS["box"], showfliers=False, ax=ax,
    )  # fmt: skip
    # No jitter: every run draws the same points, so the same data give the same file.
    sns.stripplot(
        data.frame, **order, palette=dict.fromkeys(data.groups, theme.ink), dodge=True,
        jitter=False, size=MARKER_SIZE, linewidth=LINE_WIDTHS["marker_edge"], legend=False, ax=ax,
    )  # fmt: skip
    top = float(data.frame["value"].to_numpy().max())
    for i, feature in enumerate(features):
        label = data.caption.label(feature)
        ax.text(i, top, label, ha="center", va="bottom", fontsize=SIZES["annotation"])
    low, high = ax.get_ylim()
    ax.set_ylim(top=high + (high - low) * HEADROOM)
    ax.set_xlabel(data.feature_label)
    ax.set_ylabel(f"{data.value_label} ({data.unit})")
    # Above the axes: inside, it collides with the p/r labels over the boxes.
    ax.legend(title=None, loc="lower left", bbox_to_anchor=(0, 1), ncols=2, borderaxespad=0)
    sns.despine(fig=fig)
    place_caption(fig, data.caption.text(), width, height)
    return fig


def interactive(data: BoxData) -> Any:
    """Render the plotly figure; hover shows every variable a reader would ask about.

    Args:
        data: The contract.

    Returns:
        A plotly Figure with the light template; the report page switches themes.
    """
    import plotly.graph_objects as go  # noqa: PLC0415 - heavy; imported after the runtime profile

    stats = {t.feature: t for t in data.caption.tests}
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
        custom = [[*h, *e] for h, e in zip(rows.select(data.hover).rows(), extra, strict=True)]
        k = len(data.hover)
        hover_lines = [f"{c}: %{{customdata[{i}]}}" for i, c in enumerate(data.hover)]
        fig.add_trace(
            go.Box(
                x=rows["feature"].to_list(),
                y=rows["value"].to_list(),
                name=group,
                marker_color=ROLES[role],
                boxpoints="all",
                jitter=0,
                pointpos=0,
                customdata=custom,
                hovertemplate="<br>".join(
                    [
                        f"{data.feature_label}: %{{x}}",
                        f"group: {group}",
                        *hover_lines,
                        f"{data.value_label}: %{{y}} {data.unit}",
                        f"n: %{{customdata[{k}]}}",
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
    caption, bottom = plotly_caption(data.caption.text(), figure_inches("double", "box")[0])
    fig.update_layout(
        template=plotly_templates()["white"],
        boxmode="group",
        xaxis_title=data.feature_label,
        yaxis_title=f"{data.value_label} ({data.unit})",
        annotations=[*labels, caption],
        margin={"b": bottom},
    )
    return fig
