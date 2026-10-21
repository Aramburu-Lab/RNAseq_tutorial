"""Box plot: one value per sample, two groups, one or more features.

``static()`` and ``interactive()`` read the same ``BoxData``, so the two renderings cannot
disagree about the data or the statistics. Where a contrast is paired, ``pairs`` names the
unit measured in both groups, a donor, and both draw a line from its reference value to its
comparison value (``paired_lines``, E.1e). A feature's name may hold line breaks: each line is
drawn and measured as a line (``_lines``).
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
    STACK_GAP_PT,
    THEMES,
    Theme,
    apply_matplotlib,
    figure_inches,
    on_page,
    plotly_templates,
)

REQUIRED = ("feature", "group", "sample", "value")
# Each box's centre, left of and right of its feature's position, as seaborn dodges two
# groups; the page's boxes sit at the same places on the same numeric axis.
OFFSETS = (-0.2, 0.2)
BOX_WIDTH = 0.4
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
        pairs: The column naming the unit measured in both groups, a donor, where the
            contrast is paired; empty draws no lines.
    """

    frame: pl.DataFrame
    columns: Mapping[str, tuple[str, str]]
    feature_label: str
    value_label: str
    unit: str
    groups: tuple[str, str]
    hover: tuple[str, ...]
    caption: Caption
    pairs: str = ""

    def __post_init__(self) -> None:
        """Refuse a pairs column the frame does not hold: no line could be drawn from it."""
        if self.pairs and self.pairs not in self.frame.columns:
            raise ValueError(f"missing pairs column {self.pairs!r}")

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

    def inches(self) -> tuple[float, float]:
        """The figure files' size, which a page's own PDF and SVG of a changed view take too."""
        return _inches()

    def export_layers(self) -> dict[str, Any]:
        """Nothing: the page draws every box, point and label the figure files draw."""
        return {}


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
        The figure's data contract; a paired contrast's caller sets ``pairs`` with
        ``dataclasses.replace``.

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


def paired_lines(data: BoxData) -> dict[str, list[tuple[str, float, float]]]:
    """Each feature's (unit, reference value, comparison value), one line per paired unit.

    A unit is drawn where it holds exactly one measured value in each group of that feature: a
    second sample in one group leaves no one line to draw, so it draws none. Units in the
    order they first appear; no lines where ``pairs`` is empty.
    """
    if not data.pairs:
        return {}
    reference, comparison = data.groups
    measured = data.frame.filter(pl.col("value").is_not_null() & pl.col("value").is_not_nan()
                                 & pl.col(data.pairs).is_not_null())  # fmt: skip
    once = measured.filter(pl.len().over("feature", data.pairs, "group") == 1)
    wide = once.pivot("group", index=["feature", data.pairs], values="value",
                      aggregate_function="first", maintain_order=True)  # fmt: skip
    if reference not in wide.columns or comparison not in wide.columns:
        return {}
    both = wide.drop_nulls([reference, comparison])
    out: dict[str, list[tuple[str, float, float]]] = {}
    for feature, unit, a, b in both.select("feature", data.pairs, reference, comparison).rows():
        out.setdefault(feature, []).append((str(unit), float(a), float(b)))
    return out


def _lines(name: str) -> list[str]:
    """A feature's name as its tick prints it, one entry per line."""
    return name.split("\n")


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


@dataclass(frozen=True)
class Stack:
    """How many rows the p and r labels and the feature names take, when one cannot hold them.

    At 170 mm, eleven features' ``p=3.7e-05`` labels ran over each other and the figure was
    refused (V.3). Feature i's label goes to row ``i % labels``, each row ``label_pt`` above
    the one before, and its name to row ``i % names``, each ``name_pt`` below; the page grows
    by the extra rows, so the axes keep their height. A label in ``left`` or ``right`` would
    cross that edge of the axes and is set flush inside it. ``Stack()``, one row of each,
    draws exactly what the figure drew before stacking existed.
    """

    labels: int = 1
    label_pt: float = 0.0
    names: int = 1
    name_pt: float = 0.0
    left: frozenset[int] = frozenset()
    right: frozenset[int] = frozenset()

    def extra(self) -> float:
        """The page height the extra rows add, in inches."""
        return ((self.labels - 1) * self.label_pt + (self.names - 1) * self.name_pt) / 72


def static(data: BoxData, theme: Theme) -> Any:
    """Render the matplotlib figure for one theme, its labels stacked where one row overlaps.

    Each pass measures the drawn labels and widens ``Stack`` only, so the passes end.

    Args:
        data: The contract.
        theme: One of ``style.THEMES``.

    Returns:
        A matplotlib Figure, not attached to pyplot; the caller saves it.
    """
    stack = Stack()
    while True:
        fig, ax, labels = _draw(data, theme, stack)
        wanted = stacking(fig, ax, labels, stack)
        if wanted == stack:
            return fig
        stack = wanted


def _draw(data: BoxData, theme: Theme, stack: Stack) -> tuple[Any, Any, list[Any]]:
    """The figure laid out with ``stack``'s rows, its axes and its p and r labels."""
    import matplotlib as mpl  # noqa: PLC0415 - heavy; imported after the runtime profile
    import seaborn as sns  # noqa: PLC0415 - heavy; imported after the runtime profile
    from matplotlib.figure import Figure  # noqa: PLC0415 - heavy; imported after the profile
    from matplotlib.transforms import offset_copy  # noqa: PLC0415 - heavy

    apply_matplotlib(theme)
    features = data.frame["feature"].unique(maintain_order=True).to_list()
    width, height = _inches()
    height += stack.extra()
    # Figure, not pyplot: a GUI backend rounds the canvas to whole pixels (84.84 mm, not 85).
    fig = Figure(figsize=(width, height))
    ax = fig.subplots()
    # matplotlib's own boxplot, not seaborn's: seaborn 0.13.2 passes bxp the vert argument
    # matplotlib 3.11 deprecates. Geometry as seaborn dodges two groups: boxes 0.4 wide at
    # +-0.2 around each feature, caps half a box.
    ink = {"color": theme.ink, "linewidth": LINE_WIDTHS["box"]}
    for offset, role, group in zip(OFFSETS, ("reference", "comparison"), data.groups,
                                   strict=True):  # fmt: skip
        values = [_values(data.frame, f, group) for f in features]
        positions = [i + offset for i in range(len(features))]
        ax.boxplot(
            values, positions=positions, widths=BOX_WIDTH, capwidths=BOX_WIDTH / 2,
            orientation="vertical",
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
    _pair_lines(ax, data, features, theme)
    ax.set_xticks(range(len(features)), features)
    ax.set_xlim(-0.5, len(features) - 0.5)
    top = float(data.frame["value"].to_numpy().max())
    labels = [_label(ax, data.caption.label(f), i, top, stack, len(features))
              for i, f in enumerate(features)]  # fmt: skip
    for i, tick in enumerate(ax.xaxis.get_major_ticks()):
        if i % stack.names:
            tick.set_pad(mpl.rcParams["xtick.major.pad"] + i % stack.names * stack.name_pt)
    low, high = ax.get_ylim()
    ax.set_ylim(top=high + (high - low) * HEADROOM)
    ax.set_xlabel(data.feature_label)
    ax.set_ylabel(fit_title(axis_title(data.value_label, data.unit), height))
    # Above the axes, and above every stacked row: inside, it collides with the p/r labels.
    rows_pt = (stack.labels - 1) * stack.label_pt
    over = {"bbox_transform": offset_copy(ax.transAxes, fig=fig, y=rows_pt, units="points")}
    ax.legend(title=None, loc="lower left", bbox_to_anchor=(0, 1), ncols=2, borderaxespad=0,
              **(over if rows_pt else {}))  # fmt: skip
    sns.despine(fig=fig)
    place_caption(fig, data.caption.text(), width, height)
    refit_y_title(fig, ax, axis_title(data.value_label, data.unit))
    return fig, ax, labels


def _label(ax: Any, text: str, i: int, top: float, stack: Stack, n: int) -> Any:
    """One feature's p and r label over its boxes, in its row and inside the axes."""
    row = i % stack.labels
    # Flush with an edge, one gap inside it, off the spine.
    x, ha, shift = ((-0.5, "left", STACK_GAP_PT) if i in stack.left
                    else (n - 0.5, "right", -STACK_GAP_PT) if i in stack.right
                    else (i, "center", 0.0))  # fmt: skip
    if not row and not shift:
        return ax.text(x, top, text, ha=ha, va="bottom", fontsize=SIZES["annotation"])
    return ax.annotate(text, (x, top), xytext=(shift, row * stack.label_pt),
                       textcoords="offset points", ha=ha, va="bottom",
                       fontsize=SIZES["annotation"])  # fmt: skip


def stacking(fig: Any, ax: Any, labels: list[Any], stack: Stack) -> Stack:
    """The rows the drawn labels and names need, never fewer than ``stack`` already has.

    Measured on the laid-out figure. While one row holds every label and name apart, the
    figure stays as it was drawn before stacking existed; once one does not, a label that
    would cross an edge of the axes is set flush with it, and every two in one row keep
    ``STACK_GAP_PT`` apart.

    Args:
        fig: The laid-out figure.
        ax: Its axes.
        labels: The p and r labels, in feature order.
        stack: The rows the figure was drawn with.

    Returns:
        ``stack`` itself when it holds every label and name apart, a wider one otherwise.
    """
    from matplotlib.backends.backend_agg import FigureCanvasAgg  # noqa: PLC0415 - heavy

    renderer = FigureCanvasAgg(fig).get_renderer()  # type: ignore[no-untyped-call]
    fig.draw(renderer)
    to_pt = 72 / float(fig.dpi)
    frame = ax.get_window_extent(renderer)
    centres = [float(ax.transData.transform((i, 0))[0]) for i in range(len(labels))]
    names = [t.get_window_extent(renderer) for t in ax.get_xticklabels()]
    boxes = [t.get_window_extent(renderer) for t in labels]
    spans = [(c - b.width / 2, c + b.width / 2) for c, b in zip(centres, boxes, strict=True)]
    named = [(b.x0, b.x1) for b in names]
    if stack == Stack() and _rows(spans, 0.0) == 1 and _rows(named, 0.0) == 1:
        return stack
    gap = STACK_GAP_PT / to_pt
    left = stack.left | {i for i, (a, _) in enumerate(spans) if a < frame.x0 + gap}
    right = stack.right | {i for i, (_, b) in enumerate(spans) if b > frame.x1 - gap}
    flush = [((frame.x0 + gap, frame.x0 + gap + w) if i in left
              else (frame.x1 - gap - w, frame.x1 - gap) if i in right else (a, b))
             for i, (a, b) in enumerate(spans) for w in [b - a]]  # fmt: skip
    return Stack(
        labels=max(stack.labels, _rows(flush, gap)),
        label_pt=max(stack.label_pt, max(b.height for b in boxes) * to_pt + STACK_GAP_PT),
        names=max(stack.names, _rows(named, gap)),
        name_pt=max(stack.name_pt, max(b.height for b in names) * to_pt + STACK_GAP_PT),
        left=frozenset(left), right=frozenset(right),
    )  # fmt: skip


def _rows(spans: list[tuple[float, float]], gap: float) -> int:
    """The fewest rows that keep every two spans of one row ``gap`` apart, span i in row i % k."""
    for rows in range(1, len(spans) + 1):
        if all(spans[k][0] >= spans[j][1] + gap or spans[j][0] >= spans[k][1] + gap
               for j in range(len(spans)) for k in range(j + rows, len(spans), rows)):  # fmt: skip
            return rows
    return max(len(spans), 1)


def _pair_lines(ax: Any, data: BoxData, features: list[str], theme: Theme) -> None:
    """One line per paired unit, from its reference point to its comparison point, under both.

    In the reference grey, which reads on both themes; the page draws it in the same grey.
    """
    from matplotlib.collections import LineCollection  # noqa: PLC0415 - heavy

    lines = paired_lines(data)
    segments = [[(i + OFFSETS[0], a), (i + OFFSETS[1], b)]
                for i, f in enumerate(features) for _, a, b in lines.get(f, [])]  # fmt: skip
    if segments:
        ax.add_collection(LineCollection(segments, colors=theme.reference,
                                         linewidths=LINE_WIDTHS["line"], zorder=2))  # fmt: skip


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
    features = data.frame["feature"].unique(maintain_order=True).to_list()
    at = {f: i for i, f in enumerate(features)}
    fig = go.Figure()
    if (pairs := _pair_trace(go, data, at)) is not None:
        fig.add_trace(pairs)
    for offset, role, group in zip(OFFSETS, ("reference", "comparison"), data.groups,
                                   strict=True):  # fmt: skip
        rows = data.frame.filter(pl.col("group") == group)
        extra = [
            [
                stats[f].n_reference if role == "reference" else stats[f].n_comparison,
                stats[f].p_adj,
                stats[f].r,
                " ".join(_lines(f)),
            ]
            for f in rows["feature"]
        ]
        custom = [[*h, *map(hover_value, e)]
                  for h, e in zip(hover_rows(rows, data.hover), extra, strict=True)]  # fmt: skip
        k = len(data.hover)
        hover_lines = [f"{c}: {hover_field(f'customdata[{i}]', schema[c])}"
                       for i, c in enumerate(data.hover)]  # fmt: skip
        # A numeric axis, each box at its feature's position plus its offset as the files draw
        # it, so a donor's line (``_pair_trace``) meets its two points.
        fig.add_trace(
            go.Box(
                x=[at[f] + offset for f in rows["feature"]],
                y=rows["value"].to_list(),
                width=BOX_WIDTH,
                name=group,
                marker={"color": ROLES[role], "opacity": POINT_ALPHA, "size": PLOTLY_MARKER_PX},
                boxpoints="all",
                jitter=0,
                pointpos=0,
                customdata=custom,
                hovertemplate="<br>".join(
                    [
                        f"{data.feature_label}: %{{customdata[{k + 3}]}}",
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
            "x": at[f],
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
        boxmode="overlay",
        xaxis={"title": {"text": data.feature_label}, "tickmode": "array",
               "tickvals": list(range(len(features))), "ticktext": [_tick(f) for f in features],
               "range": [-0.5, len(features) - 0.5], "zeroline": False, "showgrid": False},
        # Room above the tallest box for the p and r labels, as the static figure leaves.
        yaxis={"title": {"text": page_y_title(title, files_room(lambda t: static(data, t)))},
               "range": [low - 0.05 * span, high + HEADROOM * span]},
        annotations=labels,
    )  # fmt: skip
    y_range = (low - 0.05 * span, high + HEADROOM * span)
    wide, narrow = least_width(data, y_range), narrow_width(data, y_range)
    if sum(narrow[0]) < sum(wide[0]):
        fig.update_layout(meta={"least_width": narrow,
                                "narrow": _narrow(data, sum(wide[0]), len(labels))})  # fmt: skip
    else:
        # Turned, names of several lines can need more room than upright (Step-1's six F1
        # contrasts, 1,184 px against 621): the upright layout stays, scrolling below it.
        fig.update_layout(meta={"least_width": wide})
    return fig


def _tick(feature: str) -> str:
    """A feature's tick label on the page: each line of its name on its own line."""
    return "<br>".join(_lines(feature))


def _pair_trace(go: Any, data: BoxData, at: Mapping[str, int]) -> Any:
    """The paired units' lines as one trace, each from its reference point to its comparison's.

    Drawn first, so under the boxes' points; in the reference grey of both themes. None where
    nothing is paired.
    """
    lines = paired_lines(data)
    if not lines:
        return None
    x: list[float | None] = []
    y: list[float | None] = []
    units: list[str | None] = []
    for feature, pairs in lines.items():
        for unit, a, b in pairs:
            x += [at[feature] + OFFSETS[0], at[feature] + OFFSETS[1], None]
            y += [a, b, None]
            units += [unit, unit, None]
    return go.Scatter(x=x, y=y, mode="lines", showlegend=False, name=data.pairs,
                      line={"width": LINE_WIDTHS["line"], "color": THEMES["white"].reference},
                      customdata=units, connectgaps=False,
                      hovertemplate=f"{data.pairs}: %{{customdata}}<extra></extra>")  # fmt: skip


def _narrow(data: BoxData, below: float, n_labels: int) -> dict[str, Any]:
    """A phone's layout, swapped in below ``below`` px: each feature's p and r under its name.

    At 375 px the kidney page's four genes needed 480 px of upright labels and the box
    scrolled inside its card (C3.4). A feature's name and its p and r line go into one tick
    label turned a quarter, whose slot needs only two lines' height, and the labels over the
    boxes hide; plotly's automargin gives the turned labels their room below the axis.
    """
    features = data.frame["feature"].unique(maintain_order=True).to_list()
    ticks = [f"{_tick(f)}<br>{data.caption.label(f, 'label_line')}" for f in features]
    relayout: dict[str, Any] = {"xaxis.tickangle": -90, "xaxis.tickmode": "array",
                                "xaxis.tickvals": list(range(len(features))),
                                "xaxis.ticktext": ticks}  # fmt: skip
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
    lines = max(1 + len(_lines(f)) for f in features)
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
            *(text_width(line, tick) for line in _lines(f)),
            *(text_width(line, note) for line in data.caption.label(f).split("\n")),
        )
        for f in features
    )
    y_ticks = max(text_width(f"{v:.0f}", tick) for v in y_range)
    scaled = len(features) * (widest + BOX_GAP_PX) + y_ticks + on_page(SIZES["axis_label"])
    return [[scaled, float(PLOTLY_MARGIN_PX["l"] + PLOTLY_MARGIN_PX["r"])]]
