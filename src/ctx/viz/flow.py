"""Flow diagram: where every read goes, stage by stage, as a Sankey.

Nodes are the fates at each stage, stacked in columns, their height the reads they hold;
links are bands whose width is the reads that move. ``static()`` draws it with vector paths,
since matplotlib has no Sankey of this kind; ``interactive()`` uses plotly's. Both read the
same ``FlowData``, which refuses a node whose reads in and out differ: every read has a fate.
"""

import textwrap
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import polars as pl

from ctx.viz.annotate import place_caption
from ctx.viz.style import (
    AVERAGE_CHAR_EM,
    FLOW,
    LINE_SPACING,
    LINE_WIDTHS,
    PLOTLY_FLOW,
    PLOTLY_HEIGHT_PX,
    POINT_ALPHA,
    SIZES,
    THEMES,
    Theme,
    apply_matplotlib,
    figure_inches,
    plotly_templates,
    tint,
)

LINK_COLUMNS = ("source", "target", "value")


@dataclass(frozen=True)
class Node:
    """One fate at one stage; its colour colours the links into it (None: reference grey)."""

    id: str
    label: str
    stage: int
    colour: str | None = None


@dataclass(frozen=True)
class FlowData:
    """The contract both renderers read.

    Attributes:
        nodes: Every node, in top-to-bottom order within its stage.
        links: One row per link: ``source`` and ``target`` node ids and a ``value``.
        columns: (unit, meaning) for every column of ``links``.
        stage_labels: A title above each stage's column, in stage order.
        caption: The caption, unwrapped.
        width: ``"single"`` or ``"double"`` column.
    """

    nodes: tuple[Node, ...]
    links: pl.DataFrame
    columns: Mapping[str, tuple[str, str]]
    stage_labels: tuple[str, ...]
    caption: str
    width: str = "double"

    def __post_init__(self) -> None:
        """Refuse a flow that loses or invents reads, or names a node it does not have."""
        if missing := [c for c in LINK_COLUMNS if c not in self.links.columns]:
            raise ValueError(f"missing columns {missing}")
        ids = {n.id for n in self.nodes}
        if unknown := {*self.links["source"], *self.links["target"]} - ids:
            raise ValueError(f"links name no node: {sorted(unknown)}")
        if self.links.select("source", "target").is_duplicated().any():
            raise ValueError("a source and target are linked twice")
        if (self.links["value"] < 0).any():
            raise ValueError("a link carries a negative value")
        stage = {n.id: n.stage for n in self.nodes}
        if bad := [(s, t) for s, t in self.links.select("source", "target").rows()
                   if stage[t] <= stage[s]]:  # fmt: skip
            raise ValueError(f"links must run to a later stage: {bad[:5]}")
        if len(self.stage_labels) != len({n.stage for n in self.nodes}):
            raise ValueError("one stage label per stage")
        into, out = self._totals("target"), self._totals("source")
        if unbalanced := [i for i in ids if i in into and i in out
                          and abs(into[i] - out[i]) > 1e-9 * max(into[i], 1)]:  # fmt: skip
            raise ValueError(f"reads in and out differ at {sorted(unbalanced)}")

    def _totals(self, side: str) -> dict[str, float]:
        """Reads into (``target``) or out of (``source``) each node."""
        return dict(self.links.group_by(side).agg(pl.col("value").sum()).rows())

    def held(self) -> dict[str, float]:
        """Reads each node holds: what comes in, or what leaves a first-stage node."""
        into, out = self._totals("target"), self._totals("source")
        return {n.id: into.get(n.id, out.get(n.id, 0.0)) for n in self.nodes}

    def table(self) -> pl.DataFrame:
        """Every link, with its nodes' labels."""
        label = {n.id: n.label for n in self.nodes}
        return self.links.with_columns(
            source_label=pl.col("source").replace_strict(label),
            target_label=pl.col("target").replace_strict(label),
        )

    def table_columns(self) -> Mapping[str, tuple[str, str]]:
        """(unit, meaning) for every column of ``table()``."""
        return {**self.columns, "source_label": ("none", "What the source node is."),
                "target_label": ("none", "What the target node is.")}  # fmt: skip

    def page_table(self) -> pl.DataFrame:
        """Every link: the page draws them all."""
        return self.table()

    def caption_text(self) -> str:
        """The caption, unwrapped."""
        return self.caption

    def _layout(self, gap: float) -> tuple[dict[str, tuple[float, float, float]], float]:
        """(x, top, height) of every node in axes fractions, and the height of one read.

        One scale for every stage, set by the fullest, so equal reads draw equally tall.
        """
        stages = sorted({n.stage for n in self.nodes})
        room = 1.0 - FLOW["label_room"]
        x = {s: room * i / max(len(stages) - 1, 1) for i, s in enumerate(stages)}
        held = self.held()
        column = {s: [n for n in self.nodes if n.stage == s] for s in stages}
        top = 1.0 - 2 * gap
        # Half a label below the lowest node: its label is centred on it, whatever its height.
        scale = min((top - gap * (len(ns) - 1) - gap / 2) / max(sum(held[n.id] for n in ns), 1e-12)
                    for ns in column.values())  # fmt: skip
        placed = {}
        for s, ns in column.items():
            y = top
            for n in ns:
                placed[n.id] = (x[s], y, held[n.id] * scale)
                y -= held[n.id] * scale + gap
        return placed, scale

    def _bands(self, placed: Mapping[str, tuple[float, float, float]],
               scale: float) -> list[tuple[float, float, float, float, float, str]]:  # fmt: skip
        """(x0, top0, x1, top1, height, target) of every link, stacked at both ends.

        A node's incoming links stack from its top in their sources' order, its outgoing
        links in their targets' order, so bands do not cross at a node.
        """
        order = {n.id: i for i, n in enumerate(self.nodes)}
        rows = self.links.select(LINK_COLUMNS).rows()
        into: dict[tuple[str, str], float] = {}
        used: dict[str, float] = {}
        for source, target, value in sorted(rows, key=lambda r: (order[r[1]], order[r[0]])):
            into[(source, target)] = used.get(target, 0.0)
            used[target] = into[(source, target)] + value * scale
        bands = []
        used = {}
        for source, target, value in sorted(rows, key=lambda r: (order[r[0]], order[r[1]])):
            x0, top0, _ = placed[source]
            x1, top1, _ = placed[target]
            offset = used.get(source, 0.0)
            used[source] = offset + value * scale
            bands.append((x0 + FLOW["node"], top0 - offset, x1, top1 - into[(source, target)],
                          value * scale, target))  # fmt: skip
        return bands

    def static(self, theme: Theme) -> Any:
        """Render the matplotlib figure for one theme: vector nodes and bands.

        Args:
            theme: One of ``style.THEMES``.

        Returns:
            A matplotlib Figure, not attached to pyplot.
        """
        from matplotlib.figure import Figure  # noqa: PLC0415 - heavy

        apply_matplotlib(theme)
        width, height = figure_inches(self.width, "flow")
        line = SIZES["annotation"] * LINE_SPACING / 72
        most = max(sum(n.stage == s for n in self.nodes) for s in {n.stage for n in self.nodes})
        # Long level names wrap into their column; when the wrapped labels need more height than
        # the figure has, it grows, as a heatmap does for its labels. Two passes suffice: the
        # wrap depends on the width alone.
        for _ in range(2):
            fig = Figure(figsize=(width, height))
            ax = fig.subplots()
            ax.set_axis_off()
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            place_caption(fig, self.caption, width, height)
            axes_inches = ax.get_position().height * height
            placed, scale = self._layout(line / axes_inches)
            labels = self._labels(placed, ax.get_position().width * width)
            lines = max(label.count("\n") + 1 for label in labels.values())
            # _layout leaves gap * (nodes + 1.5) of the axes to labels and titles.
            needed = line * lines * (most + 1.5) / FLOW["label_share"]
            if needed <= axes_inches:
                break
            height += (needed - axes_inches) * 1.1
        placed, scale = self._layout(line * lines / axes_inches)
        colour = {n.id: n.colour or theme.reference for n in self.nodes}
        for band in self._bands(placed, scale):
            ax.add_patch(_band_patch(band[:5], colour[band[5]]))
        self._draw_nodes(ax, placed, colour, theme, labels)
        return fig

    def _labels(self, placed: Mapping[str, tuple[float, float, float]],
                axes_inches: float) -> dict[str, str]:  # fmt: skip
        """Each node's label and reads, wrapped to the room before the next stage's bars.

        Long level names ("Reads \u2264 76 nt (min 35 nt)") ran into the next stage on one line.
        """
        held = self.held()
        starts = sorted({x for x, _, _ in placed.values()})
        room = {x: (nxt if nxt else 1.0) - x - FLOW["node"] - 2 * FLOW["pad"]
                for x, nxt in zip(starts, [*starts[1:], None], strict=True)}  # fmt: skip
        chars = {x: max(1, int(r * axes_inches * 72 / (SIZES["annotation"] * AVERAGE_CHAR_EM)))
                 for x, r in room.items()}  # fmt: skip
        return {n.id: textwrap.fill(f"{n.label} {held[n.id]:,.0f}", chars[placed[n.id][0]],
                                    break_long_words=False)
                for n in self.nodes}  # fmt: skip

    def _draw_nodes(self, ax: Any, placed: Mapping[str, tuple[float, float, float]],
                    colour: Mapping[str, str], theme: Theme,
                    labels: Mapping[str, str]) -> None:  # fmt: skip
        """Each node's bar and its label with its reads; each stage's title above it."""
        from matplotlib.patches import Rectangle  # noqa: PLC0415 - heavy

        for n in self.nodes:
            x, top, h = placed[n.id]
            ax.add_patch(Rectangle((x, top - h), FLOW["node"], h, facecolor=colour[n.id],
                                   linewidth=LINE_WIDTHS["band_edge"]))  # fmt: skip
            ax.text(x + FLOW["node"] + FLOW["pad"], top - h / 2, labels[n.id], va="center",
                    ha="left", fontsize=SIZES["annotation"], color=theme.ink,
                    linespacing=LINE_SPACING)  # fmt: skip
        stages = sorted({n.stage for n in self.nodes})
        for s, label in zip(stages, self.stage_labels, strict=True):
            x = min(placed[n.id][0] for n in self.nodes if n.stage == s)
            ax.text(x, 1.0, label, va="top", ha="left", fontsize=SIZES["axis_label"],
                    color=theme.ink)  # fmt: skip

    def interactive(self) -> Any:
        """Render plotly's Sankey; hover gives every link's and node's reads.

        Returns:
            A plotly Figure with the light template; the page switches themes.
        """
        import plotly.graph_objects as go  # noqa: PLC0415 - heavy; imported after the profile

        light = THEMES["white"]
        index = {n.id: i for i, n in enumerate(self.nodes)}
        colour = {n.id: n.colour or light.reference for n in self.nodes}
        rows = self.links.select(LINK_COLUMNS).rows()
        fig = go.Figure(go.Sankey(
            arrangement="snap", valueformat=",.0f",
            node={"label": [n.label for n in self.nodes], "pad": PLOTLY_FLOW["pad_px"],
                  "thickness": PLOTLY_FLOW["node_px"], "color": [colour[n.id] for n in self.nodes],
                  "line": {"width": LINE_WIDTHS["band_edge"]}},
            link={"source": [index[s] for s, _, _ in rows],
                  "target": [index[t] for _, t, _ in rows],
                  "value": [v for *_, v in rows],
                  "color": [tint(colour[t], POINT_ALPHA) for _, t, _ in rows]},
        ))  # fmt: skip
        fig.update_layout(template=plotly_templates()["white"], height=PLOTLY_HEIGHT_PX["flow"])
        return fig


def _band_patch(band: tuple[float, ...], colour: str) -> Any:
    """A link as a filled band: two cubic curves between the nodes, joined at each end."""
    from matplotlib.patches import PathPatch  # noqa: PLC0415 - heavy
    from matplotlib.path import Path  # noqa: PLC0415 - heavy

    x0, top0, x1, top1, h = band
    mid = (x0 + x1) / 2
    vertices = [(x0, top0), (mid, top0), (mid, top1), (x1, top1), (x1, top1 - h),
                (mid, top1 - h), (mid, top0 - h), (x0, top0 - h), (x0, top0)]  # fmt: skip
    codes = [Path.MOVETO, *[Path.CURVE4] * 3, Path.LINETO, *[Path.CURVE4] * 3, Path.CLOSEPOLY]
    return PathPatch(Path(vertices, codes), facecolor=colour, alpha=POINT_ALPHA,
                     linewidth=LINE_WIDTHS["band_edge"])  # fmt: skip
