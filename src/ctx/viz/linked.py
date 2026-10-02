"""Linked brushing: one scatter of the samples beside a parallel-sets view of their labels.

On the page, a lasso on the scatter or a tap on a category of the parallel sets paints those
samples in the picked colour, in both views at once, and a live line counts what is painted
(``templates/_brush.js.j2``, whose logic ``tests/test_page_script.py`` runs under node). The
scatter sits beside the reading guide and the parallel sets under both, across the figure
card: two plotly figures, ``interactive()`` and ``below()`` (Rytis 2026-10-01). The figure
files draw the parallel sets as vector bands (``FlowData``); painted, the page draws both
views into its own PDF and SVG. ``static()``, ``interactive()`` and ``below()`` read the same
``LinkedData``.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import polars as pl

from ctx.viz.axes import axis_title
from ctx.viz.flow import FlowData, Node
from ctx.viz.metrics import text_width
from ctx.viz.scatter import FLOAT_HOVER
from ctx.viz.style import (
    DATA_COLOURS,
    MAX_CATEGORIES,
    PAGE_PX_PER_PT,
    PLOTLY_HEIGHT_PX,
    PLOTLY_MARGIN_PX,
    PLOTLY_MARKER_PX,
    PLOTLY_PARCATS,
    POINT_ALPHA,
    SIZES,
    THEMES,
    Theme,
    figure_inches,
    on_page,
    plotly_steps,
    plotly_templates,
)

# A page's PDF and SVG of a painted view draw the scatter over the parallel sets, as shares
# of the plot's height.
SCATTER_DOMAIN, PARCATS_DOMAIN = (0.58, 1.0), (0.0, 0.42)


@dataclass(frozen=True)
class LinkedData:
    """The contract both renderers read.

    Attributes:
        frame: One row per sample.
        columns: (unit, meaning) for every column of ``frame``.
        x: Column on the scatter's x axis.
        y: Column on the scatter's y axis.
        x_label: x-axis title, without the unit.
        y_label: y-axis title, without the unit.
        dimensions: Category columns of the parallel sets, left to right.
        dimension_labels: A title per dimension.
        colour_dimension: The dimension whose levels colour the static bands into them.
        colours: Level of ``colour_dimension`` to hex colour.
        caption: The caption, unwrapped.
        names: A name per paint colour, in ``DATA_COLOURS`` order, for the picker and summary.
        summary: The live line per paint colour, with ``{n}``, ``{colour}`` and ``{detail}``.
        empty: The line while nothing is painted.
        detail: The dimension the live line breaks each colour down by.
        hover: Columns shown on hover of a point.
        width: ``"single"`` or ``"double"`` column.
    """

    frame: pl.DataFrame
    columns: Mapping[str, tuple[str, str]]
    x: str
    y: str
    x_label: str
    y_label: str
    dimensions: tuple[str, ...]
    dimension_labels: tuple[str, ...]
    colour_dimension: str
    colours: Mapping[str, str]
    caption: str
    names: tuple[str, ...]
    summary: str
    empty: str
    detail: str
    hover: tuple[str, ...] = ()
    width: str = "double"

    def __post_init__(self) -> None:
        """Refuse a contract the renderers could only draw wrongly."""
        named = [self.x, self.y, *self.dimensions, *self.hover]
        if missing := [c for c in named if c not in self.frame.columns]:
            raise ValueError(f"missing columns {missing}")
        if len(self.dimension_labels) != len(self.dimensions):
            raise ValueError("one label per dimension")
        if {self.colour_dimension, self.detail} - set(self.dimensions):
            raise ValueError("the colour and detail dimensions must be dimensions")
        if len(self.colours) > MAX_CATEGORIES:
            raise ValueError(f"{len(self.colours)} colours; at most {MAX_CATEGORIES}")
        if stray := set(self.frame[self.colour_dimension].unique()) - set(self.colours):
            raise ValueError(f"{self.colour_dimension}: levels {sorted(stray)} have no colour")
        if len(self.names) != len(DATA_COLOURS):
            raise ValueError(f"one name per data colour: {len(DATA_COLOURS)}")
        if self.frame.select(self.dimensions).null_count().sum_horizontal().item():
            raise ValueError("a dimension has a missing value; name it as a level")

    def table(self) -> pl.DataFrame:
        """Every sample."""
        return self.frame

    def table_columns(self) -> Mapping[str, tuple[str, str]]:
        """(unit, meaning) for every column of ``table()``."""
        return self.columns

    def page_table(self) -> pl.DataFrame:
        """Every sample: the page script counts them by ``detail`` as it paints."""
        return self.frame

    def caption_text(self) -> str:
        """The caption, unwrapped."""
        return self.caption

    def flow(self) -> FlowData:
        """The parallel sets as a flow: a node per level of each dimension, samples as reads."""
        nodes, links = [], []
        for stage, column in enumerate(self.dimensions):
            for level in self.frame[column].unique(maintain_order=True).to_list():
                colour = self.colours.get(level) if column == self.colour_dimension else None
                nodes.append(Node(f"{stage}:{level}", str(level), stage, colour))
        for stage, (left, right) in enumerate(zip(self.dimensions, self.dimensions[1:],
                                                  strict=False)):  # fmt: skip
            pairs = self.frame.group_by(left, right, maintain_order=True).len()
            links += [(f"{stage}:{a}", f"{stage + 1}:{b}", float(n)) for a, b, n in pairs.rows()]
        return FlowData(
            nodes=tuple(nodes),
            links=pl.DataFrame(links, schema=["source", "target", "value"], orient="row"),
            columns={"source": ("none", "Dimension:level the samples leave."),
                     "target": ("none", "Dimension:level they reach."),
                     "value": ("samples", "Samples with both levels.")},
            stage_labels=self.dimension_labels, caption=self.caption, width=self.width,
        )  # fmt: skip

    def static(self, theme: Theme) -> Any:
        """The parallel sets, drawn as vector bands; see ``flow()``."""
        return self.flow().static(theme)

    def inches(self) -> tuple[float, float]:
        """The size of a page's own PDF and SVG of a painted view: both views and a caption."""
        return figure_inches(self.width, "linked")

    def controls(self) -> dict[str, Any]:
        """What the page script needs: the traces, the paint colours and the live line.

        The scatter is trace ``scatter`` of the figure beside the reading guide, the parallel
        sets trace ``parcats`` of the one under it.
        """
        return {
            "kind": "brush", "scatter": 0, "parcats": 0, "reference": THEMES["white"].reference,
            "colours": list(DATA_COLOURS), "names": list(self.names), "summary": self.summary,
            "empty": self.empty, "detail": self.detail,
            "above": list(SCATTER_DOMAIN), "below": list(PARCATS_DOMAIN),
        }  # fmt: skip

    def interactive(self) -> Any:
        """The scatter; the page script paints it and the parallel sets of ``below()``.

        Returns:
            A plotly Figure with the light template; the page switches themes.
        """
        import plotly.graph_objects as go  # noqa: PLC0415 - heavy; imported after the profile

        reference = THEMES["white"].reference
        n = self.frame.height
        x_title = axis_title(self.x_label, self.columns[self.x][0])
        y_title = axis_title(self.y_label, self.columns[self.y][0])
        template = "<br>".join([
            f"{x_title}: %{{x{FLOAT_HOVER}}}", f"{y_title}: %{{y{FLOAT_HOVER}}}",
            *(f"{c}: %{{customdata[{i}]}}" for i, c in enumerate(self.hover)),
        ]) + "<extra></extra>"  # fmt: skip
        fig = go.Figure(go.Scatter(
            x=self.frame[self.x].to_list(), y=self.frame[self.y].to_list(), mode="markers",
            showlegend=False, hovertemplate=template,
            customdata=self.frame.select(self.hover).rows() if self.hover else None,
            marker={"color": [reference] * n, "size": PLOTLY_MARKER_PX, "opacity": POINT_ALPHA},
        ))  # fmt: skip
        fig.update_layout(
            template=plotly_templates()["white"], height=PLOTLY_HEIGHT_PX["linked"],
            dragmode="lasso",
            xaxis={"title": {"text": x_title}},
            yaxis={"title": {"text": y_title}},
        )  # fmt: skip
        return fig

    def below(self) -> Any:
        """The parallel sets, one sample a line, grey until painted; drawn across the card.

        Returns:
            A plotly Figure with the light template; the page switches themes.
        """
        import plotly.graph_objects as go  # noqa: PLC0415 - heavy; imported after the profile

        palette = [THEMES["white"].reference, *DATA_COLOURS]
        fig = go.Figure(go.Parcats(
            arrangement="freeform", hoveron="color",
            dimensions=[{"label": label, "values": self.frame[c].cast(pl.String).to_list()}
                        for c, label in zip(self.dimensions, self.dimension_labels, strict=True)],
            line={"color": [0] * self.frame.height, "cmin": -0.5, "cmax": len(palette) - 0.5,
                  "colorscale": plotly_steps(palette), "shape": "hspline"},
        ))  # fmt: skip

        margin = self.sides(on_page(SIZES["tick"]), PLOTLY_MARGIN_PX)
        fig.update_layout(template=plotly_templates()["white"], height=PLOTLY_HEIGHT_PX["parcats"],
                          margin=margin, meta={"least_width": self._least(margin)})  # fmt: skip
        return fig

    def _widest(self, column: str, size: float) -> float:
        """The widest level of a dimension, set at ``size``."""
        return max(text_width(str(v), size) for v in self.frame[column].unique())

    def sides(self, size: float, margin: Mapping[str, float]) -> dict[str, float]:
        """Left and right margins that hold the outer dimensions' levels, set at ``size``.

        plotly sets the first dimension's levels left of its bars and the last's right of
        them, outside the plot, where a margin too narrow cut them off.
        """
        return {"l": self._widest(self.dimensions[0], size) + margin["l"],
                "r": self._widest(self.dimensions[-1], size) + margin["r"]}  # fmt: skip

    def _least(self, margin: Mapping[str, float]) -> list[list[float]]:
        """The narrowest plot that keeps its titles and levels a gap apart, at any zoom.

        The gap, ``PLOTLY_PARCATS`` gap_px, lies between neighbouring dimension titles and
        between each level and the bar before it. plotly spaces the bars evenly, so every step
        takes the widest one needs. Each title is centred on its bar, and a middle dimension's
        levels end left of its bar (the last one's start right of it, in the margin). At 375 px
        a phone squeezed four titles into 191 px and they ran into each other; the page draws
        the parallel sets no narrower than this and scrolls them sideways. Text, the gap and
        the margins grow with a page's zoom, plotly's bars, insets and label offsets do not:
        opened alone at twice the type, the demo's sets need 1,392 px, not twice 786 (Rytis
        2026-10-02, the width capped), where Chromium measured 32.2 px from the nearest level to
        the bar before it (2026-10-02).

        Returns:
            [scaled, fixed] px pairs, one per title step and per middle dimension's levels;
            at zoom s the plot needs the largest s * scaled + fixed.
        """
        title, tick, keep = on_page(SIZES["axis_label"]), on_page(SIZES["tick"]), PLOTLY_PARCATS
        titles = [text_width(label, title) for label in self.dimension_labels]
        needs = [
            ((left + right) / 2, 0.0) for left, right in zip(titles[:-1], titles[1:], strict=True)
        ]
        needs += [(self._widest(c, tick), keep["label_px"] + keep["bar_px"])
                  for c in self.dimensions[1:-1]]  # fmt: skip
        steps, sides = len(titles) - 1, margin["l"] + margin["r"]
        bars = 2 * keep["inset_px"] + keep["bar_px"]
        return [[sides + steps * (text + keep["gap_px"]), bars + steps * plotly]
                for text, plotly in needs]  # fmt: skip

    def export_layers(self) -> dict[str, Any]:
        """The parallel sets' outer margins at the figure files' sizes, in pt.

        A painted view's PDF and SVG draw the scatter's layout; without these, levels longer
        than its margin ran off the page (review, 2026-10-01).
        """
        pt = {side: px / PAGE_PX_PER_PT for side, px in PLOTLY_MARGIN_PX.items()}
        return {"sides": self.sides(SIZES["tick"], pt)}
