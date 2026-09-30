"""Linked brushing: one scatter of the samples beside a parallel-sets view of their labels.

On the page, a lasso on the scatter or a tap on a category of the parallel sets paints those
samples in the picked colour, in both views at once, and a live line counts what is painted
(``templates/_brush.js.j2``, whose logic ``tests/test_page_script.py`` runs under node). The
PDF and SVG draw the parallel sets as vector bands (``FlowData``); the scatter is its own
figure. ``static()`` and ``interactive()`` read the same ``LinkedData``.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import polars as pl

from ctx.viz.axes import axis_title
from ctx.viz.flow import FlowData, Node
from ctx.viz.scatter import FLOAT_HOVER
from ctx.viz.style import (
    DATA_COLOURS,
    MAX_CATEGORIES,
    PLOTLY_HEIGHT_PX,
    PLOTLY_MARKER_PX,
    POINT_ALPHA,
    THEMES,
    Theme,
    plotly_steps,
    plotly_templates,
)

# The scatter over the top of the plot, the parallel sets under it, as shares of its height.
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

    def controls(self) -> dict[str, Any]:
        """What the page script needs: the traces, the paint colours and the live line."""
        return {
            "kind": "brush", "scatter": 0, "parcats": 1, "reference": THEMES["white"].reference,
            "colours": list(DATA_COLOURS), "names": list(self.names), "summary": self.summary,
            "empty": self.empty, "detail": self.detail,
        }  # fmt: skip

    def interactive(self) -> Any:
        """The scatter over the parallel sets; the page script paints both.

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
        palette = [reference, *DATA_COLOURS]
        fig.add_trace(go.Parcats(
            domain={"y": list(PARCATS_DOMAIN)}, arrangement="freeform", hoveron="color",
            dimensions=[{"label": label, "values": self.frame[c].cast(pl.String).to_list()}
                        for c, label in zip(self.dimensions, self.dimension_labels, strict=True)],
            line={"color": [0] * n, "cmin": -0.5, "cmax": len(palette) - 0.5,
                  "colorscale": plotly_steps(palette), "shape": "hspline"},
        ))  # fmt: skip
        fig.update_layout(
            template=plotly_templates()["white"], height=PLOTLY_HEIGHT_PX["linked"],
            dragmode="lasso",
            xaxis={"title": {"text": x_title}},
            yaxis={"title": {"text": y_title}, "domain": list(SCATTER_DOMAIN)},
        )  # fmt: skip
        return fig
