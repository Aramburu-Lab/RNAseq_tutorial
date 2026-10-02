"""One figure as a page prints it: its plotly div, its controls, its files and its rows.

``prepare`` gives the figure macro (``templates/_figure.html.j2``) everything it prints for a
``FigureBlock``; ``builder`` renders the pages around it. Opened alone, a figure is drawn at
``style.ALONE_SCALE`` (Rytis 2026-10-01), no wider than its box (Rytis 2026-10-02).
"""

import json
import math
from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

from markupsafe import Markup

from ctx.io.naming import figure_files
from ctx.report import export, text
from ctx.viz.style import PAGE_PLOT_PX, PHONE_PLOT_PX, THEMES, zoom

if TYPE_CHECKING:
    from ctx.report.builder import FigureBlock

PLOTLY_CONFIG = {
    "displaylogo": False,
    "responsive": True,
    # No toolbar: over a phone-width plot it covers the legend, and drag-to-zoom and
    # double-click-to-reset work without it. A PNG export would carry no statistics anyway:
    # the caption sits beside the plot, not in it.
    "displayModeBar": False,
}


def json_for_script(value: Any) -> Markup:
    """JSON safe inside a ``<script>`` element: every ``<`` is escaped as a unicode escape.

    ``allow_nan=False``: a browser's JSON.parse rejects NaN and Infinity, which Python writes
    by default; rows go through ``finite`` first.
    """
    return Markup(json.dumps(value, default=str, allow_nan=False).replace("<", "\\u003c"))


def finite(row: Mapping[str, Any]) -> dict[str, Any]:
    """A row with every NaN or infinite number as null, which is what they mean on a page."""
    return {k: None if isinstance(v, float) and not math.isfinite(v) else v for k, v in row.items()}


def figure_html(figure: Any, div_id: str, scale: float = 1.0, part: str = "beside") -> Markup:
    """The plotly div and its drawing script, without plotly.js itself.

    At a ``scale`` other than 1, opened alone, every length of the figure is ``scale`` times its
    page size, and its width ``scale`` times the width the page gives a plot there (``part``)
    or the width of its box, whichever is less (Rytis 2026-10-02): plotly redraws it at its
    box's width as the window changes. It is never narrower than ``scale`` times a phone's
    plot, the narrowest a page lays a figure out for, nor than a ``layout.meta.least_width``
    of [scaled, fixed] px pairs names: the largest ``scale`` * scaled + fixed. Narrower than
    that, its box scrolls sideways under the hint.
    """
    import plotly.io as pio  # noqa: PLC0415 - heavy; imported after the runtime profile

    figure = figure.to_plotly_json()
    width, least = "100%", 0.0
    if scale != 1:
        figure = zoom(figure, scale)
        width = f"min(100%, {PAGE_PLOT_PX[part] * scale:.0f}px)"
        # At twice the type a phone's plot held none of the band map's row labels (2026-10-02).
        least = PHONE_PLOT_PX * scale
    # The parallel sets' titles ran into each other at 375 px; they name their least width.
    needs = (figure["layout"].get("meta") or {}).get("least_width")
    if needs:
        least = max(least, *(scale * a + b for a, b in needs))
    if least:
        width = f"max({width}, {least:.0f}px)"
    html = pio.to_html(
        figure, include_plotlyjs=False, full_html=False, div_id=div_id, config=PLOTLY_CONFIG,
        default_width=width,
    )  # fmt: skip
    return Markup(html)


def files(block: "FigureBlock") -> list[tuple[str, str, str]]:
    """(label, theme or "", file name) of every file a reader may open beside the figure."""
    data = text.word("ui", "figures.data")
    if not block.stem:
        return [(data, "", block.data_file)]
    names = figure_files(block.stem)
    rendered = [
        (ext.upper(), theme, name)
        for name in names
        for theme in THEMES
        for ext in ("pdf", "svg")
        if name == f"{block.stem}_{theme}.{ext}"
    ]
    return [
        (data, "", block.data_file),
        (text.word("ui", "figures.columns"), "", f"{block.stem}_README.md"),
        *rendered,
    ]


def prepare(
    block: "FigureBlock", index: int, number: int, level: int, scale: float = 1.0
) -> dict[str, Any]:
    """Everything the figure macro prints for one figure; opened alone, at ``scale``."""
    controls = export.zoomed(block.controls, scale) if block.controls else None
    exporting = bool(controls and controls.get("export"))
    # Every character the view's PDF can set, which the page's embedded fonts must keep.
    drawn = [block.figure, *([block.below] if block.below is not None else [])]
    carried = export.characters(*(f.to_dict() for f in drawn), controls) if exporting else None
    below = block.below
    if below is not None:
        below = figure_html(below, f"ctx-below-{index}", scale, "across")
    return {
        "html": figure_html(block.figure, f"ctx-figure-{index}", scale),
        "below": below,
        "caption": block.caption,
        "rows": json_for_script([finite(r) for r in block.rows]),
        "data_file": block.data_file,
        "index": index,
        "number": number,
        "level": level,
        "text": block.text,
        "title": block.text.title,
        "tab": block.text.tab or block.text.title or text.fill("ui", "figure_n", n=number),
        "page": f"{block.stem}.html" if block.stem else "",
        "files": files(block),
        "controls": json_for_script(controls) if controls else None,
        "kind": block.controls["kind"] if block.controls else None,
        "export": exporting,
        "characters": carried,
        "defaults": defaults(block.controls) if block.controls else None,
    }


def defaults(controls: Mapping[str, Any]) -> dict[str, Any]:
    """What the controls show before the reader moves them: the default rule, no paint, or all.

    A filter's live line is ``ui.toml`` [filter] summary, filled here with the numbers of every
    bar and passed unfilled for the page script to fill again (``ctxFilter.summary``).
    """
    if controls["kind"] == "brush":
        return {"colours": list(zip(controls["colours"], controls["names"], strict=True)),
                "empty": controls["empty"]}  # fmt: skip
    if controls["kind"] == "filter":
        return {"filters": list(controls["filters"].items()),
                "summary": text.fill("ui", "filter.summary", **controls["start"]),
                "template": text.raw("ui", "filter.summary"),
                "none": text.word("ui", "filter.none")}  # fmt: skip
    return {
        "fc": f"{controls['lfc_min']:g}",
        "bases": list(controls["bases"].items()),
        "alpha": f"{controls['alpha']:g}",
        "presets": sorted({*controls["presets"], f"{controls['min_cutoff']:g}"}, key=float),
        "summary": controls["default_summary"],
    }
