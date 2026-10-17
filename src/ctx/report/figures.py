"""One figure as a page prints it: its plotly div, its controls, its files and its rows.

``prepare`` gives the figure macro (``templates/_figure.html.j2``) everything it prints for a
``FigureBlock``; ``builder`` renders the pages around it. Opened alone, a figure is drawn at
``style.ALONE_SCALE`` (Rytis 2026-10-01), no wider than its box (Rytis 2026-10-02), and at the
module page's size where its box is narrower than that scale of a phone's plot (Rytis
2026-10-07).
"""

import json
import math
import string
from collections.abc import Mapping, Sequence
from dataclasses import replace
from typing import TYPE_CHECKING, Any

from markupsafe import Markup, escape

from ctx.io.naming import figure_files
from ctx.report import export, marks, text
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


def columnar(rows: Sequence[Mapping[str, Any]]) -> dict[str, list[Any]]:
    """A figure's page rows as column names once and one list of values per row.

    Rows as objects repeat every column name on every row: on M5's breast-cohort volcano the
    names were 39 % of its 4,985,574 bytes and put the page over its cap. ``ctxRows`` in the
    page script turns this back into one object per row, so every reader of the rows is
    unchanged.
    """
    columns = list(rows[0]) if rows else []
    return {"columns": columns, "data": [list(finite(r).values()) for r in rows]}


def _width(scale: float, part: str, needs: list[list[float]]) -> str:
    """The CSS width of a plot's box at ``scale``: its card's, capped opened alone.

    Opened alone, ``scale`` times the width a module page gives a plot there (``part``) or the
    box's, whichever is less (Rytis 2026-10-02); never narrower than a ``layout.meta.least_width``
    of [scaled, fixed] px pairs names, the largest ``scale`` * scaled + fixed: the parallel
    sets' titles ran into each other at 375 px. Narrower, the box scrolls under the hint.
    """
    width = "100%" if scale == 1 else f"min(100%, {PAGE_PLOT_PX[part] * scale:.0f}px)"
    least = max((scale * a + b for a, b in needs), default=0.0)
    return f"max({width}, {least:.0f}px)" if least else width


def _plain(node: Any) -> Any:
    """A plotly figure's JSON as the page receives it: arrays as lists, NaN as null."""
    from plotly.utils import PlotlyJSONEncoder  # noqa: PLC0415 - heavy

    return json.loads(json.dumps(node, cls=PlotlyJSONEncoder))


def _differ(page: Any, zoomed: Any, path: str, out: dict[str, Any]) -> None:
    """Every value ``zoom`` changed, at its plotly attribute path, as the page holds it.

    A list of objects (annotations, shapes) is walked item by item, ``annotations[0].font``; any
    other list is one value, as a marker's sizes are.
    """
    if isinstance(page, dict) and isinstance(zoomed, dict):
        for key, value in page.items():
            _differ(value, zoomed.get(key), f"{path}.{key}" if path else key, out)
    elif (isinstance(page, list) and isinstance(zoomed, list) and len(page) == len(zoomed)
          and page and all(isinstance(v, dict) for v in page)):  # fmt: skip
        for i, (a, b) in enumerate(zip(page, zoomed, strict=True)):
            _differ(a, b, f"{path}[{i}]", out)
    elif page != zoomed:
        out[path] = page


def page_layout(page: Mapping[str, Any], zoomed: Mapping[str, Any], scale: float, part: str
                ) -> dict[str, Any]:  # fmt: skip
    """What turns a figure opened alone back into the module page's, below the floor.

    Rytis, 2026-10-07: below the phone width a figure opened alone is drawn as the page draws
    it. Wherever its box is narrower than ``scale`` times a phone's plot, where the doubled
    figure was held to that width and scrolled, the page script swaps this in through the
    narrow mechanism (``_page.js.j2``, fit): every value ``zoom`` changed, taken back from the
    page's own JSON, so no second list of lengths exists; the box's width as the module page
    sets it; and a heatmap's own narrow layout at the page's size. The template and ``meta``
    are not in it: the script swaps the page's templates in, and ``meta.zoom`` tells the export
    the page's zoom.

    Args:
        page: The figure's JSON at the module page's size.
        zoomed: The same at ``scale``.
        scale: The zoom opened alone.
        part: Where the plot sits, as ``PAGE_PLOT_PX`` names it.

    Returns:
        ``below``, the floor in px; ``width``, the box's CSS width; ``relayout``, attribute
        path to the page's value; ``restyle``, [values by path, traces] pairs, one per set of
        traces that share their paths; and ``narrow``, the page's own, where it has one.
    """
    page, zoomed = _plain(page), _plain(zoomed)
    skip = ("template", "meta")
    relayout: dict[str, Any] = {}
    _differ({k: v for k, v in page["layout"].items() if k not in skip}, zoomed["layout"], "",
            relayout)  # fmt: skip
    traces: dict[tuple[str, ...], list[tuple[int, dict[str, Any]]]] = {}
    for i, (a, b) in enumerate(zip(page["data"], zoomed["data"], strict=True)):
        changed: dict[str, Any] = {}
        _differ(a, b, "", changed)
        if changed:
            traces.setdefault(tuple(changed), []).append((i, changed))
    restyle = [[{k: [c[k] for _, c in group] for k in keys}, [i for i, _ in group]]
               for keys, group in traces.items()]  # fmt: skip
    meta = page["layout"].get("meta") or {}
    out = {"below": PHONE_PLOT_PX * scale, "width": _width(1.0, part, meta.get("least_width", [])),
           "relayout": relayout, "restyle": restyle}  # fmt: skip
    return out | ({"narrow": meta["narrow"]} if "narrow" in meta else {})


def figure_html(figure: Any, div_id: str, scale: float = 1.0, part: str = "beside") -> Markup:
    """The plotly div and its drawing script, without plotly.js itself.

    At a ``scale`` other than 1, opened alone, every length of the figure is ``scale`` times its
    page size, in a box ``_width`` caps, and plotly redraws it at its box's width as the window
    changes. Its layout's ``meta.zoom`` is ``scale``, and ``meta.page`` the page's own layout,
    which the page script swaps in where the box is narrower than ``scale`` times a phone's
    plot (``page_layout``, Rytis 2026-10-07).
    """
    import plotly.io as pio  # noqa: PLC0415 - heavy; imported after the runtime profile

    figure = figure.to_plotly_json()
    needs = (figure["layout"].get("meta") or {}).get("least_width") or []
    if scale != 1:
        zoomed = zoom(figure, scale)
        meta = zoomed["layout"].get("meta") or {}
        zoomed["layout"]["meta"] = {**meta, "zoom": scale,
                                    "page": page_layout(figure, zoomed, scale, part)}  # fmt: skip
        figure = zoomed
    html = pio.to_html(
        figure, include_plotlyjs=False, full_html=False, div_id=div_id, config=PLOTLY_CONFIG,
        default_width=_width(scale, part, needs),
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


def lead_html(lead: Mapping[str, Any]) -> Markup:
    """A lead that follows its filter, each recomputed number in a span the page refills.

    The span is named for its placeholder (``data-fill``); the words around it print as any
    lead does, a backticked span as code and a glossary mark as its display text.
    """
    values = {**lead["fixed"], **lead["start"]}
    out = []
    for literal, name, _, _ in string.Formatter().parse(lead["template"]):
        out.append(marks.html(literal))
        if name is None:
            continue
        shown = escape(values[name])
        out.append(Markup(f'<span data-fill="{name}">{shown}</span>')
                   if name in lead["start"] else shown)  # fmt: skip
    return Markup("").join(out)


def prepare(
    block: "FigureBlock", index: int, number: int, level: int, scale: float = 1.0
) -> dict[str, Any]:
    """Everything the figure macro prints for one figure; opened alone, at ``scale``.

    Opened alone a ``view`` select is dropped: its options are drawn at the module page's
    size, which a figure zoomed alone cannot take (D-051: a coupling a standalone page cannot
    honour is dropped).
    """
    if scale != 1 and block.controls and block.controls.get("kind") == "view":
        block = replace(block, controls=None)
    controls = dict(block.controls) if block.controls else None
    exporting = bool(controls and controls.get("export"))
    # Every character the view's PDF can set, which the page's embedded fonts must keep.
    drawn = [block.figure, *([block.below] if block.below is not None else [])]
    carried = export.characters(*(f.to_dict() for f in drawn), controls) if exporting else None
    below = block.below
    if below is not None:
        below = figure_html(below, f"ctx-below-{index}", scale, "across")
    following = controls.get("lead") if controls else None
    return {
        "html": figure_html(block.figure, f"ctx-figure-{index}", scale),
        "below": below,
        "lead": lead_html(following) if following else marks.html(block.text.lead),
        "caption": block.caption,
        "rows": json_for_script(columnar(block.rows)),
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
        "select": controls.get("select") if controls else None,
    }


# What a form prints for the kinds that only pass on parts of their settings: a view's
# option names (D-051), M6's live ORA's words and backgrounds (GS.5), M9's checkbox sets and
# its list from a link (E.1a).
_SHAPED = {
    "view": lambda c: {"options": [o["label"] for o in c["options"]], "label": c["label"],
                       "axes": c.get("axes", [])},  # fmt: skip
    "ora": lambda c: {"words": c["words"], "backgrounds": c["backgrounds"],
                      "featured": c["featured"]},
    "subset": lambda c: {"words": c["words"], "options": c["options"],
                         "scales": c.get("scales", []), "summary": c["summary"]},
    "genelist": lambda c: {"words": c["words"], "contrasts": c["contrasts"],
                           "featured": c["featured"]},
}  # fmt: skip


def defaults(controls: Mapping[str, Any]) -> dict[str, Any]:
    """What the controls show before the reader moves them: the default rule, no paint, or all.

    A filter's live line is ``ui.toml`` [filter] summary, or summary_count where the figure
    prints no median, filled here with the numbers of every bar and passed unfilled for the
    page script to fill again (``ctxFilter.summary``). Each select is labelled in words.
    """
    if controls["kind"] in ("survival", "follow"):
        # The survival panel's own defaults, built by ctx.clinical.page.live; a figure that
        # follows the panel has none.
        return dict(controls.get("defaults", {}))
    if (shaped := _SHAPED.get(controls["kind"])) is not None:
        return shaped(controls)
    if controls["kind"] == "pairs":
        # The select's label: the figure's own where it gives one, M8's F8 naming miRNA pairs.
        return {
            "options": list(enumerate(controls["options"])),
            "summary": controls["summary"],
            "label": controls.get("label"),
        }
    if controls["kind"] == "brush":
        return {"colours": list(zip(controls["colours"], controls["names"], strict=True)),
                "empty": controls["empty"]}  # fmt: skip
    if controls["kind"] == "filter":
        line = "filter.summary" if controls.get("summary_median", True) else "filter.summary_count"
        labels = controls.get("filter_labels", {})
        return {"filters": [(c, levels, labels.get(c, c))
                            for c, levels in controls["filters"].items()],
                "summary": text.fill("ui", line, **controls["start"]),
                "template": text.raw("ui", line),
                "none": text.word("ui", "filter.none")}  # fmt: skip
    return {
        "fc": f"{controls['lfc_min']:g}",
        "bases": list(controls["bases"].items()),
        "alpha": f"{controls['alpha']:g}",
        "presets": sorted({*controls["presets"], f"{controls['min_cutoff']:g}"}, key=float),
        "summary": controls["default_summary"],
    }
