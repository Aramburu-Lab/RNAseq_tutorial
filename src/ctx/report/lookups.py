"""Views the page looks up (D-073): each option light and with its own files; and a grid.

The controls are the ``view`` kind's, as ``views.viewed`` builds them, with three changes.
The page swaps in its own theme's template whenever it draws an option (``viewing`` in
``_page_controls.js.j2``), so an option's copy of the light template, 8,376 of a p-value
histogram's 9,678 B, is dropped. Every option after the first carries the export settings its
PDF and SVG take, ``export``: its own caption, size and thinned points, where the build's
caption under another view's bars would state numbers the bars do not show, and the build's
thinned points would land in another view's traces; an option whose plot box differs from
the build's also carries its ``width`` (``views.own_widths``). And a grid names one option by
a select per setting, as M5's lesson counts by adjustment, cutoff and fold change: ``axes``
lists each setting's label, choices and the build's choice; ``grid`` the option of each combination,
keyed by ``key``, whose page twin is ``ctxGrid.key`` (``_grid.js.j2``).
"""

import json
from collections.abc import Mapping, Sequence
from dataclasses import replace
from typing import Any

from ctx.report import export, marks, text
from ctx.report.builder import FigureBlock, FigureText
from ctx.report.views import own_widths


def _plain(figure: Any) -> dict[str, Any]:
    """A plotly figure as the page receives it: arrays as lists, NaN as null."""
    from plotly.utils import PlotlyJSONEncoder  # noqa: PLC0415 - heavy

    loaded: dict[str, Any] = json.loads(json.dumps(figure.to_plotly_json(), cls=PlotlyJSONEncoder))
    return loaded


def key(choices: Sequence[int]) -> str:
    """The grid key of one combination: each axis's choice position, joined by commas."""
    return ",".join(str(c) for c in choices)


def _exported(data: Any, stem: str) -> dict[str, Any] | None:
    """The export settings of one option's PDF and SVG; None for a contract with no export."""
    if not hasattr(data, "export_layers"):
        return None
    return export.settings(data.inches(), data.caption, stem, data.export_layers())


def lookup(block: FigureBlock, options: Sequence[tuple[str, Any, FigureText]], label: str,
           axes: Sequence[tuple[str, Sequence[str]]] = (),
           at: Sequence[Sequence[int]] = ()) -> FigureBlock:  # fmt: skip
    """A view select over ``options``: (name, contract, words) each, the build's first.

    Args:
        block: The figure as built.
        options: The views; the first is the block's own figure, whose words stay as built.
        label: The select's label; with ``axes``, the form's.
        axes: (label, choices) of each setting of a grid; empty for one select.
        at: With ``axes``, each option's choice position on every axis, in option order.

    Returns:
        The block with its controls.

    Raises:
        ValueError: With fewer than two options; a first option that is not the block's; or a
            grid whose positions do not name each combination of its axes once.
    """
    if len(options) < 2:  # noqa: PLR2004 - a choice needs two
        raise ValueError("a view needs two or more options")
    if options[0][1].caption_text() != block.caption:
        raise ValueError("the first option is the figure as built")
    shown = []
    for i, (name, data, words) in enumerate(options):
        drawn = _plain(data.interactive())
        drawn["layout"] = {k: v for k, v in drawn["layout"].items() if k != "template"}
        shown.append({"label": name, "figure": drawn, "title": words.title,
                      "lead": str(marks.html(words.lead)),
                      "caption": str(marks.html(data.caption_text())),
                      **({"export": _exported(data, block.stem)} if i else {})})  # fmt: skip
    own_widths(shown)
    controls: dict[str, Any] = {
        "kind": "view", "label": label, "options": shown,
        "moved": text.raw("ui", "view.moved"), "build": options[0][0],
    }  # fmt: skip
    if (first := _exported(options[0][1], block.stem)) is not None:
        controls["export"] = first
    if axes:
        controls.update(grid_of(axes, at, len(options)))
    return replace(block, controls=controls)


def grid_of(axes: Sequence[tuple[str, Sequence[str]]], at: Sequence[Sequence[int]],
            n: int) -> Mapping[str, Any]:  # fmt: skip
    """The ``axes`` and ``grid`` entries of a grid's controls; see ``lookup``.

    Raises:
        ValueError: If the positions are not ``n`` distinct combinations covering the axes.
    """
    sizes = [len(choices) for _, choices in axes]
    if len(at) != n or any(len(p) != len(sizes) for p in at):
        raise ValueError("a grid names every option by one choice on every axis")
    if any(not 0 <= c < s for p in at for c, s in zip(p, sizes, strict=True)):
        raise ValueError("a grid position lies outside its axis")
    grid = {key(p): i for i, p in enumerate(at)}
    total = 1
    for s in sizes:
        total *= s
    if len(grid) != n or n != total:
        raise ValueError(f"a grid of {total} combinations names {len(grid)} of {n} options")
    return {"axes": [{"label": a, "choices": list(c), "start": at[0][i]}
                     for i, (a, c) in enumerate(axes)],
            "grid": grid}  # fmt: skip
