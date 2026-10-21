"""The ``view`` control kind (D-051): a select between named alternatives Python computed.

Lookup only: each option is a figure contract drawn in Python, with its own title, lead and
caption, and the page swaps it in. The first option is the figure as built, which the PDF and
SVG files draw; another one's files come from the page through the export (D-033), with the
live line naming the view. M5's p-value histogram by contrast and M6's gene set enrichment by
contrast use it (work-order items 6 and 7).
"""

import json
from collections.abc import Sequence
from dataclasses import replace
from typing import Any

from ctx.report import export, marks, text
from ctx.report.builder import FigureBlock, FigureText
from ctx.report.figures import plot_width


def _plain(figure: Any) -> dict[str, Any]:
    """A plotly figure as the page receives it: arrays as lists, NaN as null."""
    from plotly.utils import PlotlyJSONEncoder  # noqa: PLC0415 - heavy

    loaded: dict[str, Any] = json.loads(json.dumps(figure.to_plotly_json(), cls=PlotlyJSONEncoder))
    return loaded


def own_widths(shown: Sequence[dict[str, Any]]) -> None:
    """Give each option after the first whose plot box differs from the first's its ``width``.

    The page sets the shown option's (``ctx-view`` in ``_page_fit.js.j2``), and the first's
    for an option without one: a view of 4 then 11 features kept the 4's least width and the
    11 turned names overlapped 10 times (V.4b). Options of one width carry nothing new.
    """
    built = plot_width(shown[0]["figure"])
    for option in shown[1:]:
        if (width := plot_width(option["figure"])) != built:
            option["width"] = width


def viewed(block: FigureBlock, options: Sequence[tuple[str, Any, FigureText]],
           label: str) -> FigureBlock:  # fmt: skip
    """A figure block with a select between ``options``: (name, contract, words) each.

    The first option must be the block's own figure. Its words stay as built; another
    option's title, lead and caption replace them while it is shown. An option whose plot box
    differs from the first's, as a box of 11 features against 4, carries its own width; one
    whose files differ in caption, size or layers carries their export settings, so its PDF
    and SVG are the view as set (V.4b). An option equal to the build in all three carries none.

    Raises:
        ValueError: With fewer than two options, or a first option that is not the block's.
    """
    if len(options) < 2:  # noqa: PLR2004 - a choice needs two
        raise ValueError("a view needs two or more options")
    first = options[0][1]
    if first.caption_text() != block.caption:
        raise ValueError("the first option is the figure as built")
    shown = [{"label": name, "figure": _plain(data.interactive()), "title": words.title,
              "lead": str(marks.html(words.lead)), "caption": str(marks.html(data.caption_text()))}
             for name, data, words in options]  # fmt: skip
    own_widths(shown)
    controls: dict[str, Any] = {
        "kind": "view", "label": label, "options": shown,
        "moved": text.raw("ui", "view.moved"), "build": options[0][0],
    }  # fmt: skip
    if hasattr(first, "export_layers"):
        controls["export"] = _settings(first, block.stem)
        for option, (_, data, _) in zip(shown[1:], options[1:], strict=True):
            own = _settings(data, block.stem)
            if own != controls["export"]:
                option["export"] = own
    return replace(block, controls=controls)


def _settings(data: Any, stem: str) -> dict[str, Any]:
    """One option's PDF and SVG settings, from the size and layers its own files take."""
    return export.settings(data.inches(), data.caption, stem, data.export_layers())
