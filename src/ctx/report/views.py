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


def _plain(figure: Any) -> dict[str, Any]:
    """A plotly figure as the page receives it: arrays as lists, NaN as null."""
    from plotly.utils import PlotlyJSONEncoder  # noqa: PLC0415 - heavy

    loaded: dict[str, Any] = json.loads(json.dumps(figure.to_plotly_json(), cls=PlotlyJSONEncoder))
    return loaded


def viewed(block: FigureBlock, options: Sequence[tuple[str, Any, FigureText]],
           label: str) -> FigureBlock:  # fmt: skip
    """A figure block with a select between ``options``: (name, contract, words) each.

    The first option must be the block's own figure. Its words stay as built; another
    option's title, lead and caption replace them while it is shown.

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
    controls: dict[str, Any] = {
        "kind": "view", "label": label, "options": shown,
        "moved": text.raw("ui", "view.moved"), "build": options[0][0],
    }  # fmt: skip
    if hasattr(first, "export_layers"):
        controls["export"] = export.settings(first.inches(), first.caption, block.stem,
                                             first.export_layers())  # fmt: skip
    return replace(block, controls=controls)
