"""The tables M4 writes so ``ctx figure`` can redraw its figures from tables alone (F.1).

The sheet as the figures name it, and the featured contrast's samples with their level, as
M5 picks them: written with the per-sample tables, before any figure, so a spec finds them.
"""

from typing import TYPE_CHECKING

from ctx.quant import blocks as blk

if TYPE_CHECKING:
    from ctx.quant.page import Context


def write(ctx: "Context", found: blk.Shown) -> None:
    """The sheet the figures read, and the contrast's samples where a contrast is featured."""
    from ctx.quant.page import described  # noqa: PLC0415 - page imports this module

    named = {"figure_sheet": found.sheet}
    if ctx.inputs.contrast is not None:
        named["contrast_samples"] = blk.contrast_samples(ctx.inputs)
    for what, frame in named.items():
        ctx.writer.table(frame, what, described(frame))
