"""The tables M8 writes so ``ctx figure`` can redraw its figures from tables alone (F.1).

The run's options as ``ctx ncrna`` was given them, with the annotation's file name the first
caption prints, and the annotation's record of each count row with its biotype group, which
no other table holds whole: the pairs read its places. Both are written before any figure, so
a figure's spec finds them.
"""

from typing import TYPE_CHECKING

import polars as pl

from ctx.de.run import Writer
from ctx.report import text

if TYPE_CHECKING:
    from ctx.ncrna.page import Context

PAGE = "ncrna"
# A list option as one cell: the biotype groups of --biotypes, in order.
SEPARATOR = ","


def options_row(ctx: "Context") -> pl.DataFrame:
    """The biotype groups, the filter's minimum count and the annotation's file name."""
    return pl.DataFrame({"biotypes": [SEPARATOR.join(ctx.biotypes)],
                         "min_count": [ctx.min_count], "gtf": [ctx.inputs.gtf.name]})  # fmt: skip


def biotypes_of(options: pl.DataFrame) -> tuple[str, ...]:
    """The biotype groups ``options_row`` recorded, in the run's order."""
    return tuple(options["biotypes"][0].split(SEPARATOR))


def write(writer: Writer, ctx: "Context") -> None:
    """Both tables, their columns described in ``content/ncrna.toml``."""
    named = {
        "figure_options": (options_row(ctx), text.raw(PAGE, "tables.figure_options")),
        "gene_biotypes": (ctx.inputs.genes, text.raw(PAGE, "tables.gene_biotypes")),
    }
    for what, (frame, words) in named.items():
        writer.table(frame, what, {c: (words["units"].get(c, "none"), words["columns"][c])
                                   for c in frame.columns})  # fmt: skip
