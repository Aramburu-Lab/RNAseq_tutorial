"""The tables M3 writes so ``ctx figure`` can redraw its figures from tables alone (F.1).

Each holds what a figure drew and no other table holds: the reads per outcome of every bar,
as the aligner's reports gave them, the gene sums each scatter draws, the sheet groups the
bars filter by, and what the page states beside them. They are written after the run's own
tables and before any figure, so a figure's spec finds them.
"""

from collections.abc import Mapping
from typing import TYPE_CHECKING

import polars as pl

from ctx.de.run import Writer
from ctx.process import Aligner
from ctx.qc.registry import Registry
from ctx.report import text

if TYPE_CHECKING:
    from ctx.process.page import Frames, Shown

PAGE = "process"
OUTCOME_SCHEMA = {"label": pl.String, "outcome": pl.String, "reads": pl.Int64}


def outcome_rows(reads: Mapping[str, Mapping[str, int]]) -> pl.DataFrame:
    """Each bar's reads per outcome and its input reads, in the order the bars stand."""
    rows = [(label, k, v) for label, r in reads.items() for k, v in r.items()]
    return pl.DataFrame(rows, schema=OUTCOME_SCHEMA, orient="row")


def outcomes_from(frame: pl.DataFrame) -> dict[str, dict[str, int]]:
    """``outcome_rows`` read back as ``page.Frames`` holds them."""
    found: dict[str, dict[str, int]] = {}
    for label, outcome, reads in frame.rows():
        found.setdefault(label, {})[outcome] = reads
    return found


def options_row(shown: "Shown", registry: Registry) -> pl.DataFrame:
    """What the page states beside the tables: aligner, featured sample, genome, registry."""
    return pl.DataFrame({
        "aligner": [str(shown.aligner)], "featured": [shown.featured],
        "genome_name": [shown.genome_name], "genome_genes": [shown.genome_genes],
        "partial": [shown.partial], "registry": [str(registry.path)],
    })  # fmt: skip


def aligner_of(options: pl.DataFrame) -> Aligner:
    """The aligner ``options_row`` recorded."""
    return Aligner(options["aligner"][0])


def write(writer: Writer, found: "Frames", shown: "Shown", registry: Registry) -> None:
    """Every table above, its columns described in ``content/process.toml``."""
    named = {
        "figure_options": options_row(shown, registry),
        "outcomes": outcome_rows(found.outcomes),
        "setting_outcomes": outcome_rows(found.settings),
        "sums": found.agreements.sums,
        "groups": shown.groups,
    }
    for what, frame in named.items():
        units = text.raw(PAGE, f"tables.{what}.units")
        writer.table(frame, what, {c: (units.get(c, "none"),
                                       text.word(PAGE, f"tables.{what}.columns.{c}"))
                                   for c in frame.columns})  # fmt: skip
