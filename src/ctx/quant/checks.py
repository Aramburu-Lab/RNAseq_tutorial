"""What must hold before M4 writes a table, a figure or a page: each check lists every offender.

- Every fragment has one fate: Salmon's mapped, decoy and below-score fragments never exceed
  what it processed, so the no-hit remainder is never negative.
- Length scaling keeps each sample's reads: a length-scaled column sums to Salmon's mapped
  fragments (tximport scales to the sum of the estimated counts, which is ``num_mapped`` up
  to quant.sf's rounding).
- ctx's gene tables are nf-core's: the port of ``ctx.quant.tximport`` on every ``quant.sf``
  equals each merged gene table R wrote, gene for gene.
- The ceiling: the matrix counts no gene the annotation lacks, so no more genes than it holds.
"""

from collections.abc import Sequence
from pathlib import Path

import numpy as np
import polars as pl

from ctx.quant import CountsFrom, QuantError, outdir
from ctx.quant.tximport import GeneTables

# quant.sf prints NumReads to three decimals, so the estimated counts tximport scales to can
# miss num_mapped by half a thousandth per transcript in the index, and no more.
ROUNDING = 0.0005
# R writes its tables with 15 significant digits.
TABLE_TOLERANCE = 1e-9
SHOWN = 5
# Each merged gene table and the port's table it must equal.
NF_CORE_TABLES = {"gene_counts": CountsFrom.NO, "gene_counts_length_scaled":
                  CountsFrom.LENGTH_SCALED_TPM, "gene_counts_scaled": CountsFrom.SCALED_TPM,
                  "gene_tpm": "tpm", "gene_lengths": "lengths"}  # fmt: skip


def _names(frame: pl.DataFrame) -> str:
    """The first samples of a frame, and how many more."""
    names = frame["sample"].to_list()
    return f"{names[:SHOWN]}{f' and {len(names) - SHOWN} more' if len(names) > SHOWN else ''}"


def fate_problems(fates: pl.DataFrame) -> list[str]:
    """Samples whose mapped, decoy and below-score fragments exceed those processed."""
    lost = fates.filter(pl.col("no_hit") < 0)
    return [f"fragments with two fates, more than Salmon processed: {_names(lost)}"] \
        if lost.height else []  # fmt: skip


def sum_problems(sums: pl.DataFrame, targets: pl.DataFrame) -> list[str]:
    """Samples whose length-scaled column sum is not Salmon's mapped fragments.

    Args:
        sums: ``tables.sums``.
        targets: sample and targets, the transcripts in each sample's index.
    """
    off = sums.join(targets, on="sample").filter(
        (pl.col("length_scaled_sum") - pl.col("mapped")).abs() > ROUNDING * pl.col("targets")
    )
    return [f"length-scaled counts that do not sum to the mapped fragments, beyond quant.sf's "
            f"rounding: {_names(off)}"] if off.height else []  # fmt: skip


def table_values(ours: GeneTables, what: CountsFrom | str) -> np.ndarray:
    """The port's table that a merged table holds."""
    if isinstance(what, CountsFrom):
        return ours.counts_from(what)
    values: np.ndarray = getattr(ours, what)
    return values


def tximport_problems(ours: GeneTables, folder: Path) -> list[str]:
    """Each merged gene table of nf-core's that ctx's port does not reproduce.

    Args:
        ours: The port's tables, from every sample's ``quant.sf``.
        folder: The salmon folder holding ``salmon.merged.*.tsv``.
    """
    problems = []
    for table, what in NF_CORE_TABLES.items():
        theirs = outdir.merged(folder, table, outdir.GENE_KEYS, ours.samples)
        if theirs["gene_id"].to_list() != list(ours.genes):
            problems.append(f"{table}: its genes are not the port's, in the port's order")
            continue
        expected, found = theirs.select(ours.samples).to_numpy(), table_values(ours, what)
        if not np.allclose(found, expected, rtol=TABLE_TOLERANCE, atol=TABLE_TOLERANCE):
            worst = float(np.max(np.abs(found - expected) / np.maximum(np.abs(expected), 1)))
            problems.append(f"{table}: the port differs from nf-core's by up to {worst:.3g}")
    return problems


def ceiling_problems(counted: Sequence[str], annotated: Sequence[str]) -> list[str]:
    """Genes the matrix counts that the annotation lacks; none may exceed the reference."""
    known = set(annotated)
    if stray := [g for g in counted if g not in known]:
        return [f"{len(stray)} counted genes are not in the annotation: {stray[:SHOWN]}"]
    return []


def check(problems: Sequence[str]) -> None:
    """Stop with every failed check together.

    Raises:
        QuantError: If any check failed.
    """
    if problems:
        raise QuantError("ctx quant checks failed; nothing written:\n  " + "\n  ".join(problems))
