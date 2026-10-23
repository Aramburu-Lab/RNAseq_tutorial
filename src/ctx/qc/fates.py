"""Read fates from an nf-core/rnaseq outdir: what Trim Galore did with each read, then Salmon.

nf-core's own MultiQC wrote the tables read here (``multiqc_io.nfcore``, D-116); ctx runs no
MultiQC over the outdir and parses no tool log. A paired sample appears as its read-1 and
read-2 rows (``<sample>_1``, ``<sample>_2``; ``<sample> Read 1`` in General Statistics) and a
merged row under its own name; the sample takes read 1, which counts pairs, as Salmon does.
Every read must have a fate: Trim Galore's fates add up to its total, and Salmon processed no
more than the reads Trim Galore wrote, since a step between them (BBSplit, SortMeRNA) may
remove reads; either failing stops the run, naming the samples.
"""

from collections.abc import Iterable

import polars as pl

from ctx.qc import QcError, multiqc_io
from ctx.qc.registry import Registry

# nf-core's tables M2 reads: Trim Galore, Salmon alone, and FastQC on the raw reads, which
# nf-core's MultiQC config anchors as fastqc_raw (multiqc.log of the test profile, 2026-10-23).
MODULES = ("trim_galore", "salmon", "fastqc_fastqc_raw")
TRIM = ("total_reads", "reads_written", "reads_too_short", "reads_too_long", "reads_too_many_n",
        "reads_discarded_untrimmed")  # fmt: skip
SALMON = ("num_processed", "num_mapped", "num_decoy_fragments")
# Read 1 of a pair (rank 0), then a single-end sample's own row (1); a pair's merged row (2)
# and its read 2 (3) are never read. A paired sample's own name holds its merged row only in
# the tables that merge mates (General Statistics, FastQC); Salmon and STAR, which count
# fragments, write it once under that name, so it is read there at rank 1.
FIRST_UNUSED_RANK = 2


def ranks(samples: Iterable[str]) -> dict[str, tuple[str, int]]:
    """Every name nf-core's MultiQC gives a sample, with the rank ``multiqc_io.by_unit`` reads.

    ``<sample>_1`` and ``<sample> Read 1`` are read 1 of a pair, ``<sample>_2`` and
    ``<sample> Read 2`` read 2, ``<sample>`` a single-end sample or a pair's merged or
    fragment row.
    """
    out: dict[str, tuple[str, int]] = {}
    for sample in samples:
        out |= {f"{sample}_1": (sample, 0), f"{sample} Read 1": (sample, 0), sample: (sample, 1),
                f"{sample}_2": (sample, 3), f"{sample} Read 2": (sample, 3)}  # fmt: skip
    return out


def first_names(samples: Iterable[str]) -> list[str]:
    """The names of read 1 or a single-end sample's one file, as FastQC's plots name them."""
    return [n for n, (_, rank) in ranks(samples).items() if rank < FIRST_UNUSED_RANK]


def checked(run: multiqc_io.MultiqcRun, samples: list[str]) -> multiqc_io.MultiqcRun:
    """nf-core's MultiQC tables, every sample in each table M2 reads, or refused naming them.

    Raises:
        QcError: If a table is missing, a row names no sample, or a sample has no row.
    """
    names = {name: sample for name, (sample, _) in ranks(samples).items()}
    for module in MODULES:
        found = multiqc_io.map_names(multiqc_io.table(run, module)["Sample"], names)
        if missing := sorted(set(samples) - set(found.values())):
            raise QcError(f"no {module} report for {missing}")
    return run


def values(run: multiqc_io.MultiqcRun, samples: list[str], registry: Registry) -> pl.DataFrame:
    """One row per sample: every trim, quant and raw-read FastQC registry metric, by its column.

    Read 1 of a pair stands for the sample.
    """
    metrics = [m for m in registry.metrics.values()
               if m.stage in ("reads", "trim", "quant") and m.module is not None]  # fmt: skip
    frame = pl.DataFrame({"sample": samples})
    for name in dict.fromkeys(m.table for m in metrics):
        table = multiqc_io.table(run, name)
        wanted = [m for m in metrics if m.table == name]
        if missing := [m.column for m in wanted if m.column not in table.columns]:
            raise QcError(f"nf-core's MultiQC {name} table has no {missing}")
        numbers = [m.column for m in wanted if m.kind in ("numeric", "strandedness")]
        table = multiqc_io.numeric(table, numbers).select("Sample", *(m.column for m in wanted))
        part = multiqc_io.by_unit(table, ranks(samples), "sample", FIRST_UNUSED_RANK)
        frame = frame.join(part, on="sample", how="left", maintain_order="left")
    return frame


def fates(run: multiqc_io.MultiqcRun, samples: list[str]) -> pl.DataFrame:
    """Each sample's reads by fate: Trim Galore's, then Salmon's for the reads it processed.

    Salmon's reads not processed are those a step between trimming and Salmon removed.

    Raises:
        QcError: If a sample's Trim Galore fates do not add up to its total, or Salmon
            processed more than the reads Trim Galore wrote.
    """
    trim = multiqc_io.numeric(multiqc_io.table(run, "trim_galore"), TRIM)
    salmon = multiqc_io.numeric(multiqc_io.table(run, "salmon"), SALMON)
    frame = multiqc_io.by_unit(trim.select("Sample", *TRIM), ranks(samples), "sample",
                               FIRST_UNUSED_RANK).join(
        multiqc_io.by_unit(salmon.select("Sample", *SALMON), ranks(samples), "sample",
                           FIRST_UNUSED_RANK), on="sample", how="left",
        maintain_order="left")  # fmt: skip
    frame = frame.with_columns(
        not_assigned=pl.col("num_processed") - pl.col("num_mapped") - pl.col("num_decoy_fragments"),
        removed_before_salmon=pl.col("reads_written") - pl.col("num_processed"),
    )
    lost = frame.filter(pl.col("total_reads") != pl.sum_horizontal(TRIM[1:]))["sample"]
    gained = frame.filter(pl.col("num_processed") > pl.col("reads_written"))["sample"]
    if lost.len() or gained.len():
        raise QcError(f"reads without a fate: Trim Galore's fates miss its total for "
                      f"{lost.to_list()}; Salmon processed more than the reads written for "
                      f"{gained.to_list()}")  # fmt: skip
    return frame.select(pl.col("sample"), pl.exclude("sample").cast(pl.Int64))
