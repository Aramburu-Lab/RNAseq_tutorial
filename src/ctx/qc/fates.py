"""Read fates from an nf-core/rnaseq outdir: what Trim Galore did with each read, then Salmon.

MultiQC reads the copied Trim Galore JSON reports and Salmon folders (D-025); nothing here
parses a tool log. A paired sample appears as its read-1 and read-2 reports and a summed
row; the sample takes read 1, which counts pairs, as Salmon does. Every read must have a fate:
Trim Galore's fates add up to its total, and the reads it wrote are the fragments Salmon
processed; either failing stops the run, naming the samples.
"""

from collections.abc import Iterable
from pathlib import Path

import polars as pl

from ctx.qc import QcError, multiqc_io
from ctx.qc.registry import Registry

MODULES = ("trim_galore", "salmon")
TRIM = ("total_reads", "reads_written", "reads_too_short", "reads_too_long", "reads_too_many_n",
        "reads_discarded_untrimmed")  # fmt: skip
SALMON = ("num_processed", "num_mapped", "num_decoy_fragments")
# Read 1 of a pair (rank 0), then a single-end sample's own row (1); a pair's summed row (2)
# and its read 2 (3) are never read.
FIRST_UNUSED_RANK = 2


def ranks(samples: Iterable[str]) -> dict[str, tuple[str, int]]:
    """Every name MultiQC gives a sample, with the rank ``multiqc_io.by_unit`` reads it at."""
    out: dict[str, tuple[str, int]] = {}
    for sample in samples:
        out |= {f"{sample}_trimmed_1": (sample, 0), sample: (sample, 1),
                f"{sample}_trimmed": (sample, 2), f"{sample}_trimmed_2": (sample, 3)}  # fmt: skip
    return out


def run_multiqc(
    rnaseq_dir: Path, dest: Path, config: Path, samples: list[str]
) -> multiqc_io.MultiqcRun:
    """MultiQC on the outdir's Trim Galore and Salmon folders; every sample must have both.

    Checked on MultiQC's own tables, not its sources file: MultiQC 1.35's Salmon module lists
    no sources there (measured on the rnaseq fixture), so a missing report shows as a
    missing row.
    """
    folders = [rnaseq_dir / "trimgalore", rnaseq_dir / "salmon"]
    if missing := [str(f) for f in folders if not f.is_dir()]:
        raise QcError(f"the rnaseq outdir has no {missing}; was it copied with step 7?")
    run = multiqc_io.run(folders, dest, config=config, modules=MODULES)
    names = {name: sample for name, (sample, _) in ranks(samples).items()}
    for module in MODULES:
        found = multiqc_io.map_names(multiqc_io.table(run, module)["Sample"], names)
        if missing := sorted(set(samples) - set(found.values())):
            raise QcError(f"no {module} report for {missing}")
    return run


def values(run: multiqc_io.MultiqcRun, samples: list[str], registry: Registry) -> pl.DataFrame:
    """One row per sample: every trim and quant registry metric, by its column."""
    metrics = [m for m in registry.metrics.values()
               if m.stage in ("trim", "quant") and m.table == "general_stats"]  # fmt: skip
    table = multiqc_io.table(run, "general_stats")
    if missing := [m.column for m in metrics if m.column not in table.columns]:
        raise QcError(f"MultiQC's general_stats table has no {missing}")
    numbers = [m.column for m in metrics if m.kind in ("numeric", "strandedness")]
    table = multiqc_io.numeric(table, numbers).select("Sample", *(m.column for m in metrics))
    return multiqc_io.by_unit(table, ranks(samples), "sample", FIRST_UNUSED_RANK)


def fates(run: multiqc_io.MultiqcRun, samples: list[str]) -> pl.DataFrame:
    """Each sample's reads by fate: Trim Galore's, then Salmon's for the reads it wrote.

    Raises:
        QcError: If a sample's Trim Galore fates do not add up to its total, or Salmon
            processed other than the reads Trim Galore wrote.
    """
    trim = multiqc_io.numeric(multiqc_io.table(run, "trim_galore"), TRIM)
    salmon = multiqc_io.numeric(multiqc_io.table(run, "salmon"), SALMON)
    frame = multiqc_io.by_unit(trim.select("Sample", *TRIM), ranks(samples), "sample",
                               FIRST_UNUSED_RANK).join(
        multiqc_io.by_unit(salmon.select("Sample", *SALMON), ranks(samples), "sample",
                           FIRST_UNUSED_RANK), on="sample", how="left")  # fmt: skip
    frame = frame.with_columns(
        not_assigned=pl.col("num_processed") - pl.col("num_mapped") - pl.col("num_decoy_fragments")
    )
    lost = frame.filter(pl.col("total_reads") != pl.sum_horizontal(TRIM[1:]))["sample"]
    unequal = frame.filter(pl.col("num_processed") != pl.col("reads_written"))["sample"]
    if lost.len() or unequal.len():
        raise QcError(f"reads without a fate: Trim Galore's fates miss its total for "
                      f"{lost.to_list()}; Salmon processed other than the reads written for "
                      f"{unequal.to_list()}")  # fmt: skip
    return frame.select(pl.col("sample"), pl.exclude("sample").cast(pl.Int64))
