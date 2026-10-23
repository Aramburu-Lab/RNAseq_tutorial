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
# nf-core/rnaseq's FastQC on the raw reads, where the outdir's copy holds it (P2Q.9): MultiQC
# then judges every sample's reads by FastQC's own verdicts, not the teaching runs' alone.
RAW_FASTQC = Path("fastqc") / "raw"
TRIM =("total_reads", "reads_written", "reads_too_short", "reads_too_long", "reads_too_many_n",
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
                f"{sample}_trimmed": (sample, 2), f"{sample}_trimmed_2": (sample, 3),
                **raw_ranks(sample)}  # fmt: skip
    return out


def nfcore_ranks(samples: Iterable[str]) -> dict[str, tuple[str, int]]:
    """Every name nf-core's own MultiQC gives a sample, ranked as ``ranks`` ranks (D-116).

    ``<sample>_1`` and ``<sample> Read 1`` are read 1 of a pair, ``<sample>_2`` and
    ``<sample> Read 2`` read 2, ``<sample>`` a single-end sample or a pair's merged or
    fragment row.
    """
    out: dict[str, tuple[str, int]] = {}
    for sample in samples:
        out |= {f"{sample}_1": (sample, 0), f"{sample} Read 1": (sample, 0), sample: (sample, 1),
                f"{sample}_2": (sample, 3), f"{sample} Read 2": (sample, 3)}  # fmt: skip
    return out


def raw_ranks(sample: str) -> dict[str, tuple[str, int]]:
    """The names nf-core's raw-read FastQC reports take in MultiQC, ranked as ``ranks`` ranks.

    ``<sample>_raw`` for a single-end sample, ``<sample>_raw_1`` and ``_raw_2`` for a pair.
    """
    return {f"{sample}_raw_1": (sample, 0), f"{sample}_raw": (sample, 1),
            f"{sample}_raw_2": (sample, 3)}  # fmt: skip


def has_raw_fastqc(rnaseq_dir: Path) -> bool:
    """Whether the outdir's copy holds nf-core's FastQC reports on the raw reads."""
    return any((rnaseq_dir / RAW_FASTQC).glob("*_fastqc.zip"))


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
    raw = has_raw_fastqc(rnaseq_dir)
    modules = (*MODULES, "fastqc") if raw else MODULES
    run = multiqc_io.run([*folders, *[rnaseq_dir / RAW_FASTQC] * raw], dest, config=config,
                         modules=modules)  # fmt: skip
    names = {name: sample for name, (sample, _) in ranks(samples).items()}
    for module in modules:
        found = multiqc_io.map_names(multiqc_io.table(run, module)["Sample"], names)
        if missing := sorted(set(samples) - set(found.values())):
            raise QcError(f"no {module} report for {missing}")
    return run


def values(run: multiqc_io.MultiqcRun, samples: list[str], registry: Registry) -> pl.DataFrame:
    """One row per sample: every trim and quant registry metric, by its column.

    Where MultiQC read the raw-read FastQC reports too, FastQC's metrics join them, its
    general statistics and its module verdicts, read 1 of a pair standing for the sample.
    """
    raw = (run.data / "multiqc_fastqc.txt").is_file()
    metrics = [m for m in registry.metrics.values()
               if (m.stage in ("trim", "quant") and m.table == "general_stats")
               or (raw and m.module == "fastqc")]  # fmt: skip
    frame = pl.DataFrame({"sample": samples})
    for name in dict.fromkeys(m.table for m in metrics):
        table = multiqc_io.table(run, name)
        wanted = [m for m in metrics if m.table == name]
        if missing := [m.column for m in wanted if m.column not in table.columns]:
            raise QcError(f"MultiQC's {name} table has no {missing}")
        numbers = [m.column for m in wanted if m.kind in ("numeric", "strandedness")]
        table = multiqc_io.numeric(table, numbers).select("Sample", *(m.column for m in wanted))
        part = multiqc_io.by_unit(table, ranks(samples), "sample", FIRST_UNUSED_RANK)
        frame = frame.join(part, on="sample", how="left", maintain_order="left")
    return frame


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
                           FIRST_UNUSED_RANK), on="sample", how="left",
        maintain_order="left")  # fmt: skip
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
