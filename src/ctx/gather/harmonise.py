"""M1c: turn one checked fetchngs outdir into the harmonised sample sheet and its page."""

from datetime import date
from pathlib import Path

import polars as pl

from ctx import runtime
from ctx.gather.fetchngs import Outdir, read
from ctx.io.naming import MODULES, stem
from ctx.io.readme_companion import write_table
from ctx.report.builder import ModulePage, ParameterCard, Reproduction, SelfCheck, write

MODULE = "gather"
ENA_SOURCE = "ENA Portal API returnFields, result=read_run, retrieved 2026-09-26"
# (unit, meaning) per column; ENA meanings quote ENA's own field description.
ENA_COLUMNS = {
    "study_accession": ("none", f"study accession number ({ENA_SOURCE})"),
    "sample_accession": ("none", f"sample accession number ({ENA_SOURCE})"),
    "library_layout": ("none", f"sequencing library layout ({ENA_SOURCE})"),
    "library_selection": (
        "none",
        "method used to select or enrich the material being sequenced "
        f"({ENA_SOURCE}); carried to M8",
    ),
    "library_source": ("none", f"source material being sequenced ({ENA_SOURCE})"),
    "library_strategy": (
        "none",
        f"sequencing technique intended for the library ({ENA_SOURCE}); carried to M8",
    ),
    "instrument_model": ("none", f"instrument model used in sequencing experiment ({ENA_SOURCE})"),
    "read_count": ("reads", f"number of reads ({ENA_SOURCE})"),
    "base_count": ("base pairs", f"number of base pairs ({ENA_SOURCE})"),
    "fastq_bytes": (
        "bytes",
        f"size (in bytes) of FASTQ files, semicolon-separated per file ({ENA_SOURCE})",
    ),
    "scientific_name": ("none", f"scientific name of an organism ({ENA_SOURCE})"),
    "sample_title": ("none", f"brief sample title ({ENA_SOURCE}); free text, input to curation"),
    "experiment_title": ("none", f"brief experiment title ({ENA_SOURCE})"),
}
OWN_COLUMNS = {
    "sample": (
        "none",
        "Experiment accession: fetchngs's sample, shared by every run of one experiment.",
    ),
    "run_accession": ("none", f"accession number ({ENA_SOURCE}); one row per run."),
    "run_index": ("none", "Order of this run within its sample, by numeric accession, from 1."),
    "n_runs": (
        "runs",
        "Runs sharing this sample: technical replicates, combined where FASTQs are consumed.",
    ),
    "single_end": ("none", "True when ENA reports the library layout as SINGLE."),
    "fastq_1": ("none", "First FASTQ, relative to the fetchngs outdir; empty when not downloaded."),
    "status_1": (
        "none",
        "verified, no_md5 (sra-tools path) or not_downloaded; mismatch and missing stop gather.",
    ),
    "fastq_2": ("none", "Second FASTQ of a paired run, relative to the outdir; empty otherwise."),
    "status_2": ("none", "As status_1, for fastq_2; empty for single-end runs."),
    "strandedness": (
        "none",
        "auto, as fetchngs writes it; M4 fills it from Salmon's lib_format_counts.json.",
    ),
    "ena_files": (
        "files",
        "FASTQ files ENA lists for the run; above 2, FTP keeps only the last two.",
    ),
}
CURATION = ("donor", "cell_type", "condition", "replicate", "batch")
CURATION_COLUMNS = {
    c: ("none", f"{c}, from the dataset's curation file (Phase 2); empty until then.")
    for c in CURATION
}
COLUMNS = {**OWN_COLUMNS, **ENA_COLUMNS, **CURATION_COLUMNS}
INT32_MAX = 2**31 - 1
# fetchngs's FTP path keeps only the last two ENA files of a run (sra_runinfo_to_ftp.py:70-71).
FTP_KEEPS = 2


def sheet(outdir: Outdir) -> pl.DataFrame:
    """The harmonised sheet: one row per run, the chosen columns, curation columns empty.

    Args:
        outdir: A checked fetchngs outdir.

    Returns:
        The sheet, in ``COLUMNS`` order. ``read_count`` and ``base_count`` are Int32 where
        every value fits, Int64 otherwise: base counts of deep runs pass 2**31.
    """
    runs = outdir.runs
    frame = runs.select([c for c in COLUMNS if c in runs.columns]).with_columns(
        [pl.lit("").alias(c) for c in CURATION]
    )
    counts = frame.select(pl.col("read_count", "base_count").cast(pl.Int64))
    # Downcast only after checking the range (CONVENTIONS.md, Stack), column by column.
    fitted = [
        counts[c].cast(pl.Int32)
        if int(counts.select(pl.col(c).max()).item() or 0) <= INT32_MAX
        else counts[c]
        for c in counts.columns
    ]
    return frame.with_columns(fitted).select(list(COLUMNS))


def _facts(outdir: Outdir, frame: pl.DataFrame) -> tuple[str, ...]:
    """What the outdir holds, counted, one sentence each."""
    statuses = frame.select(pl.concat_list("status_1", "status_2").explode().alias("s")).filter(
        pl.col("s") != ""
    )
    counted = ", ".join(f"{k} {v}" for k, v in sorted(statuses["s"].value_counts().iter_rows()))
    multi = frame.filter(pl.col("n_runs") > 1)["sample"].unique().sort().to_list()
    wide = frame.filter(pl.col("ena_files") > FTP_KEEPS)["run_accession"].to_list()
    return (
        f"{frame.height} runs in {frame['sample'].n_unique()} experiments, from nf-core/fetchngs "
        f"{outdir.provenance['fetchngs']} (commit {outdir.provenance['fetchngs_commit']}), "
        "Nextflow "
        f"{outdir.provenance['nextflow']}.",
        f"FASTQ files by status: {counted or 'none'}.",
        f"Experiments with more than one run: {', '.join(multi) or 'none'}.",
        f"Runs with more than two ENA files: {', '.join(wide) or 'none'}.",
        f"FASTQ files in fastq/ that no row names: {', '.join(outdir.unlisted_fastq) or 'none'}.",
        "Curation columns (donor, cell_type, condition, replicate, batch) are empty: the curation "
        "file arrives in Phase 2.",
    )


QUESTION = (
    "Which runs did the download bring, are the files intact, and what do we know about "
    "each sample?"
)
CONCEPT = (
    "nf-core/fetchngs downloads reads and writes a samplesheet, but it is not the end of the "
    "job. Files can be missing or corrupted, one experiment can hold several runs, and public "
    "metadata is free text. gather checks the outdir against itself and against the md5 "
    "sums, and writes one sheet every later module reads."
)
PRACTICE = (
    "Run gather after every fetch. A run marked no_md5 came through sra-tools and was not "
    "checked against ENA; several runs under one sample are technical replicates of one "
    "library, not biological replicates."
)
THINK = (
    SelfCheck(
        "Why does gather resolve FASTQ files by name instead of the path in samplesheet.csv?",
        "fetchngs writes the path as it was given on its own machine; after a copy to a laptop "
        "or from Zenodo that path points nowhere.",
    ),
    SelfCheck(
        "Why is a run with no md5 file reported, not rejected?",
        "Runs with no FASTQ at ENA are fetched through sra-tools, which writes no md5; the file "
        "is still usable, but it was not checked against a published checksum.",
    ),
)
MISTAKES = (
    "Treating the runs of one experiment as biological replicates.",
    "Re-running fetchngs into an existing outdir: it rebuilds the samplesheet from the new "
    "input only.",
)


def run(  # noqa: PLR0913 - one call from the CLI, every argument is a named option or context
    settings: runtime.Settings,
    cards: tuple[ParameterCard, ...],
    command: str,
    dataset: str,
    fetchngs_dir: Path | None,
) -> list[Path]:
    """Check the outdir, then write the sheet and its page under ``outputs/01_gather``.

    Args:
        settings: The resolved runtime.
        cards: Parameter cards for every option of this command.
        command: The command line, as the reader should retype it.
        dataset: Slug naming the dataset in output names.
        fetchngs_dir: The fetchngs outdir; unset means ``<root>/data/raw/<dataset>``.

    Returns:
        The page, then the sheet files.
    """
    day = date.today()
    sheet_stem = stem(MODULE, f"{dataset}_samplesheet", 1, day=day)  # validates the slug first
    page_stem = stem(MODULE, f"{dataset}_overview", 1, day=day)
    source = (fetchngs_dir or settings.root / "data" / "raw" / dataset).resolve()
    checked = read(source)  # raises before anything is written
    frame = sheet(checked)
    out = settings.root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}"
    releases = {"fetchngs outdir": str(source), **checked.provenance}
    files = write_table(frame, out, sheet_stem, columns=COLUMNS, data_releases=releases)
    usage = runtime.measure(settings.started)
    budget = runtime.budget(MODULE)
    page = ModulePage(
        title=f"M1c Gather: {dataset}",
        question=QUESTION,
        concept=CONCEPT,
        parameters=cards,
        results=_facts(checked, frame),
        practice=PRACTICE,
        think=THINK,
        mistakes=MISTAKES,
        reproduce=Reproduction(
            command=command,
            wall_seconds=usage.wall_seconds,
            peak_rss_self=usage.peak_rss_self,
            peak_rss_children=usage.peak_rss_children,
            profile=f"{settings.profile}, {settings.threads} threads",
            budget=f"Track B: under {budget.minutes} min and {budget.ram_gb} GB",
            versions=runtime.versions(),
        ),
    )
    return [write(page, out / f"{page_stem}.html"), *files]
