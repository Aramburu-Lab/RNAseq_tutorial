"""M1c: turn one checked fetchngs outdir into the harmonised sample sheet and its page."""

import hashlib
from dataclasses import replace
from datetime import date
from pathlib import Path

import numpy as np
import polars as pl

from ctx import runtime
from ctx.dataset import Dataset, DatasetError, load
from ctx.gather import rnaseq
from ctx.gather.curation import PROBLEM_COLUMNS, CurationError, curate
from ctx.gather.fetchngs import Outdir, read
from ctx.io.naming import MODULES, stem
from ctx.io.readme_companion import write_table
from ctx.report.builder import (
    FigureText,
    Group,
    KeyNumber,
    ModulePage,
    ParameterCard,
    Reproduction,
    Result,
    ResultTable,
    SelfCheck,
    write,
)
from ctx.viz.bar import BarData
from ctx.viz.save import figure_block, write_figure
from ctx.viz.style import category_colours

MODULE = "gather"
# Below this, "0.0 million reads" would hide the depth: print the reads themselves.
DEPTH_IN_MILLIONS = 100_000
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
DERIVED_COLUMNS = {
    "replicate": (
        "none",
        "Order of this sample among the runs of one donor, cell type and condition: by the "
        "dataset file's [order], then by run accession; from 1.",
    ),
    "batch": ("none", "Batch from the dataset file's read-length bins ([derived.batch])."),
    "read_length": ("nt", "base_count / read_count / mates (2 for paired runs); never rounded."),
}
COLUMNS = {**OWN_COLUMNS, **ENA_COLUMNS, **CURATION_COLUMNS}
INT32_MAX = 2**31 - 1
# fetchngs's FTP path keeps only the last two ENA files of a run (sra_runinfo_to_ftp.py:70-71).
FTP_KEEPS = 2


def columns(spec: Dataset | None) -> dict[str, tuple[str, str]]:
    """(unit, meaning) of every sheet column: fixed columns, then the dataset's extras."""
    if spec is None:
        return COLUMNS
    extras = {c.name: (c.unit, c.meaning) for c in spec.columns}
    fixed = {**COLUMNS, **{k: v for k, v in extras.items() if k in CURATION}}
    fixed |= {k: DERIVED_COLUMNS[k] for k in ("replicate", "batch")}
    return {**fixed, **{k: v for k, v in extras.items() if k not in CURATION},
            "read_length": DERIVED_COLUMNS["read_length"]}  # fmt: skip


def sheet(outdir: Outdir, curated: pl.DataFrame | None = None) -> pl.DataFrame:
    """The harmonised sheet: one row per run, the chosen columns, curation filled if given.

    Args:
        outdir: A checked fetchngs outdir.
        curated: ``curation.curate()`` output, or None to leave the curation columns empty.

    Returns:
        The sheet. ``read_count`` and ``base_count`` are Int32 where every value fits,
        Int64 otherwise: base counts of deep runs pass 2**31.
    """
    runs = outdir.runs
    frame = runs.select([c for c in COLUMNS if c in runs.columns]).with_columns(
        [pl.lit("").alias(c) for c in CURATION]
    )
    if curated is not None:
        frame = frame.drop(CURATION).join(
            curated, on="run_accession", how="left", maintain_order="left"
        )
    counts = frame.select(pl.col("read_count", "base_count").cast(pl.Int64))
    # Downcast only after checking the range (CONVENTIONS.md, Stack), column by column.
    fitted = [
        counts[c].cast(pl.Int32)
        if int(counts.select(pl.col(c).max()).item() or 0) <= INT32_MAX
        else counts[c]
        for c in counts.columns
    ]
    extras = [c for c in frame.columns if c not in COLUMNS]
    return frame.with_columns(fitted).select(*COLUMNS, *extras)


def _statuses(frame: pl.DataFrame) -> dict[str, int]:
    """FASTQ files per status, over both mates."""
    statuses = frame.select(pl.concat_list("status_1", "status_2").explode().alias("s")).filter(
        pl.col("s") != ""
    )
    return dict(sorted(statuses["s"].value_counts().iter_rows()))


def _keys(outdir: Outdir, frame: pl.DataFrame) -> tuple[KeyNumber, ...]:
    """Runs, experiments and files: the first three things to check after a download."""
    statuses = _statuses(frame)
    detail = ", ".join(f"{n} {k}" for k, n in statuses.items()) or "none"
    return (
        KeyNumber("Runs", f"{frame.height:,}"),
        KeyNumber("Experiments", f"{frame['sample'].n_unique():,}"),
        KeyNumber("FASTQ files", f"{sum(statuses.values()):,}", detail),
        KeyNumber("Unlisted FASTQs", f"{len(outdir.unlisted_fastq)}", "in fastq/, in no row"),
    )


def _facts(outdir: Outdir, frame: pl.DataFrame, curated: bool) -> ResultTable:
    """What the outdir holds, checked, one row per check."""
    multi = frame.filter(pl.col("n_runs") > 1)["sample"].unique().sort().to_list()
    wide = frame.filter(pl.col("ena_files") > FTP_KEEPS)["run_accession"].to_list()
    counted = ", ".join(f"{k} {n}" for k, n in _statuses(frame).items())
    rows = [
        ("Source", f"nf-core/fetchngs {outdir.provenance['fetchngs']} (commit "
                   f"{outdir.provenance['fetchngs_commit']}), Nextflow "
                   f"{outdir.provenance['nextflow']}"),
        ("FASTQ files by status", counted or "none"),
        ("Experiments with more than one run", ", ".join(multi) or "none"),
        ("Runs with more than two ENA files", ", ".join(wide) or "none"),
        ("FASTQ files that no row names", ", ".join(outdir.unlisted_fastq) or "none"),
    ]  # fmt: skip
    if not curated:
        rows.append(("Curation", "Not curated: without --dataset-file, donor, cell type, "
                                 "condition, replicate and batch stay empty."))  # fmt: skip
    return ResultTable("Download checks", ("Check", "Result"), tuple(rows))


QUESTION = (
    "Which runs did the download bring, are the files intact, and what do we know about "
    "each sample?"
)
CONCEPT = (
    "nf-core/fetchngs downloads the reads and writes a sample sheet, but a finished download "
    "can still be incomplete. Files go missing or arrive corrupted, one experiment can hold "
    "several runs, and public metadata is free text. gather checks the download against its "
    "own records and the published checksums, then writes the one sample sheet every later "
    "module reads."
)
TERMS = (
    ("Run", "One sequencing run, one accession such as SRR or ERR, with one or two FASTQ files."),
    ("Experiment", "The library a run was sequenced from; one experiment can hold several runs."),
    ("no_md5", "Fetched through sra-tools, which publishes no checksum to verify against."),
)
PRACTICE = (
    "Run gather after every download and resolve what its checks report before any analysis. "
    "From then on, treat its sample sheet as the only record of what each sample is."
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


def _composition_text(data: BarData) -> FigureText:
    """Title, lead and reading guide of the composition bars, from the counts drawn."""
    groups, batches = data.frame["group"].n_unique(), data.frame["batch"].unique().sort()
    reading = [
        "Each bar is one cell type and condition; its length counts the runs in it.",
        "Hover a bar for the individuals behind it; one individual can give several cell "
        "types, so individuals do not add up across bars.",
    ]
    if batches.len() > 1:
        reading.insert(
            1,
            f"Colours split each bar by batch ({', '.join(batches)}): a group "
            "that sits in one batch alone cannot be told apart from that batch.",
        )
    sizes = data.frame.group_by("group").agg(pl.col("n_samples").sum())
    smallest = sizes.sort("n_samples", "group").row(0)
    return FigureText(
        title="Samples by cell type and condition",
        lead=f"{groups} groups; the smallest, {smallest[0]}, holds {smallest[1]} runs, "
        "which limits what any contrast on it can find.",
        reading=tuple(reading),
    )  # fmt: skip


def _library_text(data: BarData) -> FigureText:
    """Title, lead and reading guide of the library sizes, from the depths drawn."""
    median = float(np.median(data.frame["read_count"].to_numpy()))
    depth = f"{median / 1e6:.1f} million" if median >= DEPTH_IN_MILLIONS else f"{median:,.0f}"
    batched = data.frame["batch"].n_unique() > 1
    return FigureText(
        title="Sequencing depth per run",
        lead=f"The median run holds {depth} reads.",
        reading=(
            "Each bar is one run; its height is the reads ENA reports for it, in millions.",
            "Runs are grouped by batch and sorted by depth; colour marks the batch."
            if batched
            else "Runs are sorted by depth.",
            "A run far below the others"
            + (" in its batch" if batched else "")
            + " may be a failed library: DESeq2's size factors correct depth, not quality.",
        ),
    )


def _composition(frame: pl.DataFrame, spec: Dataset) -> BarData:
    """Samples per cell type and condition, stacked by batch, individuals on hover."""
    groups = (
        frame.with_columns(batch=pl.col("batch").replace("", "all"))
        .group_by("cell_type", "condition", "batch")
        .agg(n_samples=pl.len().cast(pl.Int64), n_individuals=pl.col("donor").n_unique())
        .sort("cell_type", "condition", "batch")
        .with_columns(group=pl.format("{} {}", "cell_type", "condition"))
    )
    levels = [b for b in (*spec.batch_bins, "all") if b in set(groups["batch"])]
    columns = {
        "cell_type": ("none", "Cell type, curated."),
        "condition": ("none", "Condition, curated."),
        "batch": ("none", "Batch from read length; all when the dataset declares none."),
        "n_samples": ("runs", "Runs in this cell type, condition and batch."),
        "n_individuals": ("individuals", "Distinct donors among those runs."),
        "group": ("none", "Cell type and condition, the bar's label."),
    }
    caption = (
        f"{frame.height} runs from {frame['donor'].n_unique()} individuals, by cell type, "
        "condition and batch; counts, no test. An individual can give several cell types, "
        "so individuals do not add up across bars."
    )
    return BarData(
        frame=groups, columns=columns, category="group", value="n_samples",
        category_label="Cell type and condition", value_label="Runs", series="batch",
        colours=category_colours(levels), caption=caption, horizontal=True, stacked=True,
        hover=("n_individuals",), width="double",
    )  # fmt: skip


def _library_sizes(frame: pl.DataFrame, spec: Dataset) -> BarData:
    """Reads per run, grouped by batch, largest first."""
    runs = frame.with_columns(batch=pl.col("batch").replace("", "all")).sort(
        "batch", "read_count", descending=[False, True]
    )
    levels = [b for b in (*spec.batch_bins, "all") if b in set(runs["batch"])]
    medians = "; ".join(
        f"{b}: n = {n} runs, median {m:,.0f} reads"
        for b, n, m in runs.group_by("batch", maintain_order=True)
        .agg(pl.len(), pl.col("read_count").median())
        .rows()
    )
    columns = {
        "run_accession": COLUMNS["run_accession"],
        "read_count": COLUMNS["read_count"],
        "million_reads": ("million reads", "read_count / 1,000,000, the bar's length."),
        "batch": ("none", "Batch from read length; all when the dataset declares none."),
        "cell_type": ("none", "Cell type, curated."),
        "condition": ("none", "Condition, curated."),
    }
    return BarData(
        frame=runs.with_columns(
            pl.col("read_count").cast(pl.Int64), million_reads=pl.col("read_count") / 1e6
        ).select(*columns),
        columns=columns, category="run_accession", value="million_reads",
        category_label="Run, by batch and depth", value_label="Reads", series="batch",
        colours=category_colours(levels),
        caption=f"Reads per run as ENA reports them, no test. {medians}.",
        stacked=True, hover=("cell_type", "condition"), width="double",
    )  # fmt: skip


def _curation_facts(frame: pl.DataFrame, spec: Dataset) -> ResultTable:
    """What curation joined and derived, counted."""
    assert spec.sdrf is not None  # noqa: S101 - curation ran, so the dataset has an SDRF
    lengths = ", ".join(
        f"{v:g} nt {n}" for v, n in frame["read_length"].value_counts().sort("read_length").rows()
    )
    batches = ", ".join(
        f"{b or 'all'} {n}" for b, n in frame["batch"].value_counts().sort("batch").rows()
    )
    checked = ", ".join(a for a, _ in spec.sdrf.cross_check)
    repeats = frame.filter(pl.col("replicate") > 1).height
    joined = (
        f"Curation from {spec.path.name}: each of {frame.height} runs joined to one SDRF row "
        f"on {spec.sdrf.join[0]}, in both directions"
        + (f"; {checked} agrees for every run" if checked else "")
        + "; every value is in the dataset file's value maps."
    )
    rows = (
        ("Join", joined),
        ("Read length", lengths),
        ("Batch", batches),
        ("Repeated runs", f"{repeats} {'run' if repeats == 1 else 'runs'} beyond the first "
                          "of a donor, cell type and condition"),
    )  # fmt: skip
    return ResultTable("Curation", ("Check", "Result"), rows)


def _counts_tables(
    counts: rnaseq.Counts,
) -> list[tuple[str, pl.DataFrame, dict[str, tuple[str, str]]]]:
    """(description, frame, column dictionary) of the counts and genes tables."""
    source = f"nf-core/rnaseq {counts.provenance['nf-core/rnaseq']} {rnaseq.MATRIX_FILE}"
    per_sample = {
        s: ("reads", f"Length-scaled count for sample {s}, rounded half to even ({source}).")
        for s in counts.matrix.columns[1:]
    }
    gene_id = ("none", "Gene identifier from the annotation the index was built with.")
    return [
        ("counts", counts.matrix, {"gene_id": gene_id, **per_sample}),
        ("genes", counts.genes, {"gene_id": gene_id,
                                 "gene_name": ("none", "Gene symbol from the same annotation.")}),
    ]  # fmt: skip


def _counts_facts(counts: rnaseq.Counts) -> ResultTable:
    """What the counts step read and did, counted."""
    genes, samples = counts.matrix.height, counts.matrix.width - 1
    kind = counts.matrix.dtypes[1] if samples else "none"
    rows = (
        ("Matrix", f"{genes:,} genes by {samples} samples, the length-scaled gene matrix"),
        ("Source", f"nf-core/rnaseq {counts.provenance['nf-core/rnaseq']}, Salmon "
                   f"{counts.provenance['salmon']}"),
        ("Rounded", f"{counts.fractional:,} of {genes * samples:,} values were fractional and "
                    "were rounded half to even"),
        ("Stored as", f"{kind}, {counts.matrix.estimated_size():,} bytes in memory"),
    )  # fmt: skip
    return ResultTable("Counts", ("Check", "Result"), rows)


def _write_problems(
    error: CurationError, out: Path, dataset: str, releases: dict[str, str]
) -> Path:
    """Write the problem table beside where the sheet would have gone; no sheet is written."""
    problems = pl.DataFrame(error.problems, schema=list(PROBLEM_COLUMNS), orient="row")
    problem_stem = stem(MODULE, f"{dataset}_curation_problems", 1, day=date.today())
    written = write_table(
        problems, out, problem_stem, columns=PROBLEM_COLUMNS, data_releases=releases
    )
    return written[0]


def _curate(
    settings: runtime.Settings, spec: Dataset, checked: Outdir, out: Path, releases: dict[str, str]
) -> pl.DataFrame:
    """Curate the runs with the dataset's SDRF, adding what was read to ``releases``.

    Raises:
        DatasetError: If the SDRF has not been fetched under the root.
        CurationError: After writing the problem table.
    """
    assert spec.sdrf is not None  # noqa: S101 - run() refuses a dataset file without [sdrf]
    sdrf = settings.root / spec.resource(spec.sdrf.resource).dest
    if not sdrf.is_file():
        raise DatasetError(f"{sdrf} is missing; run ctx setup --dataset-file {spec.path}")
    releases |= {
        "dataset file": str(spec.path),
        f"sha256:{spec.path.name}": hashlib.sha256(spec.path.read_bytes()).hexdigest(),
        f"sha256:{sdrf.name}": hashlib.sha256(sdrf.read_bytes()).hexdigest(),
    }
    try:
        return curate(checked.runs, sdrf.read_text(), spec)
    except CurationError as error:
        _write_problems(error, out, spec.slug, releases)
        raise


def run(  # noqa: PLR0913 - one call from the CLI, every argument is a named option or context
    settings: runtime.Settings,
    cards: tuple[ParameterCard, ...],
    command: str,
    dataset: str,
    *,
    fetchngs_dir: Path | None = None,
    dataset_file: Path | None = None,
    rnaseq_dir: Path | None = None,
) -> list[Path]:
    """Check the outdir, curate it if a dataset file is given, and write the sheet and page.

    Args:
        settings: The resolved runtime.
        cards: Parameter cards for every option of this command.
        command: The command line, as the reader should retype it.
        dataset: Slug naming the dataset in output names.
        fetchngs_dir: The fetchngs outdir; unset means ``<root>/data/raw/<dataset>``.
        dataset_file: The dataset file; unset writes the sheet with empty curation columns.
        rnaseq_dir: An nf-core/rnaseq outdir whose gene counts to check and write; unset
            writes no counts.

    Returns:
        The page, then the sheet files, then the figure files.

    Raises:
        DatasetError: If the dataset file is unusable, names another slug, or its SDRF has
            not been fetched under the root.
        CurationError: After writing the problem table, if curation found any problem.
        rnaseq.RnaseqError: If the rnaseq outdir cannot be used; nothing is written.
    """
    day = date.today()
    sheet_stem = stem(MODULE, f"{dataset}_samplesheet", 1, day=day)  # validates the slug first
    page_stem = stem(MODULE, f"{dataset}_overview", 1, day=day)
    spec = load(dataset_file) if dataset_file else None
    if spec and spec.slug != dataset:
        raise DatasetError(f"{dataset_file} is for {spec.slug!r}, not --dataset {dataset!r}")
    if spec and spec.sdrf is None:
        raise DatasetError(f"{dataset_file}: no [sdrf] table; nothing to curate from")
    source = (fetchngs_dir or settings.root / "data" / "raw" / dataset).resolve()
    checked = read(source)  # raises before anything is written
    out = settings.root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}"
    releases = {"fetchngs outdir": str(source), **checked.provenance}
    curated = _curate(settings, spec, checked, out, releases) if spec else None
    frame = sheet(checked, curated)
    # Read before anything is written, so a bad outdir leaves no sheet behind.
    counts = rnaseq.read(rnaseq_dir.resolve(), set(frame["sample"])) if rnaseq_dir else None
    files = write_table(frame, out, sheet_stem, columns=columns(spec), data_releases=releases)
    results: list[Result] = [_facts(checked, frame, curated is not None)]
    if counts is not None:
        counts_releases = {**releases, "rnaseq outdir": str(rnaseq_dir), **counts.provenance}
        for name, table, described in _counts_tables(counts):
            table_stem = stem(MODULE, f"{dataset}_{name}", 1, day=day)
            files += write_table(table, out, table_stem, columns=described,
                                 data_releases=counts_releases)  # fmt: skip
        results.append(_counts_facts(counts))
    if spec and curated is not None:
        results.append(_curation_facts(frame, spec))
        composition, depth = _composition(frame, spec), _library_sizes(frame, spec)
        figures = []
        for name, data, words in (
            ("composition", composition, _composition_text(composition)),
            ("library_sizes", depth, _library_text(depth)),
        ):
            figure_stem = stem(MODULE, f"{dataset}_{name}", 1, day=day)
            text = replace(words, back=f"{page_stem}.html")
            files += write_figure(data, out, figure_stem, title=text.title,
                                  data_releases=releases, text=text)  # fmt: skip
            figures.append(figure_block(data, figure_stem, text))
        results.append(Group("Samples", tuple(figures)))
    usage = runtime.measure(settings.started)
    budget = runtime.budget(MODULE)
    page = ModulePage(
        title=f"M1 Gather: {dataset}",
        question=QUESTION,
        concept=CONCEPT,
        parameters=cards,
        results=tuple(results),
        practice=PRACTICE,
        think=THINK,
        mistakes=MISTAKES,
        key_numbers=_keys(checked, frame),
        terms=TERMS,
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
