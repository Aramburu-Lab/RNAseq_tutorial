"""M2 QC (`ctx qc`): every sample judged by the registry, and the judgement tested.

The run reads the gather outputs, never the fetchngs layout (CONVENTIONS), measures each
stage it has inputs for, judges every metric the registry holds against its bands and its
stratum, writes one status row per run, and asks whether a metric, or a flag, differs between
the groups of a comparison M5 makes. Nothing is excluded.
"""

import json
from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from datetime import date
from pathlib import Path
from typing import Any

import polars as pl

from ctx import runtime
from ctx.dataset import Dataset, DatasetError, load
from ctx.de.run import Writer, newest
from ctx.gather import gdc
from ctx.io.naming import MODULES, stem
from ctx.qc import Fastp, QcError, checkpoint, fastqc, fates, multiqc_io, reads
from ctx.qc import registry as reg
from ctx.qc.blocks import (
    Reads,
    Stage,
    fastqc_figures,
    fate_figures,
    lines,
    overview_groups,
    read_figures,
    strips,
)
from ctx.report import text
from ctx.report.builder import (
    Group,
    KeyNumber,
    ModulePage,
    ParameterCard,
    Reproduction,
    Result,
    ResultTable,
    write,
)

MODULE = "qc"
PAGE = MODULE
GTF_RESOURCE = "annotation_gtf"
# The page lists this many confounding tests, smallest adjusted p first; the CSV holds all.
SHOWN_TESTS = 10
ALPHA = 0.05
# M3's tables M2 judges on the align stage (M3.12): `ctx process` writes one row per sample.
ALIGN_TABLES = ("alignment", "strandedness")


@dataclass(frozen=True)
class Options:
    """What `ctx qc` was asked to do, beside the dataset.

    Attributes:
        dataset_file: The dataset file: contrasts, level orders and resources.
        registry: The registry file; None resolves it (``registry.resolve``).
        strata: Sheet columns whose shared values make a stratum, for the MAD bands.
        peers: Sheet columns whose shared values make a peer group.
        mad: The MAD multiplier; None keeps the registry's.
        rnaseq_dir: An nf-core/rnaseq outdir whose Trim Galore and Salmon reports to read.
        reads: Run FastQC, fastp and seqkit on the sheet's FASTQ files (M2a).
        fetchngs_dir: The fetchngs outdir the sheet's FASTQ paths are relative to; None is
            ``<root>/data/raw/<dataset>``.
        fastp: fastp's filter for the reads and the grid.
        grid_run: The run the -q × -l grid is drawn on; None takes the first read.
    """

    dataset_file: Path
    registry: Path | None = None
    strata: tuple[str, ...] = ("cell_type",)
    peers: tuple[str, ...] = ("cell_type", "condition")
    mad: float | None = None
    rnaseq_dir: Path | None = None
    reads: bool = False
    fetchngs_dir: Path | None = None
    fastp: Fastp = field(default_factory=Fastp)
    grid_run: str | None = None


def load_registry(options: Options) -> reg.Registry:
    """The registry, with the MAD multiplier the reader chose."""
    registry = reg.load(reg.resolve(options.registry))
    return registry if options.mad is None else replace(registry, mad_multiplier=options.mad)


def gene_types(root: Path, spec: Dataset) -> tuple[pl.DataFrame | None, str]:
    """The annotation's gene types, or None and why: no such resource, or not fetched."""
    found = [r for r in spec.resources if r.name == GTF_RESOURCE]
    if not found:
        return None, text.fill(PAGE, "notes.no_annotation", resource=GTF_RESOURCE)
    path = root / found[0].dest
    if not path.is_file():
        return None, text.fill(PAGE, "notes.not_fetched", path=str(path))
    return checkpoint.gene_types(path), ""


def counts_stage(
    settings: runtime.Settings, spec: Dataset, gathered: Mapping[str, pl.DataFrame],
    options: Options,
) -> tuple[checkpoint.Measured | None, tuple[str, ...]]:  # fmt: skip
    """The count checkpoint, and what was left out: all of it where gather wrote no counts."""
    if "counts" not in gathered:
        return None, (text.word(PAGE, "notes.no_counts"),)
    types, note = gene_types(settings.root, spec)
    measured = checkpoint.sample_metrics(gathered["counts"], gathered["samplesheet"],
                                         gathered["genes"], types, options.peers,
                                         settings.threads)  # fmt: skip
    return measured, (note,) if note else ()


def gathered_tables(root: Path, dataset: str) -> dict[str, Path]:
    """Gather's newest sheet, and its counts and genes where it wrote them (with --rnaseq-dir)."""
    folder = root / "outputs" / f"{dict(MODULES)['gather']:02d}_gather"
    inputs = {"samplesheet": newest(folder, dataset, "samplesheet")}
    for what in ("counts", "genes"):
        try:
            inputs[what] = newest(folder, dataset, what)
        except DatasetError:
            continue
    return inputs


def read_gathered(inputs: Mapping[str, Path]) -> dict[str, pl.DataFrame]:
    """Gather's tables, the sheet through the GDC door and the counts cut to its samples.

    A GDC sheet lists every file; ``gdc.chosen``, the door M5 reads it through, keeps one per
    case and sample type, so both modules count the same samples (TODO S.1: 517 LUAD tumours
    here against M5's 515). The counts are read for the kept samples' columns alone.
    """
    sheet = gdc.chosen(pl.read_parquet(inputs["samplesheet"]))
    samples = sheet.unique("sample", keep="first", maintain_order=True)["sample"].to_list()
    gathered = {"samplesheet": sheet}
    if "counts" in inputs:
        gathered["counts"] = pl.read_parquet(inputs["counts"], columns=["gene_id", *samples])
    if "genes" in inputs:
        gathered["genes"] = pl.read_parquet(inputs["genes"])
    return gathered


def align_tables(root: Path, dataset: str) -> dict[str, Path]:
    """The newest `ctx process` tables of the dataset that the align stage reads, if any."""
    folder = root / "outputs" / f"{dict(MODULES)['process']:02d}_process"
    found = {}
    for what in ALIGN_TABLES:
        try:
            found[what] = newest(folder, dataset, what)
        except DatasetError:
            continue
    return found


def align_stage(tables: Mapping[str, Path], registry: reg.Registry) -> pl.DataFrame | None:
    """M3's per-sample alignment values, by the registry's align-stage columns (M3.12).

    Each table gives the columns of the registry's align metrics it holds; None without any.
    """
    frame = None
    for what, path in tables.items():
        wanted = [m.column for m in registry.stage("align") if m.table == what]
        table = pl.read_parquet(path)
        if missing := [c for c in wanted if c not in table.columns]:
            raise QcError(f"{path.name} has no {missing}, which the registry's align stage reads")
        part = table.select("sample", *[pl.col(c).cast(pl.Float64) for c in wanted])
        frame = part if frame is None else frame.join(part, on="sample", how="full",
                                                      coalesce=True)  # fmt: skip
    return frame


def rnaseq_stage(
    writer: Writer, rnaseq_dir: Path, samples: list[str], registry: reg.Registry
) -> tuple[pl.DataFrame, pl.DataFrame, fastqc.Plots | None]:
    """The pipeline's reports through MultiQC: values per sample, the fates, FastQC's plots.

    The MultiQC config the registry implies is written first, beside the tables, so
    MultiQC's own report colours its table by the same bands. FastQC's plots are None where
    the outdir's copy holds no raw-read FastQC reports.
    """
    mqc = fates.run_multiqc(rnaseq_dir, writer.out / "multiqc", multiqc_config(writer, registry),
                            samples)  # fmt: skip
    plots = None
    if fates.has_raw_fastqc(rnaseq_dir):
        names = {name: s for name, (s, _) in fates.ranks(samples).items()}
        first = [n for n, (_, r) in fates.ranks(samples).items() if r < fates.FIRST_UNUSED_RANK]
        plots = fastqc.read(mqc, names, first, "sample")
    return fates.values(mqc, samples, registry), fates.fates(mqc, samples), plots


def multiqc_config(writer: Writer, registry: reg.Registry) -> Path:
    """The MultiQC config the registry implies, written once beside the tables."""
    config = writer.out / f"{writer.stem('multiqc_config')}.yml"
    if config not in writer.files:
        writer.out.mkdir(parents=True, exist_ok=True)
        settings = reg.multiqc_config(registry, text.word(PAGE, "multiqc.comment"))
        config.write_text(json.dumps(settings, indent=1) + "\n")
        writer.files.append(config)
    return config


def reads_stage(writer: Writer, settings: runtime.Settings, dataset: str, sheet: pl.DataFrame,
                registry: reg.Registry, options: Options) -> Reads:  # fmt: skip
    """M2a: FastQC, fastp and seqkit on the sheet's FASTQ files, and fastp's grid on one run.

    MultiQC reads the reads into their own report, ``multiqc/reads/``, apart from the
    pipeline's: FastQC names samples by file, and one table would mix the two naming rules.
    """
    outdir = options.fetchngs_dir or settings.root / "data" / "raw" / dataset
    files, skipped = reads.fastq_files(sheet, outdir)
    logs = writer.out / "multiqc" / "reads" / "logs"
    folders = reads.run_tools(files, logs, threads=settings.threads, fastp=options.fastp)
    config = multiqc_config(writer, registry)
    mqc = multiqc_io.run(folders, writer.out / "multiqc" / "reads", config=config,
                         modules=reads.MODULES)  # fmt: skip
    grid_run = options.grid_run or str(files["run_accession"][0])
    paths = files.filter(pl.col("run_accession") == grid_run).sort("mate")["path"].to_list()
    if not paths:
        raise QcError(f"--grid-run {grid_run} is not a run whose reads were read")
    grid = reads.fastp_grid(paths, logs / "grid", writer.out / "multiqc" / "grid",
                            grid=(reads.GRID_Q, reads.GRID_L), fastp=options.fastp,
                            threads=settings.threads)  # fmt: skip
    names = {reads.stem(p): r for r, p in files.select("run_accession", "path").rows()}
    first = [reads.stem(p) for p in files.filter(pl.col("mate") == 1)["path"]]
    return Reads(reads.run_values(mqc, files, registry), skipped,
                 fastqc.read(mqc, names, first, "run_accession"), grid, grid_run)  # fmt: skip


def judge(values: pl.DataFrame, sheet: pl.DataFrame, found: Reads | None,
          registry: reg.Registry, options: Options, notes: tuple[str, ...]) -> Stage:  # fmt: skip
    """Every metric judged: per sample, and the reads per run, both within their strata."""
    samples = sheet.unique("sample", keep="first", maintain_order=True)
    bands, flags = checkpoint.judge(values, samples, registry, options.strata)
    if found is None:
        return Stage(values, bands, flags, notes)
    by_run = found.values.drop("sample").rename({"run_accession": "sample"})
    runs = sheet.select(
        pl.col("run_accession").alias("sample"), *options.strata, checkpoint.LENGTH_GROUP
    )
    run_bands, run_flags = checkpoint.judge(by_run, runs, registry, options.strata)
    return Stage(values, pl.concat([bands, run_bands], how="diagonal_relaxed"), flags, notes,
                 found.values, run_flags)  # fmt: skip


def confounding(stage: Stage, sheet: pl.DataFrame, spec: Dataset,
                registry: reg.Registry) -> pl.DataFrame:  # fmt: skip
    """The confounding tests on every sample-level value and flag, reads included.

    A sample takes its first run's read metrics and read flags, as M5's selection keeps one
    run per individual.
    """
    values, flags = stage.values, stage.flags
    if stage.run_values is not None and stage.run_flags is not None:
        # A metric measured on the sample, from the pipeline's FastQC, keeps that value.
        own = stage.run_values.drop("sample", *(c for c in values.columns if c != "sample"),
                                    strict=False)  # fmt: skip
        values = values.join(by_first_run(own, sheet, "run_accession"), on="sample",
                             how="left")  # fmt: skip
        flags = pl.concat([flags, by_first_run(stage.run_flags, sheet, "sample")
                           .select(flags.columns)])  # fmt: skip
    return checkpoint.confounding(values, flags, sheet, spec.contrasts, spec.order, registry)


def by_first_run(run_frame: pl.DataFrame, sheet: pl.DataFrame, key: str) -> pl.DataFrame:
    """A per-run table keyed by sample, each sample taking the run the sheet lists first."""
    first = sheet.unique("sample", keep="first", maintain_order=True).select("run_accession",
                                                                             "sample")  # fmt: skip
    return (
        run_frame.rename({key: "run_accession"})
        .join(first, on="run_accession", how="inner")
        .drop("run_accession")
    )


def columns_of(frame: pl.DataFrame, registry: reg.Registry) -> dict[str, tuple[str, str]]:
    """(unit, meaning) of every column of a values table, from the registry's own words."""
    by_column = {m.column: (m.unit, m.meaning) for m in registry.metrics.values()}
    own = {c: (u, text.word(PAGE, f"columns.{c}")) for c, u in
           (("sample", "none"), ("run_accession", "none"))}  # fmt: skip
    return {c: own.get(c) or by_column[c] for c in frame.columns}


def _write_tables(writer: Writer, values: pl.DataFrame, frames: Mapping[str, pl.DataFrame],
                  registry: reg.Registry) -> None:  # fmt: skip
    """The values table, its columns described by the registry, then every other table.

    Each other table's columns, with their units, are described in ``content/qc.toml``; a
    column it does not describe stops the run.
    """
    writer.table(values, "metrics", columns_of(values, registry))
    for what, frame in frames.items():
        units = text.raw(PAGE, f"tables.{what}.units")
        writer.table(frame, what, {c: (units.get(c, "none"),
                                       text.word(PAGE, f"tables.{what}.columns.{c}"))
                                   for c in frame.columns})  # fmt: skip


def _status_table(status: pl.DataFrame) -> ResultTable:
    """The flagged runs, as the page lists them."""
    flagged = status.filter(pl.col("qc_status") != "pass")
    rows = tuple((r, s, q, reason) for r, s, q, reason, _ in flagged.rows())
    return ResultTable(text.fill(PAGE, "status.title", flagged=flagged.height, runs=status.height),
                       tuple(text.words(PAGE, "status.columns")), rows)  # fmt: skip


def _confounding_table(tests: pl.DataFrame) -> ResultTable:
    """The confounding tests with the smallest adjusted p; the CSV holds every one."""
    tested = tests.filter(pl.col("p_adj").is_not_null()).sort("p_adj").head(SHOWN_TESTS)
    rows = tuple(
        (r["contrasts"], r["metric"], r["test"], f"{r['effect']:.2f}", f"{r['p_adj']:.2g}",
         text.word(PAGE, "confounding.tied") if r["read_length_tied"] else "")
        for r in tested.iter_rows(named=True)
    )  # fmt: skip
    total = tests["p_adj"].is_not_null().sum()
    title = text.fill(PAGE, "confounding.title", shown=len(rows), tests=total)
    columns = tuple(text.words(PAGE, "confounding.columns"))
    return ResultTable(title, columns, rows, numeric=(3, 4))


def _tally(status: pl.DataFrame, tests: pl.DataFrame) -> dict[str, int]:
    """Runs, flagged runs (excluded ones among them), tests run, and tests below ``ALPHA``."""
    return {
        "runs": status.height,
        "flagged": status.filter(pl.col("qc_status") != "pass").height,
        "excluded": status.filter(pl.col("qc_status") == "exclude").height,
        "tests": int(tests["p_adj"].is_not_null().sum()),
        "confounded": tests.filter(pl.col("p_adj") < ALPHA).height,
    }


def _keys(stage: Stage, tally: Mapping[str, int]) -> tuple[KeyNumber, ...]:
    """The numbers at the top of Results."""
    return (
        KeyNumber(text.word(PAGE, "keys.samples"), f"{stage.values.height:,}"),
        KeyNumber(text.word(PAGE, "keys.flagged"), f"{tally['flagged']:,}",
                  text.fill(PAGE, "keys.flagged_detail", **tally)),
        KeyNumber(text.word(PAGE, "keys.metrics"), f"{stage.bands['metric'].n_unique():,}"),
        KeyNumber(text.word(PAGE, "keys.confounded"), f"{tally['confounded']:,}",
                  text.fill(PAGE, "keys.confounded_detail", alpha=ALPHA)),
    )  # fmt: skip


def _sections(stage: Stage, registry: reg.Registry, tally: Mapping[str, int]) -> dict[str, Any]:
    """The page's fixed sections, those that state this run's numbers filled from it."""
    filled = {
        "input": text.fill(PAGE, "input", samples=stage.values.height, runs=tally["runs"],
                           registry=registry.path.name),
        "method": text.fill(PAGE, "method", k=f"{registry.mad_multiplier:g}"),
        "interpretation": text.fill(PAGE, "interpretation", alpha=ALPHA, **tally),
        "pitfalls": lines("pitfalls", minimum=registry.min_samples),
    }  # fmt: skip
    return text.skeleton(PAGE) | filled


def pipeline_groups(
    writer: Writer, found: tuple[pl.DataFrame | None, fastqc.Plots | None, Reads | None],
    stage: Stage, sheet: pl.DataFrame, registry: reg.Registry, options: Options,
) -> list[Group]:  # fmt: skip
    """The groups of the pipeline's reports and of the reads, where the run read them.

    ``found`` is the fates, the pipeline's FastQC plots and what the reads showed, each None
    where the run did not read it.
    """
    fate, plots, seen = found
    groups = []
    if fate is not None:
        groups.append(fate_figures(writer, fate, stage.values, sheet, registry, options))
    if plots is not None:
        units = sheet.unique("sample", keep="first", maintain_order=True)
        items = fastqc_figures(writer, plots, units, registry, options)
        groups.append(Group(text.word(PAGE, "groups.fastqc"), tuple(items)))
    if seen is not None:
        groups.append(read_figures(writer, seen, sheet, registry, options))
        run_units = sheet.select("run_accession", *options.strata, checkpoint.LENGTH_GROUP)
        run_strips = strips(writer, stage, run_units, registry, options, "run_accession")
        if run_strips:
            groups.append(Group(text.word(PAGE, "groups.run_metrics"), tuple(run_strips),
                                tabs=len(run_strips) > 1))  # fmt: skip
    return groups


def run(
    settings: runtime.Settings, cards: tuple[ParameterCard, ...], command: str, dataset: str,
    options: Options,
) -> list[Path]:  # fmt: skip
    """Judge the dataset's samples and write the M2 tables and page.

    Returns:
        The page, then every table.

    Raises:
        DatasetError: If the dataset file does not match the slug, or gather has not run.
        RegistryError: If the registry cannot be found or read.
    """
    spec = load(options.dataset_file)
    if spec.slug != dataset:
        raise DatasetError(f"{options.dataset_file}: slug {spec.slug!r}, not {dataset!r}")
    registry = load_registry(options)
    inputs = gathered_tables(settings.root, dataset)
    gathered = read_gathered(inputs)
    sheet = checkpoint.length_groups(gathered["samplesheet"], bool(spec.batch_bins))
    if missing := [c for c in {*options.strata, *options.peers} if c not in sheet.columns]:
        raise DatasetError(f"--strata and --peers name columns the sheet lacks: {missing}")
    measured, notes = counts_stage(settings, spec, gathered, options)
    samples = sheet.unique("sample", keep="first", maintain_order=True).select("sample")
    values = measured.values if measured else samples
    sexes = measured.sexes if measured else None
    inputs |= align_tables(settings.root, dataset)
    aligned = align_stage({w: inputs[w] for w in ALIGN_TABLES if w in inputs}, registry)
    if aligned is not None:
        values = values.join(aligned, on="sample", how="left")
    day = date.today()
    releases = {"registry": str(registry.path), **{f"input:{w}": p.name for w, p in inputs.items()}}
    writer = Writer(settings.root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}", dataset,
                    day, releases)  # fmt: skip
    fate = plots = None
    if options.rnaseq_dir:
        found, fate, plots = rnaseq_stage(writer, options.rnaseq_dir, values["sample"].to_list(),
                                          registry)  # fmt: skip
        values = values.join(found, on="sample", how="left")
    seen = (
        reads_stage(writer, settings, dataset, sheet, registry, options) if options.reads else None
    )
    stage = judge(values, sheet, seen, registry, options, notes)
    status = checkpoint.status(sheet, stage.flags, registry, stage.run_flags)
    tests = confounding(stage, sheet, spec, registry)
    frames = {"bands": stage.bands, "flags": stage.flags, "status": status, "confounding": tests}
    frames |= {"sex_check": sexes} if sexes is not None else {}
    frames |= {"fates": fate} if fate is not None else {}
    if seen is not None and stage.run_flags is not None:
        frames |= {"read_flags": stage.run_flags.rename({"sample": "run_accession"}),
                   "skipped": seen.skipped, "grid": seen.grid}  # fmt: skip
        writer.table(seen.values, "read_metrics", columns_of(seen.values, registry))
    _write_tables(writer, stage.values, frames, registry)
    results: list[Result] = [*stage.notes, _status_table(status), _confounding_table(tests)]
    results += overview_groups(writer, stage, measured, (sheet, status), registry, options)
    results += pipeline_groups(writer, (fate, plots, seen), stage, sheet, registry, options)
    usage, budget = runtime.measure(settings.started), runtime.budget(MODULE)
    tally = _tally(status, tests)
    page = ModulePage(
        title=text.fill(PAGE, "title", dataset=spec.name), **_sections(stage, registry, tally),
        parameters=cards, results=tuple(results), key_numbers=_keys(stage, tally),
        reproduce=Reproduction(
            command=command, wall_seconds=usage.wall_seconds, peak_rss_self=usage.peak_rss_self,
            peak_rss_children=usage.peak_rss_children,
            profile=f"{settings.profile}, {settings.threads} threads",
            budget=text.fill(PAGE, "budget", minutes=budget.minutes, ram_gb=budget.ram_gb),
            versions=runtime.versions(), inputs=tuple(p.name for p in inputs.values()),
        ),
    )  # fmt: skip
    out = writer.out / f"{stem(MODULE, f'{dataset}_overview', 1, day=day)}.html"
    return [write(page, out), *writer.files]
