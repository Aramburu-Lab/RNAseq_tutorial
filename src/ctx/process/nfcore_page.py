"""M3 from nf-core/rnaseq's alignment (`ctx process --rnaseq-dir`): check, tables, page (D-116).

Order of operations: every input problem is listed together before anything is written (the
dataset file, the registry, the sheet, the outdir and its tables, the annotation); then every
frame is read from nf-core's own files (``ctx.process.nfcore``); then the tables, the tables
``ctx figure`` redraws from, the figures and the page. Nothing is aligned, counted or run.
"""

from collections import Counter
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any

import polars as pl

from ctx import concordance, runtime
from ctx.dataset import Dataset, DatasetError, load
from ctx.de.results import newest
from ctx.de.run import Writer
from ctx.io.naming import MODULES, stem
from ctx.process import ProcessError, figure_tables, nfcore
from ctx.process import blocks as process_blocks
from ctx.process import nfcore_blocks as blocks
from ctx.process import tables as tab
from ctx.process.run import annotation, runs_of
from ctx.qc import figures as qc_figures
from ctx.qc import registry as reg
from ctx.qc.checkpoint import gene_types
from ctx.report import text
from ctx.report.builder import Group, KeyNumber, ModulePage, ParameterCard, Reproduction, write

PAGE = "process_nfcore"
MODULE = "process"
# Sheet columns the bars may filter by and the curves be coloured by, where they have two
# levels or more and no more than a figure has colours.
GROUPS = ("cell_type", "condition")
# The registry metric whose plot line is nf-core's stranded threshold (D-116).
STRAND_SHARE = "rseqc_pct"
# The one comparison: Salmon alone (a) against Salmon on STAR's alignment (b).
COMPARISON = "salmon"


@dataclass(frozen=True)
class Options:
    """What `ctx process --rnaseq-dir` was asked, beside the dataset.

    Attributes:
        dataset_file: The dataset file; its ``annotation_gtf`` resource, where it declares
            one, gives the gene types the counter figures colour by.
        rnaseq_dir: An nf-core/rnaseq 3.27.0 outdir aligned with STAR, as the cluster's copy
            holds it (``ctx.gather.rnaseq.ALIGNED_FILES``).
        registry: The QC registry; None resolves it.
    """

    dataset_file: Path
    rnaseq_dir: Path
    registry: Path | None = None


@dataclass(frozen=True)
class Inputs:
    """Every input, checked."""

    spec: Dataset
    sheet: pl.DataFrame
    sheet_path: Path
    outdir: nfcore.Outdir
    registry: reg.Registry
    genes: pl.DataFrame | None
    notes: tuple[str, ...]


def _out(root: Path) -> Path:
    """The module's output folder."""
    return root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}"


def _spec(options: Options, dataset: str) -> tuple[Dataset | None, list[str]]:
    """The dataset file, read, its slug the dataset's, or why not."""
    if not options.dataset_file.is_file():
        return None, [f"--dataset-file {options.dataset_file} is not a file"]
    try:
        spec = load(options.dataset_file)
    except DatasetError as error:
        return None, [str(error)]
    if spec.slug != dataset:
        return spec, [f"{options.dataset_file}: slug {spec.slug!r}, not {dataset!r}"]
    return spec, []


def _registry(options: Options) -> tuple[reg.Registry | None, list[str]]:
    """The QC registry, holding the strand share whose line the page draws, or why not."""
    try:
        registry = reg.load(reg.resolve(options.registry))
    except reg.RegistryError as error:
        return None, [str(error)]
    if STRAND_SHARE not in registry.metrics:
        return None, [f"the QC registry {registry.path} has no {STRAND_SHARE}; M3 draws its line"]
    return registry, []


def _sheet(root: Path, dataset: str) -> tuple[pl.DataFrame | None, Path | None, list[str]]:
    """Gather's newest sheet, through the doors ``process.run`` reads it through, or why not."""
    folder = root / "outputs" / f"{dict(MODULES)['gather']:02d}_gather"
    try:
        path = newest(folder, dataset, "samplesheet")
    except DatasetError as error:
        return None, None, [str(error)]
    sheet, refused = runs_of(path)
    return sheet, path, refused


def _genes(root: Path, spec: Dataset | None) -> tuple[pl.DataFrame | None, list[str], list[str]]:
    """The annotation's gene types where the dataset file declares one: (types, problems, notes).

    A declared annotation that is missing or another file is a problem; none declared is a note,
    and the counter figures draw every gene as one type.
    """
    if spec is None:
        return None, [], []
    if not any(r.name == "annotation_gtf" for r in spec.resources):
        return None, [], [text.word(PAGE, "notes.no_annotation")]
    path, problems = annotation(root, spec)
    return (gene_types(path) if path is not None else None), problems, []


def resolve(settings: runtime.Settings, dataset: str, options: Options) -> Inputs:
    """Every input checked, every problem listed together before anything is written.

    Raises:
        ProcessError: Listing every problem found.
    """
    spec, problems = _spec(options, dataset)
    registry, found = _registry(options)
    problems += found
    sheet, sheet_path, found = _sheet(settings.root, dataset)
    problems += found
    genes, found, notes = _genes(settings.root, spec)
    problems += found
    outdir = None
    if not options.rnaseq_dir.is_dir():
        problems.append(f"--rnaseq-dir {options.rnaseq_dir} is not a folder")
    elif sheet is not None:
        samples = sheet.unique("sample", keep="first", maintain_order=True)["sample"].to_list()
        try:
            outdir = nfcore.open_outdir(options.rnaseq_dir, samples)
        except ProcessError as error:
            problems.append(str(error))
    if (problems or spec is None or registry is None or sheet is None or sheet_path is None
            or outdir is None):  # fmt: skip
        raise ProcessError("ctx process cannot start:\n  " + "\n  ".join(problems))
    return Inputs(spec, sheet, sheet_path, outdir, registry, genes, tuple(notes))


def alignment_frame(outdir: nfcore.Outdir) -> pl.DataFrame:
    """One row per sample: Trim Galore's reads, STAR's outcomes and splices, flagstat's records.

    Splices per 100 uniquely mapped reads, the annotated share and primary mapped % are this
    page's arithmetic on STAR's and samtools' counts, which neither prints.

    Raises:
        ProcessError: If STAR took in more reads than Trim Galore wrote for a sample.
    """
    star, flags, trim = nfcore.star(outdir), nfcore.flagstat(outdir), nfcore.trimmed(outdir)
    frame = trim.join(star, on="sample", maintain_order="left").join(
        flags, on="sample", maintain_order="left")  # fmt: skip
    if (more := frame.filter(pl.col("input_reads") > pl.col("trim_out"))["sample"]).len():
        raise ProcessError(f"STAR took in more reads than Trim Galore wrote for {more.to_list()}")
    extra = pl.col("secondary_passed") + pl.col("supplementary_passed")
    records = pl.col("total_passed") - extra
    return frame.with_columns(
        splices_per_100_unique=pl.when(pl.col("unique") > 0).then(
            100 * pl.col("splices") / pl.col("unique")),
        splices_annotated_share=pl.when(pl.col("splices") > 0).then(
            pl.col("splices_annotated") / pl.col("splices")),
        records=pl.col("total_passed").cast(pl.Int64),
        mapped_records_pct=pl.col("mapped_passed_pct"),
        primary_mapped_pct=pl.when(records > 0).then(
            100 * (pl.col("mapped_passed") - extra) / records),
    ).drop(*nfcore.FLAGSTAT)  # fmt: skip


def outcomes(alignment: pl.DataFrame) -> dict[str, dict[str, int]]:
    """Each sample's reads per STAR outcome and its input reads, under ``"input"``."""
    return {r["sample"]: {**{k: int(r[k]) for k in nfcore.OUTCOMES}, "input": int(r["input_reads"])}
            for r in alignment.iter_rows(named=True)}  # fmt: skip


def comparison(outdir: nfcore.Outdir, genes: pl.DataFrame | None
               ) -> tuple[pl.DataFrame, tab.Agreements]:  # fmt: skip
    """Salmon alone against Salmon on STAR's alignment, per gene: the long form and agreement.

    Agreement is ``concordance.per_unit``, through ``tables.agreements``, as the laptop route's.
    """
    long = nfcore.matrices(outdir)
    named = long.select("gene_id", "gene_name").unique("gene_id", maintain_order=True)
    if genes is not None:
        named = named.join(genes.select("gene_id", "gene_type"), on="gene_id", how="left",
                           maintain_order="left")  # fmt: skip
    else:
        named = named.with_columns(gene_type=pl.lit(None, pl.String))
    plain = long.select(comparison=pl.lit(COMPARISON), sample="sample", gene_id="gene_id",
                        a="a", b="b")  # fmt: skip
    agreed = tab.agreements(plain, named)
    return plain.join(named, on="gene_id", how="left", maintain_order="left"), agreed


def frames(inputs: Inputs) -> blocks.Found:
    """Every frame of the run, read from nf-core's own files."""
    outdir = inputs.outdir
    alignment = alignment_frame(outdir)
    complexity, perfect = nfcore.complexity(outdir)
    long, agreed = comparison(outdir, inputs.genes)
    return blocks.Found(
        alignment=alignment, outcomes=outcomes(alignment), strand=nfcore.strand(outdir),
        distribution=nfcore.read_distribution(outdir), junctions=nfcore.junctions(outdir),
        coverage=nfcore.coverage(outdir), duplication=nfcore.duplication(outdir),
        complexity=complexity, perfect=perfect, biotypes=nfcore.biotypes(outdir),
        comparison=long, agreements=agreed, steps=nfcore.steps(outdir),
    )  # fmt: skip


def described(key: str, columns: list[str]) -> dict[str, tuple[str, str]]:
    """(unit, meaning) of each column from ``tables.<key>.units`` and ``tables.<key>.columns``."""
    units = text.raw(PAGE, f"tables.{key}.units")
    return {c: (units.get(c, "none"), text.word(PAGE, f"tables.{key}.columns.{c}"))
            for c in columns}  # fmt: skip


# Each table's file description and the frame it holds; a figure's stem is never one of these.
TABLES = {
    "alignment": "alignment", "strand_check": "strand", "distribution_tags": "distribution",
    "junction_counts": "junctions", "coverage_profile": "coverage",
    "duplication_curve": "duplication", "complexity_curve": "complexity",
    "complexity_perfect": "perfect", "biotype_reads": "biotypes", "steps": "steps",
}  # fmt: skip


def write_tables(writer: Writer, found: blocks.Found, shown: blocks.Shown,
                 registry: reg.Registry) -> None:  # fmt: skip
    """Every table as Parquet, CSV and README, then the tables ``ctx figure`` redraws from."""
    for what, attr in TABLES.items():
        frame = getattr(found, attr)
        writer.table(frame, what, described(what, frame.columns))
    named = {"gene_comparison": found.comparison, "concordance": found.agreements.concordance,
             "discordant": found.agreements.discordant}  # fmt: skip
    for what, frame in named.items():
        writer.table(frame, what, described(what, frame.columns))
    shared = {
        "figure_options": options_row(shown, registry),
        "outcomes": figure_tables.outcome_rows(found.outcomes),
        "sums": found.agreements.sums, "groups": shown.groups,
    }  # fmt: skip
    for what, frame in shared.items():
        writer.table(frame, what, described(what, frame.columns))


def options_row(shown: blocks.Shown, registry: reg.Registry) -> pl.DataFrame:
    """What the page states beside the tables, in the laptop route's columns and its own.

    The counter figures are the laptop route's (``process.blocks``), whose redraw reads the
    aligner, featured sample, genome and registry from here; nf-core aligned with STAR on the
    whole genome and featured no sample.
    """
    return pl.DataFrame({
        "aligner": ["star"], "featured": [""], "genome_name": [""], "genome_genes": [0],
        "partial": [False], "registry": [str(registry.path)], "colour": [shown.colour],
        "route": ["nfcore"],
    }, schema_overrides={"colour": pl.String})  # fmt: skip


def groups_of(sheet: pl.DataFrame, samples: list[str]) -> pl.DataFrame:
    """The samples and the sheet columns worth filtering their bars by (``page.groups_of``)."""
    from ctx.process.page import groups_of as laptop_groups  # noqa: PLC0415 - one rule, one place

    return laptop_groups(sheet, samples)


def results(writer: Writer, found: blocks.Found, shown: blocks.Shown,
            notes: tuple[str, ...]) -> tuple[Any, ...]:  # fmt: skip
    """The Results section in order: notes, the alignment, strandedness, nf-core's QC, counts."""
    scatters, plots = process_blocks.agreement_figures(writer, found, shown)  # type: ignore[arg-type]
    agree = found.agreements.concordance.filter(
        (pl.col("comparison") == COMPARISON) & (pl.col("unit") == concordance.ALL)
    ).row(0, named=True)
    return (
        *notes,
        Group(text.word(PAGE, "groups.align"), (
            blocks.mapping_figure(writer, found, shown), blocks.alignment_table(found),
            text.word(PAGE, "tables.alignment.page.note"))),
        Group(text.word(PAGE, "groups.strand"), (
            blocks.strand_figure(writer, found, shown),
            blocks.strand_check_figure(writer, found, shown), blocks.strand_table(found))),
        Group(text.word(PAGE, "groups.qc"), (
            blocks.distribution_figure(writer, found, shown),
            blocks.junction_figure(writer, found, shown),
            blocks.coverage_figure(writer, found, shown),
            blocks.duplication_figure(writer, found, shown),
            blocks.complexity_figure(writer, found, shown),
            blocks.biotype_figure(writer, found, shown)), tabs=True),
        Group(text.word(PAGE, "groups.counters"), (*scatters, *plots),
              lead=text.fill(PAGE, "groups.counters_lead", n=f"{agree['n_genes']:,}",
                             same=f"{agree['n_identical']:,}"), tabs=True),
        blocks.concordance_table(found),
        blocks.steps_table(found),
    )  # fmt: skip


def key_numbers(found: blocks.Found) -> tuple[KeyNumber, ...]:
    """Samples, median uniquely mapped, nf-core's commonest strand call, genes that differ."""
    alignment = found.alignment
    low = alignment.sort("unique_pct").row(0, named=True)
    median: Any = alignment["unique_pct"].median()
    reads: Any = alignment["input_reads"].median()
    call, n_call = Counter(found.strand["rseqc_inferred"]).most_common(1)[0]
    agree = found.agreements.concordance.filter(pl.col("unit") == concordance.ALL).row(
        0, named=True)  # fmt: skip
    return (
        KeyNumber(text.word(PAGE, "keys.samples"), f"{alignment.height:,}",
                  text.fill(PAGE, "keys.samples_detail", reads=f"{reads:,.0f}")),
        KeyNumber(text.word(PAGE, "keys.unique"), f"{median:.1f} %",
                  text.fill(PAGE, "keys.unique_detail", low=f"{low['unique_pct']:.1f}",
                            sample=low["sample"])),
        KeyNumber(text.word(PAGE, "keys.strand"), call,
                  text.fill(PAGE, "keys.strand_detail", n_call=n_call, n=found.strand.height)),
        KeyNumber(text.word(PAGE, "keys.differ"),
                  f"{agree['n_genes'] - agree['n_identical']:,}",
                  text.fill(PAGE, "keys.differ_detail", n=f"{agree['n_genes']:,}")),
    )  # fmt: skip


def sections(inputs: Inputs, found: blocks.Found) -> dict[str, Any]:
    """The fixed sections, those that state this run's numbers filled from it."""
    unique = found.alignment["unique_pct"]
    low: Any = unique.min()
    high: Any = unique.max()
    params = inputs.outdir.params
    values = {"samples": found.alignment.height, "fasta": Path(str(params.get("fasta"))).name,
              "gtf": Path(str(params.get("gtf"))).name,
              "low": f"{float(low):.1f}", "high": f"{float(high):.1f}",
              "release": inputs.outdir.provenance.get("nf-core/rnaseq", "")}  # fmt: skip
    filled = {k: text.fill(PAGE, k, **values) for k in ("input", "method", "interpretation")}
    return text.skeleton(PAGE) | filled


def _versions(inputs: Inputs) -> dict[str, str]:
    """The versions of ctx, then nf-core's release, Nextflow and Salmon as its outdir records."""
    kept = {k: v for k, v in inputs.outdir.provenance.items() if not k.startswith(
        ("sha256:", "params "))}  # fmt: skip
    return {**runtime.versions(), **kept, "multiqc (nf-core)": inputs.outdir.mqc.version}


def build(settings: runtime.Settings, inputs: Inputs, writer: Writer,
          cards: tuple[ParameterCard, ...], command: str) -> ModulePage:  # fmt: skip
    """Write every table and figure, and return the page."""
    found = frames(inputs)
    samples = found.alignment["sample"].to_list()
    groups = groups_of(inputs.sheet, samples)
    colour = qc_figures.colour_column(groups, GROUPS)
    shown = blocks.Shown(groups, colour, inputs.registry.metrics[STRAND_SHARE])
    write_tables(writer, found, shown, inputs.registry)
    body = results(writer, found, shown, inputs.notes)
    body = (*body, blocks.reference_table(inputs.outdir.params))
    usage, budget = runtime.measure(settings.started), runtime.budget(MODULE)
    files = [inputs.sheet_path.name, *nfcore.MATRICES.values(),
             f"{inputs.outdir.mqc.data.relative_to(inputs.outdir.path)}/"]  # fmt: skip
    return ModulePage(
        title=text.fill(PAGE, "title", dataset=inputs.spec.name), **sections(inputs, found),
        parameters=cards, results=body, key_numbers=key_numbers(found),
        reproduce=Reproduction(
            command=command, wall_seconds=usage.wall_seconds,
            peak_rss_self=usage.peak_rss_self, peak_rss_children=usage.peak_rss_children,
            profile=f"{settings.profile}, {settings.threads} threads",
            budget=text.fill(PAGE, "budget", minutes=budget.minutes, ram_gb=budget.ram_gb),
            versions=_versions(inputs), inputs=tuple(files),
            notes=(text.word(PAGE, "reproduce_note"),),
        ),
    )  # fmt: skip


def releases(inputs: Inputs) -> dict[str, str]:
    """What every README records: nf-core's release and params, each file's SHA-256, inputs."""
    return {**inputs.outdir.provenance, "multiqc (nf-core)": inputs.outdir.mqc.version,
            "input:samplesheet": inputs.sheet_path.name,
            "input:rnaseq_dir": inputs.outdir.path.name}  # fmt: skip


def run(settings: runtime.Settings, cards: tuple[ParameterCard, ...], command: str,
        dataset: str, options: Options) -> list[Path]:  # fmt: skip
    """Read nf-core's alignment of the dataset and write the M3 page and its files.

    Returns:
        The page, then every table and figure file.

    Raises:
        ProcessError: Listing every input problem, when nothing is written; or a check of
            nf-core's tables that failed, before any table, figure or page is written.
    """
    inputs = resolve(settings, dataset, options)
    day = date.today()
    writer = Writer(_out(settings.root), dataset, day, releases(inputs))
    page = build(settings, inputs, writer, cards, command)
    out = writer.out / f"{stem(MODULE, f'{dataset}_overview', 1, day=day)}.html"
    return [write(page, out), *writer.files]
