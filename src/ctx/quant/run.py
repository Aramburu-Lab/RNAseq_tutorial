"""M4 (`ctx quant`): what Salmon did with every sample's reads, and how they became gene counts.

Order of operations: every input problem is found and listed before anything runs; MultiQC
reads Salmon's logs; the port of tximport summarises every ``quant.sf`` when copy Q1 is there;
the checks of ``ctx.quant.checks`` must all hold before any table, figure or page is written.
Reads only the nf-core/rnaseq outdir and copy Q1; writes only under ``outputs/04_quant/``.
"""

from dataclasses import dataclass
from datetime import date
from pathlib import Path

import polars as pl

from ctx import runtime, tools
from ctx.dataset import Contrast, Dataset, DatasetError, load
from ctx.de.run import Writer, newest
from ctx.gather import rnaseq
from ctx.io.naming import MODULES, stem
from ctx.process.run import annotation
from ctx.qc.checkpoint import gene_types
from ctx.quant import CountsFrom, QuantError, checks, outdir, tables, tximport
from ctx.quant import page as pages
from ctx.report.builder import ParameterCard, write

MODULE = "quant"


@dataclass(frozen=True)
class Options:
    """What `ctx quant` was asked to do, beside the dataset.

    Attributes:
        dataset_file: The dataset file; its ``annotation_gtf`` gives gene types, its
            ``[de] featured`` contrast the fold-change figure.
        rnaseq_dir: The nf-core/rnaseq outdir ``ctx gather`` read.
        salmon_dir: Copy Q1, the salmon folder with every ``quant.sf``; None is
            ``<root>/data/raw/<dataset>_salmon`` where it exists, else no transcript figure.
        counts_from: Which gene counts ctx writes from its port of tximport.
    """

    dataset_file: Path
    rnaseq_dir: Path
    salmon_dir: Path | None = None
    counts_from: CountsFrom = CountsFrom.LENGTH_SCALED_TPM


@dataclass(frozen=True)
class Inputs:
    """Every input, checked.

    Attributes:
        samples: The sheet's samples, in sheet order, one row each.
        genes: The annotation's gene records (``checkpoint.gene_types``).
        counts: The length-scaled matrix as ``ctx gather`` reads it, integer.
        salmon_dir: Copy Q1, or None when it is not there.
        contrast: The featured contrast, or None.
    """

    spec: Dataset
    sheet_path: Path
    samples: pl.DataFrame
    rnaseq_dir: Path
    gtf: Path
    genes: pl.DataFrame
    counts: rnaseq.Counts
    salmon_dir: Path | None
    contrast: Contrast | None


@dataclass(frozen=True)
class Measured:
    """What the run measured before any check: Salmon's table, fates, sums, the port."""

    salmon: pl.DataFrame
    fates: pl.DataFrame
    sums: pl.DataFrame
    port: tximport.GeneTables | None


def _out(root: Path) -> Path:
    """The module's output folder."""
    return root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}"


def _spec(options: Options, dataset: str) -> tuple[Dataset | None, list[str]]:
    """The dataset file, read and checked, or why not."""
    if not options.dataset_file.is_file():
        return None, [f"--dataset-file {options.dataset_file} is not a file"]
    try:
        spec = load(options.dataset_file)
    except DatasetError as error:
        return None, [str(error)]
    if spec.slug != dataset:
        return spec, [f"{options.dataset_file}: slug {spec.slug!r}, not {dataset!r}"]
    return spec, []


def _sheet(root: Path, dataset: str) -> tuple[Path | None, pl.DataFrame | None, list[str]]:
    """The newest gather sheet and one row per sample, or why not."""
    try:
        path = newest(root / "outputs" / f"{dict(MODULES)['gather']:02d}_gather", dataset,
                      "samplesheet")  # fmt: skip
    except DatasetError as error:
        return None, None, [str(error)]
    sheet = pl.read_parquet(path).unique("sample", keep="first", maintain_order=True)
    return path, sheet, []


def _counts(folder: Path, samples: set[str]) -> tuple[rnaseq.Counts | None, list[str]]:
    """The rnaseq outdir's length-scaled matrix, checked as ``ctx gather`` checks it."""
    if not folder.is_dir():
        return None, [f"--rnaseq-dir {folder} is not a folder"]
    try:
        return rnaseq.read(folder, samples), []
    except rnaseq.RnaseqError as error:
        return None, [str(error)]


def salmon_dir(root: Path, dataset: str, given: Path | None) -> Path | None:
    """Copy Q1: the folder given, else the default where it exists, else None."""
    if given is not None:
        return given
    default = root / "data" / "raw" / f"{dataset}_salmon"
    return default if default.is_dir() else None


def q1_problems(folder: Path | None, samples: list[str]) -> list[str]:
    """Everything copy Q1 must hold that is missing: each sample's quant.sf, the tables."""
    if folder is None:
        return []
    if not folder.is_dir():
        return [f"--salmon-dir {folder} is not a folder"]
    problems = outdir.quant_problems(folder, samples)
    for what in outdir.GENE_TABLES:
        problems += outdir.table_problems(folder / outdir.MERGED.format(what), outdir.GENE_KEYS,
                                          samples)  # fmt: skip
    for what in outdir.TRANSCRIPT_TABLES:
        problems += outdir.table_problems(folder / outdir.MERGED.format(what),
                                          outdir.TRANSCRIPT_KEYS, samples)  # fmt: skip
    try:
        outdir.tx2gene(folder)
    except QuantError as error:
        problems.append(str(error))
    return problems


def _contrast(spec: Dataset) -> tuple[Contrast | None, list[str]]:
    """The dataset's featured DE contrast, or None when it names none."""
    if spec.de is None:
        return None, []
    found = [c for c in spec.contrasts if c.name == spec.de.featured]
    if not found:
        return None, [f"{spec.path.name}: [de] featured {spec.de.featured!r} is no contrast"]
    return found[0], []


def resolve(settings: runtime.Settings, dataset: str, options: Options) -> Inputs:
    """Every input checked, every problem listed together before anything runs.

    Raises:
        QuantError: Listing every problem found.
    """
    root = settings.root
    spec, problems = _spec(options, dataset)
    sheet_path, sheet, found = _sheet(root, dataset)
    problems += found
    gtf, found = annotation(root, spec) if spec else (None, [])
    problems += found
    folder = options.rnaseq_dir
    samples = sheet["sample"].to_list() if sheet is not None else []
    counts, found = _counts(folder, set(samples)) if sheet is not None else (None, [])
    problems += found
    q1 = salmon_dir(root, dataset, options.salmon_dir)
    problems += q1_problems(q1, samples)
    contrast, found = _contrast(spec) if spec else (None, [])
    problems += found
    try:
        tools.require("multiqc")
    except tools.ToolError as error:
        problems.append(str(error))
    if problems or spec is None or sheet is None or sheet_path is None or gtf is None \
            or counts is None:  # fmt: skip
        raise QuantError("ctx quant cannot start:\n  " + "\n  ".join(problems))
    return Inputs(spec, sheet_path, sheet, folder, gtf, gene_types(gtf), counts, q1, contrast)


def measure(inputs: Inputs, out: Path) -> Measured:
    """MultiQC on Salmon's logs, the fates and sums, and the port on every quant.sf."""
    samples = inputs.samples["sample"].to_list()
    run = outdir.salmon_run(inputs.rnaseq_dir, out / "multiqc" / inputs.spec.slug, samples)
    salmon = outdir.salmon_table(run)
    fates = tables.fates(salmon)
    scaled = outdir.column_sums(inputs.rnaseq_dir / "salmon", "gene_counts_length_scaled", samples)
    port = None
    if inputs.salmon_dir is not None:
        port = tximport.summarize(outdir.quants(inputs.salmon_dir, samples),
                                  outdir.tx2gene(inputs.salmon_dir))  # fmt: skip
    return Measured(salmon, fates, tables.sums(fates, scaled), port)


def verify(inputs: Inputs, measured: Measured) -> None:
    """Every check of ``ctx.quant.checks``, before anything is written.

    Raises:
        QuantError: Listing every failed check.
    """
    targets = tables.library(measured.salmon).select("sample", "targets")
    problems = checks.fate_problems(measured.fates) + checks.sum_problems(measured.sums, targets)
    problems += checks.ceiling_problems(inputs.counts.matrix["gene_id"].to_list(),
                                        inputs.genes["gene_id"].to_list())  # fmt: skip
    if measured.port is not None and inputs.salmon_dir is not None:
        problems += checks.tximport_problems(measured.port, inputs.salmon_dir)
    checks.check(problems)


def releases(inputs: Inputs) -> dict[str, str]:
    """What every README records: the outdir's provenance, the annotation and the inputs."""
    return {
        **inputs.counts.provenance,
        "annotation": inputs.gtf.name, "tool:multiqc": tools.TOOLS["multiqc"].version,
        "input:samplesheet": inputs.sheet_path.name,
        "input:salmon_dir": str(inputs.salmon_dir) if inputs.salmon_dir else "none",
    }  # fmt: skip


def run(settings: runtime.Settings, cards: tuple[ParameterCard, ...], command: str,
        dataset: str, options: Options) -> list[Path]:  # fmt: skip
    """Read Salmon's outputs, check them, and write the M4 page, tables and figures.

    Returns:
        The page, then every table and figure file.

    Raises:
        QuantError: Listing every input problem, when nothing has run; or every failed
            check, when only MultiQC's folder is written.
    """
    inputs = resolve(settings, dataset, options)
    out, day = _out(settings.root), date.today()
    measured = measure(inputs, out)
    verify(inputs, measured)
    writer = Writer(out, dataset, day, releases(inputs))
    page = pages.build(pages.Context(writer, inputs, options, measured, cards, command,
                                     settings))  # fmt: skip
    path = out / f"{stem(MODULE, f'{dataset}_overview', 1, day=day)}.html"
    return [write(page, path), *writer.files]
