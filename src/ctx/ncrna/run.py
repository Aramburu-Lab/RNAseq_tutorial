"""``ctx ncrna``: the M8 tables, figures and page of a dataset's Step-1 half.

Every input is checked, and every problem listed together, before anything is computed:
the dataset file's ``[ncrna]``, gather's newest sample sheet and counts, M5's newest results
of the featured contrast, and the annotation, its SHA-256 the one the dataset file declares.
``page`` computes the figures and words the page.
"""

from dataclasses import dataclass
from datetime import date
from pathlib import Path

import numpy as np
import polars as pl

from ctx import runtime
from ctx.dataset import Contrast, Dataset, DatasetError, load
from ctx.de.deseq import select_samples
from ctx.de.filters import MIN_COUNT
from ctx.de.run import Writer, newest
from ctx.io.naming import MODULES
from ctx.ncrna import NcrnaError, biotypes, page
from ctx.process.run import annotation
from ctx.qc.checkpoint import gene_types
from ctx.report.builder import ParameterCard

MODULE = "ncrna"


@dataclass(frozen=True)
class Options:
    """Every option of ``ctx ncrna`` but the dataset.

    Attributes:
        dataset_file: The dataset file, with ``[ncrna]`` and the annotation resource.
        biotypes: The biotype groups drawn and tested, each its own BH stratum.
        min_count: The filter's reads at the median library size, for the filter figure.
    """

    dataset_file: Path
    biotypes: tuple[str, ...] = biotypes.TESTED
    min_count: int = MIN_COUNT


@dataclass(frozen=True)
class Inputs:
    """What the page is computed from, every table read once.

    Attributes:
        spec: The dataset file.
        contrast: The featured contrast.
        sheet: Gather's sample sheet, one row per sample.
        used: The contrast's samples, as M5 chose them, with ``group``.
        genes: The annotation's record of each count row, in count order, with ``biotype``.
        matrix: The counts, genes by the sheet's samples.
        results: M5's results of the featured contrast.
        gtf: The annotation file.
        read: The input files, for the READMEs.
    """

    spec: Dataset
    contrast: Contrast
    sheet: pl.DataFrame
    used: pl.DataFrame
    genes: pl.DataFrame
    matrix: np.ndarray
    results: pl.DataFrame
    gtf: Path
    read: tuple[Path, ...]


def _newest(folder: Path, dataset: str, what: str, problems: list[str]) -> Path | None:
    """The newest table of a kind, or None with the problem recorded."""
    try:
        return newest(folder, dataset, what)
    except DatasetError as error:
        problems.append(str(error))
        return None


def _spec(options: Options, dataset: str, problems: list[str]) -> tuple[Dataset, Contrast] | None:
    """The dataset file and its featured contrast, or None with the problems recorded."""
    try:
        spec = load(options.dataset_file)
    except DatasetError as error:
        problems.append(str(error))
        return None
    if spec.slug != dataset:
        problems.append(f"{options.dataset_file} is for {spec.slug!r}, not {dataset!r}")
    if spec.ncrna is None:
        problems.append(f"{options.dataset_file}: M8 needs [ncrna]")
        return None
    contrast = next(c for c in spec.contrasts if c.name == spec.ncrna.featured)
    return spec, contrast


def resolve(root: Path, dataset: str, options: Options) -> Inputs:
    """Every input read and checked.

    Raises:
        NcrnaError: Listing every problem found, before anything is computed.
    """
    problems: list[str] = []
    found = _spec(options, dataset, problems)
    gtf, missing = annotation(root, found[0]) if found else (None, [])
    problems += missing
    gathered = root / "outputs" / f"{dict(MODULES)['gather']:02d}_gather"
    de_dir = root / "outputs" / f"{dict(MODULES)['de']:02d}_de"
    sheet_path = _newest(gathered, dataset, "samplesheet", problems)
    counts_path = _newest(gathered, dataset, "counts", problems)
    results_path = (_newest(de_dir, dataset, f"{found[1].name}_results", problems)
                    if found else None)  # fmt: skip
    if problems or not found or gtf is None or sheet_path is None or counts_path is None \
            or results_path is None:  # fmt: skip
        raise NcrnaError("ctx ncrna cannot start:\n  " + "\n  ".join(problems))
    spec, contrast = found
    sheet = pl.read_parquet(sheet_path).unique("sample", keep="first", maintain_order=True)
    counts = pl.read_parquet(counts_path)
    if absent := [s for s in sheet["sample"] if s not in counts.columns]:
        raise NcrnaError(f"{counts_path.name} lacks {len(absent)} samples of the sheet, "
                         f"for example {absent[:3]}")  # fmt: skip
    used, _ = select_samples(sheet, contrast, spec.order)
    genes = biotypes.annotate(counts["gene_id"], gene_types(gtf))
    matrix = counts.select(sheet["sample"].to_list()).to_numpy()
    return Inputs(spec, contrast, sheet, used, genes, matrix, pl.read_parquet(results_path),
                  gtf, (sheet_path, counts_path, results_path, gtf))  # fmt: skip


def run(settings: runtime.Settings, cards: tuple[ParameterCard, ...], command: str,
        dataset: str, options: Options) -> list[Path]:  # fmt: skip
    """Build the M8 tables, figures and page.

    Args:
        settings: The resolved runtime.
        cards: Parameter cards for every option of this command.
        command: The command line, as the reader should retype it.
        dataset: The dataset's slug.
        options: The command's options.

    Returns:
        The page, then every table and figure file.

    Raises:
        NcrnaError: If an input is missing or does not fit the others.
    """
    chosen = biotypes.parse(options.biotypes)
    inputs = resolve(settings.root, dataset, options)
    out = settings.root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}"
    releases = {"annotation": inputs.gtf.name, "dataset file": str(inputs.spec.path),
                **{f"input:{p.name}": p.name for p in inputs.read[:3]}}  # fmt: skip
    writer = Writer(out, dataset, date.today(), releases)
    context = page.Context(inputs=inputs, biotypes=chosen, min_count=options.min_count,
                           settings=settings, command=command, cards=cards)  # fmt: skip
    return page.build(writer, context)
