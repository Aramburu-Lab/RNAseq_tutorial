"""``ctx clinical``: the M10 tables, figures and page of a GDC dataset, from gather's tables alone.

The patient table (``harmonise``), each patient's tumour through ``gdc.chosen()``; expression
from gather's counts, stored and filtered by ``score``; the panel the page scores live
(``panel``); ``page`` computes the analyses, writes the tables and figures, and words the page.
"""

from dataclasses import dataclass
from datetime import date
from pathlib import Path

import numpy as np
import polars as pl

from ctx import runtime
from ctx.clinical import ClinicalError, harmonise, page, panel, score
from ctx.dataset import ENDPOINTS, Dataset, DatasetError, load
from ctx.de.filters import log2cpm
from ctx.de.run import Writer, newest
from ctx.io.naming import MODULES
from ctx.report import genes, text
from ctx.report.builder import ParameterCard

MODULE = "clinical"
GATHER_TABLES = ("samplesheet", "counts", "genes", "cdr", "cases", "molecular_tests",
                 "treatments", "annotations")  # fmt: skip
COVARIATES = ("age", "stage")
CODING = "protein_coding"


@dataclass(frozen=True)
class Options:
    """Every option of ``ctx clinical`` but the dataset.

    Attributes:
        dataset_file: The dataset file, with ``[gdc]`` and ``[clinical]``.
        endpoint: OS, DSS, DFI or PFI; None takes the dataset file's.
        split: median or tertiles.
        horizon: Months up to which RMST is compared.
        covariates: Patient-table columns the adjusted Cox model adds.
        strata: Patient-table columns given their own baseline hazard.
        genes_from: A gene list file in place of the dataset file's build list.
        min_log2cpm: The expression filter's median log2 CPM.
    """

    dataset_file: Path
    endpoint: str | None = None
    split: str = "median"
    horizon: float = 60.0
    covariates: tuple[str, ...] = COVARIATES
    strata: tuple[str, ...] = ()
    genes_from: Path | None = None
    min_log2cpm: float = 1.0


def gathered(root: Path, dataset: str) -> dict[str, pl.DataFrame]:
    """Gather's newest table of each kind for the dataset.

    Raises:
        DatasetError: If one is missing; ``ctx gather --gdc-dir`` writes them all.
    """
    out = root / "outputs" / f"{dict(MODULES)['gather']:02d}_gather"
    return {name: pl.read_parquet(newest(out, dataset, name)) for name in GATHER_TABLES}


def expression(tables: dict[str, pl.DataFrame], patients: pl.DataFrame, min_log2cpm: float
               ) -> tuple[pl.DataFrame, np.ndarray, tuple[int, int]]:  # fmt: skip
    """The patients' tumours' log2 CPM, stored as int16, kept by the expression filter.

    Returns:
        The kept genes, their stored values (genes by patients, in table order), and how many
        the filter dropped for low median and for no variance.
    """
    sizes = dict(tables["samplesheet"].select("sample", "library_size").rows())
    samples = patients["sample"].to_list()
    matrix = tables["counts"].select(samples).to_numpy()
    values = score.stored(log2cpm(matrix, np.array([sizes[s] for s in samples])))
    kept, flat = score.expressed(values, min_log2cpm)
    dropped = (int((~kept & ~flat).sum()), int(flat.sum()))
    return tables["genes"].filter(pl.Series(kept)), values[kept], dropped


def build_list(spec: Dataset, genes_from: Path | None, table: genes.GeneTable
               ) -> tuple[list[int], str]:  # fmt: skip
    """The rows of the genes scored by default, and where the list came from.

    Raises:
        ClinicalError: If the list names a gene ambiguously, or none the expressed genes hold.
    """
    assert spec.clinical is not None  # noqa: S101 - run() refuses a dataset without [clinical]
    if genes_from is not None:
        found, source = genes.read_list(genes_from, table).genes, genes_from.name
    else:
        found = genes.parse("\n".join(spec.clinical.genes), table)
        source = text.word("clinical", "live.build")
    if found.ambiguous or not found.found:
        named = found.ambiguous or found.missing
        raise ClinicalError(f"{source}: no expressed gene named unambiguously: {named}")
    return [table.row[n] for n in found.found], source


def run(settings: runtime.Settings, cards: tuple[ParameterCard, ...], command: str,
        dataset: str, options: Options) -> list[Path]:  # fmt: skip
    """Build the M10 tables, figures and page.

    Args:
        settings: The resolved runtime.
        cards: Parameter cards for every option of this command.
        command: The command line, as the reader should retype it.
        dataset: The dataset's slug.
        options: The command's options.

    Returns:
        The page, then every table and figure file.

    Raises:
        DatasetError: If the dataset file names another slug or lacks [gdc] or [clinical].
        ClinicalError: If gather's tables cannot make the analysis, or an option names
            what the patient table lacks.
    """
    spec = load(options.dataset_file)
    if spec.slug != dataset:
        raise DatasetError(f"{options.dataset_file} is for {spec.slug!r}, not {dataset!r}")
    if spec.gdc is None or spec.clinical is None:
        raise DatasetError(f"{options.dataset_file}: M10 needs [gdc] and [clinical]")
    endpoint = options.endpoint or spec.clinical.endpoint
    if endpoint not in ENDPOINTS:
        raise ClinicalError(f"--endpoint {endpoint!r}: one of {ENDPOINTS}")
    if options.split not in score.SPLITS:
        raise ClinicalError(f"--split {options.split!r}: one of {score.SPLITS}")
    if bad := [c for c in options.covariates if c not in COVARIATES]:
        raise ClinicalError(f"--covariates {bad}: the model takes {COVARIATES}")
    if bad := [c for c in options.strata if c not in page.CATEGORICAL]:
        raise ClinicalError(f"--cox-strata {bad}: one of {page.CATEGORICAL}")
    tables = gathered(settings.root, dataset)
    clinical = harmonise.build(tables["samplesheet"], tables["cdr"], tables["cases"],
                               tables["molecular_tests"], tables["treatments"],
                               tables["annotations"], tumour=spec.gdc.tumour)  # fmt: skip
    gene_table, values, dropped = expression(tables, clinical.patients, options.min_log2cpm)
    lookup = genes.table(gene_table["gene_id"].to_list(), gene_table["gene_name"].to_list())
    chosen, source = build_list(spec, options.genes_from, lookup)
    coding = (gene_table["gene_type"] == CODING).to_numpy()
    carried = panel.choose(values, coding, np.array(chosen))
    out = settings.root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}"
    releases = {"GDC": spec.gdc.release, "dataset file": str(spec.path),
                "TCGA-CDR": "Liu et al. 2018 Cell, Table S1, as gather copied it"}  # fmt: skip
    writer = Writer(out, dataset, date.today(), releases)
    context = page.Context(
        spec=spec, options=options, endpoint=endpoint, clinical=clinical, values=values,
        gene_table=gene_table, carried=carried, chosen=chosen,
        source=source, filtered=dropped, settings=settings, command=command, cards=cards,
    )  # fmt: skip
    return page.build(writer, context)
