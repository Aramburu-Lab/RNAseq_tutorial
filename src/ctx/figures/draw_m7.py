"""M7's figures redrawn from the tables ``ctx proteo`` wrote, by M7's own block functions.

The run is read back as the page held it when it drew (``ctx.proteo.figure_tables``): its
options, its genes, and what each figure reads besides. The class test and the fold-change
concordance are M7's own metric functions run again on the genes read back, its bootstrap
seeded as the run's. Each drawer calls the page's block function that drew the figure, with a
writer that writes nothing (``Recorder``), so the figure and its words are M7's own; a block
that draws several figures draws them all, and the one asked for is kept.
"""

import time
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ctx import runtime
from ctx.dataset import Dataset
from ctx.figures.kind import OneContrast
from ctx.figures.kinds_m7 import OneGene
from ctx.figures.recorder import Recorder
from ctx.proteo import figure_tables, page
from ctx.proteo import run as m7
from ctx.report.builder import FigureText

FOLDER = "07_proteo"


def _context(inputs: Mapping[str, Path], spec: Dataset, measured: m7.Measured) -> page.Context:
    """The page's context: the dataset file, the run's options and what was read back.

    No block a drawer calls reads the runtime, the command or the cards, which the page alone
    prints.
    """
    options = figure_tables.options_from(inputs["proteo_options"], inputs["dataset_file"])
    settings = runtime.Settings("laptop", Path(), 1, None, time.time())
    return page.Context(spec, options, measured, settings, "", ())


def _missing(inputs: Mapping[str, Path], spec: Dataset, what: str) -> tuple[Any, FigureText]:
    """F1 and F2, drawn by ``page._missingness``."""
    ctx = _context(inputs, spec, figure_tables.measured_from(inputs))
    recorder = Recorder(FOLDER, spec.slug)
    page._missingness(recorder, ctx)
    return recorder.one(what)


def quantified(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """Genes with a protein value in each sample."""
    return _missing(inputs, spec, "quantified")


def coverage(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """The quantified share per protein, a curve per RNA abundance tertile."""
    return _missing(inputs, spec, "coverage")


def _correlation(inputs: Mapping[str, Path], spec: Dataset, what: str,
                 contrast: str | None = None) -> tuple[Any, FigureText]:  # fmt: skip
    """The correlation figures, drawn by ``page._across``; F4 of one contrast where named."""
    measured = figure_tables.measured_from(inputs)
    ctx = _context(inputs, spec, measured)
    assert spec.proteo is not None  # noqa: S101 - the kind is made only under [proteo]
    m7.correlation_metrics(measured, spec.proteo, ctx.options.method)
    if spec.proteo.design == "tumour_normal":
        m7.normal_metrics(measured, spec.proteo, ctx.options.method)
    measured.extra["folds"] = [] if contrast is None else [figure_tables.folds_from(
        inputs["contrast_folds"], inputs["contrast_summary"], spec, contrast)]  # fmt: skip
    recorder = Recorder(FOLDER, spec.slug)
    page._across(recorder, ctx)
    return recorder.one(what)


def rho_by_class(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """Per-gene rho, a cumulative curve per gene class."""
    return _correlation(inputs, spec, "rho_by_class")


def rho_by_class_normal(inputs: Mapping[str, Path], _: Any,
                        spec: Dataset) -> tuple[Any, FigureText]:  # fmt: skip
    """The normals' per-gene rho, a cumulative curve per gene class."""
    return _correlation(inputs, spec, "rho_by_class_normal")


def worked(inputs: Mapping[str, Path], p: OneGene, spec: Dataset) -> tuple[Any, FigureText]:
    """One worked gene's RNA against protein, sample by sample."""
    return _correlation(inputs, spec, f"worked_{p.gene.lower()}")


def paired_zscores(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """The drawn genes' RNA row above their protein row, z-scored."""
    return _correlation(inputs, spec, "paired_zscores")


def contrast(inputs: Mapping[str, Path], p: OneContrast, spec: Dataset) -> tuple[Any, FigureText]:
    """One contrast's RNA against protein fold change, at the build's call rule."""
    return _correlation(inputs, spec, f"contrast_{p.contrast}", p.contrast)


def _studies(inputs: Mapping[str, Path], spec: Dataset, what: str) -> tuple[Any, FigureText]:
    """The cross-study figures, drawn by ``page._cross``.

    Each protein sample's group is its contrast level, so the groups count both levels' n.
    """
    measured = figure_tables.measured_from(inputs)
    ctx = _context(inputs, spec, measured)
    assert spec.proteo is not None  # noqa: S101 - the kind is made only under [proteo]
    m7.fold_metrics(measured)
    declared = next(c for c in spec.contrasts if c.name == spec.proteo.contrast)
    measured.extra |= {"contrast": declared, "n_tested": measured.groups.count(declared.tested),
                       "n_reference": measured.groups.count(declared.reference)}  # fmt: skip
    recorder = Recorder(FOLDER, spec.slug)
    page._cross(recorder, ctx)
    return recorder.one(what)


def fold_change(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """RNA against protein fold change per gene, a density per class."""
    return _studies(inputs, spec, "fold_change")


def shift_by_class(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """Protein minus RNA fold change, a cumulative curve per class."""
    return _studies(inputs, spec, "shift_by_class")
