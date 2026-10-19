"""M5's figures redrawn from the tables ``ctx de`` wrote, by M5's own functions.

A fit is read back from two tables: the contrast's results, whose ``gene_name`` M5 joined on
from gather's genes, and the size factors, which list the samples the fit used with their
donor and group; the dispersion plot reads the featured fit's dispersions too. Nothing is
refitted. The PCA computes M5's one VST again from gather's counts, and the named-gene boxes
normalise gather's counts by the size factors, as M5 does.
"""

import os
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import polars as pl

from ctx import runtime
from ctx.dataset import Dataset
from ctx.de import figures as m5_figures
from ctx.de import run as m5
from ctx.de.deseq import Fit
from ctx.de.page import (
    box_text,
    calls_text,
    calls_text_at,
    dispersion_text,
    every,
    figure,
    lfc_text,
    ma_text,
    phist_text,
    volcano_text,
)
from ctx.de.stats import ALPHA
from ctx.figures.kind import NoFigureError, OneContrast
from ctx.figures.kinds_m5 import Lesson, OnePca, TwoContrasts
from ctx.gather import gdc
from ctx.report.builder import FigureText
from ctx.viz.box import BoxData


def read_fit(
    results: Path, factors: Path, spec: Dataset, contrast: str, dispersions: Path | None = None
) -> tuple[Fit, pl.DataFrame]:
    """One contrast's fit as M5 held it when it drew, and the genes' symbols.

    Returns:
        The fit (no dropped runs, which no figure reads; dispersions where their table is
        given), and ``gene_id`` with ``gene_name``, as the figure functions join it.
    """
    table = pl.read_parquet(results)
    rows = pl.read_parquet(factors).filter(pl.col("contrast") == contrast)
    declared = next(c for c in spec.contrasts if c.name == contrast)
    fitted = Fit(
        contrast=declared, results=table.drop("gene_name"),
        samples=rows.select("sample", "donor", "group"), dropped=(),
        dispersions=pl.read_parquet(dispersions) if dispersions else pl.DataFrame(),
        size_factors=rows.select("sample", "size_factor"),
    )  # fmt: skip
    return fitted, table.select("gene_id", "gene_name")


def _one(inputs: Mapping[str, Path], p: OneContrast, spec: Dataset) -> tuple[Fit, pl.DataFrame]:
    """The contrast's fit and its genes."""
    return read_fit(inputs["results"], inputs["size_factors"], spec, p.contrast)


def volcano(inputs: Mapping[str, Path], p: OneContrast, spec: Dataset) -> tuple[Any, FigureText]:
    """The volcano, every tested gene drawn, live calls in the page."""
    fitted, genes = _one(inputs, p, spec)
    said = every(spec.contrasts, spec.columns)[p.contrast]
    return m5.volcano(fitted, genes, said=said), volcano_text(said)


def ma(inputs: Mapping[str, Path], p: OneContrast, spec: Dataset) -> tuple[Any, FigureText]:
    """Mean against fold change, raw and shrunk where the fit has both."""
    fitted, genes = _one(inputs, p, spec)
    said = every(spec.contrasts, spec.columns)[p.contrast]
    shrunk = bool(fitted.results["log2FoldChange_shrunk"].is_not_null().any())
    return m5._ma(fitted, genes, said), ma_text(said, shrunk)


def phist(inputs: Mapping[str, Path], p: OneContrast, spec: Dataset) -> tuple[Any, FigureText]:
    """The p-value histogram with its flat level."""
    fitted, _ = _one(inputs, p, spec)
    said = every(spec.contrasts, spec.columns)[p.contrast]
    data, flat, low = m5.p_histogram(fitted, said)
    with_p = fitted.results["pvalue"].drop_nulls().len()
    first = int(data.frame["genes"][0])
    return data, phist_text(said, with_p, low, flat, first, 1 / m5.P_BINS)


def lfc(inputs: Mapping[str, Path], p: TwoContrasts, spec: Dataset) -> tuple[Any, FigureText]:
    """One contrast's raw fold changes against another's."""
    x, genes = read_fit(inputs["results_x"], inputs["size_factors"], spec, p.x)
    y, _ = read_fit(inputs["results_y"], inputs["size_factors"], spec, p.y)
    said = every(spec.contrasts, spec.columns)
    data, rho = m5_figures.lfc_scatter(x, y, genes, (said[p.x], said[p.y]))
    return data, lfc_text(said[p.x].label, said[p.y].label, rho)


def calls(inputs: Mapping[str, Path], p: Lesson, spec: Dataset) -> tuple[Any, FigureText]:
    """One lesson's up and down calls at the setting the parameters name."""
    assert spec.de is not None  # noqa: S101 - the kind's check refuses a file without [de]
    names = spec.de.lessons[p.lesson]
    fits = {n: read_fit(inputs[f"results_{n}"], inputs["size_factors"], spec, n)[0]
            for n in names}  # fmt: skip
    said = every(spec.contrasts, spec.columns)
    lesson = {n: said[n] for n in names}
    sizes = "; ".join(
        figure("calls", "size", label=w.label, sizes=m5._sizes(fits[n], w),
               design=fits[n].contrast.design)
        for n, w in lesson.items()
    )  # fmt: skip
    setting = m5_figures.Setting(p.basis, p.alpha, p.lfc_min)
    labels = {n: w.label for n, w in lesson.items()}
    data = m5_figures.calls_bars(p.lesson, lesson, fits, sizes, setting)
    if setting == m5_figures.BUILD:
        return data, calls_text(p.lesson, labels, fits)
    counts = {n: sum(m5_figures.counted(fits[n], setting)) for n in names}
    return data, calls_text_at(p.lesson, labels, counts)


def _threads() -> int:
    """The threads the runtime profile set for this process, as M5's fit takes them."""
    return int(os.environ.get(runtime.THREAD_VARIABLES[0], "1"))


def pca(inputs: Mapping[str, Path], p: OnePca, spec: Dataset) -> tuple[Any, FigureText]:
    """One PCA of the blind VST, M5's own, on the samples the dataset file selects."""
    assert spec.de is not None  # noqa: S101 - the kind's check refuses a file without [de]
    declared = next(q for q in spec.de.pca if q.name == p.pca)
    sheet = gdc.chosen(pl.read_parquet(inputs["samplesheet"]))
    counts = pl.read_parquet(inputs["counts"])
    data, _, text, _ = m5._pca_figure(declared, sheet, counts, _threads(), spec.columns)
    return data, text


def box(inputs: Mapping[str, Path], p: OneContrast, spec: Dataset) -> tuple[Any, FigureText]:
    """The named genes' normalised counts in one contrast, MWU in the caption.

    Raises:
        NoFigureError: When the annotation holds none of the named genes.
    """
    fitted, _ = _one(inputs, p, spec)
    genes = pl.read_parquet(inputs["genes"])
    symbols = spec.de.box_genes if spec.de else ()
    said = every(spec.contrasts, spec.columns)[p.contrast]
    data, _ = m5._box(fitted, pl.read_parquet(inputs["counts"]), genes, symbols, said)
    if data is None:
        raise NoFigureError(f"{p.contrast}: no named gene is in the annotation, so M5 draws none")
    assert isinstance(data, BoxData)  # noqa: S101 - _box returns one
    tests = data.caption.tests
    return data, box_text(said, len(tests), sum(t.p_adj < ALPHA for t in tests))


def dispersion(inputs: Mapping[str, Path], p: OneContrast, spec: Dataset) -> tuple[Any, FigureText]:
    """Gene-wise and final dispersions with the fitted trend, from the table M5 wrote."""
    fitted, _ = read_fit(inputs["results"], inputs["size_factors"], spec, p.contrast,
                         inputs["dispersions"])  # fmt: skip
    said = every(spec.contrasts, spec.columns)[p.contrast]
    return m5._dispersion(fitted, said), dispersion_text(said)
