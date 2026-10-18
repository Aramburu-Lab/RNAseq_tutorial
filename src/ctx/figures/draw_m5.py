"""M5's figures redrawn from the tables ``ctx de`` wrote, by M5's own functions.

A fit is read back from two tables: the contrast's results, whose ``gene_name`` M5 joined on
from gather's genes, and the size factors, which list the samples the fit used with their
donor and group. Nothing is refitted; the dispersions are not written, so the dispersion plot
has no spec.
"""

from collections.abc import Mapping
from pathlib import Path
from typing import Any

import polars as pl

from ctx.dataset import Dataset
from ctx.de import figures as m5_figures
from ctx.de import run as m5
from ctx.de.deseq import Fit
from ctx.de.page import (
    calls_text,
    calls_text_at,
    every,
    figure,
    lfc_text,
    ma_text,
    phist_text,
    volcano_text,
)
from ctx.figures.kinds import Lesson, OneContrast, TwoContrasts
from ctx.report.builder import FigureText


def read_fit(
    results: Path, factors: Path, spec: Dataset, contrast: str
) -> tuple[Fit, pl.DataFrame]:
    """One contrast's fit as M5 held it when it drew, and the genes' symbols.

    Returns:
        The fit (no dispersions, no dropped runs: no figure with a spec reads them), and
        ``gene_id`` with ``gene_name``, as the figure functions join it.
    """
    table = pl.read_parquet(results)
    rows = pl.read_parquet(factors).filter(pl.col("contrast") == contrast)
    declared = next(c for c in spec.contrasts if c.name == contrast)
    fitted = Fit(
        contrast=declared, results=table.drop("gene_name"),
        samples=rows.select("sample", "donor", "group"), dropped=(),
        dispersions=pl.DataFrame(), size_factors=rows.select("sample", "size_factor"),
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
