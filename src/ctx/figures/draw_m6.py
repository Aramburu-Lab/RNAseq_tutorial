"""M6's figures redrawn from the tables M6 and M5 wrote, by M6's own functions.

Nothing is recomputed: the GSEA and ORA tables are the ones the page drew, so no permutation
reruns. A redrawn ORA figure carries no live panel, which needs the gene-set file and every
contrast's results; the page's other views stay on the module page.
"""

from collections.abc import Mapping
from pathlib import Path
from typing import Any

import polars as pl

from ctx.dataset import Dataset
from ctx.de.page import every
from ctx.de.stats import ALPHA
from ctx.figures.kinds import NoFigureError, OneContrast, OneSet
from ctx.function import run as m6
from ctx.function.page import gsea_text, nes_text, ora_text
from ctx.report.builder import FigureText


def gsea(inputs: Mapping[str, Path], p: OneContrast, spec: Dataset) -> tuple[Any, FigureText]:
    """The contrast's GSEA bars, the sets with the smallest FDR."""
    table = pl.read_parquet(inputs["gsea"])
    said = every(spec.contrasts, spec.columns)[p.contrast]
    ranked = m6._gsea_bars(said, table)
    drawn = (int((ranked.frame["fdr"] < ALPHA).sum()), ranked.frame.height)
    return ranked, gsea_text(said, int((table["fdr"] < ALPHA).sum()), table.height, drawn)


def ora(inputs: Mapping[str, Path], p: OneContrast, spec: Dataset) -> tuple[Any, FigureText]:
    """The contrast's ORA bars over its tested genes.

    Raises:
        NoFigureError: When no set passes, where M6 prints a sentence instead of the figure.
    """
    table = pl.read_parquet(inputs["ora"])
    background = m6._symbols(pl.read_parquet(inputs["results"]))
    said = every(spec.contrasts, spec.columns)[p.contrast]
    bars = m6._ora_bars(said, table, len(background))
    if bars is None:
        raise NoFigureError(
            f"{p.contrast}: no gene set passes FDR < {ALPHA}, so M6 draws no ORA bars"
        )
    passing = dict(table.filter(pl.col("fdr") < ALPHA)["direction"].value_counts().rows())
    return bars, ora_text(said, len(background), passing.get("up", 0), passing.get("down", 0))


def nes(inputs: Mapping[str, Path], p: OneSet, spec: Dataset) -> tuple[Any, FigureText]:
    """One set's NES in every contrast that tested it.

    Raises:
        NoFigureError: When no contrast tested the set.
    """
    tables = {c.name: pl.read_parquet(inputs[f"gsea_{c.name}"]) for c in spec.contrasts}
    tested = [n for n, t in tables.items() if p.gene_set in set(t["gene_set"])]
    if not tested:
        raise NoFigureError(f"{p.gene_set}: no contrast tested it, so M6 draws no NES bars")
    data = m6._featured_bars(p.gene_set, tables, every(spec.contrasts, spec.columns))
    return data, nes_text(p.gene_set, len(tested))
