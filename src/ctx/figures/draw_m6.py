"""M6's figures redrawn from the tables M6 and M5 wrote, by M6's own functions.

The GSEA, ORA and NES bars read the tables the page drew, so no permutation reruns. The
running score needs a set's positions along the ranking, which no table holds: it runs M6's own
GSEA again on M5's results and the gene-set file, seed and permutations as M6 ran them, since
the FDR pools every set's null. The background lesson's build view tests the up calls against
the tested genes and the gene-set union, which needs no annotation. A redrawn ORA figure
carries no live panel, which needs every contrast's results; the page's other views stay on
the module page.
"""

from collections.abc import Mapping
from pathlib import Path
from typing import Any

import polars as pl

from ctx.dataset import Dataset
from ctx.de.page import every
from ctx.de.stats import ALPHA
from ctx.figures.kind import NoFigureError, OneContrast
from ctx.figures.kinds_m6 import OneSet, SetInContrast
from ctx.function import background as m6_background
from ctx.function import lite
from ctx.function import run as m6
from ctx.function.page import gsea_text, nes_text, ora_text, running_text, trap_text
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


def running(inputs: Mapping[str, Path], p: SetInContrast, spec: Dataset) -> tuple[Any, FigureText]:
    """One set's running enrichment score, from M6's GSEA of the contrast run again.

    Raises:
        NoFigureError: When the contrast did not test the set, where M6 prints a sentence.
    """
    results = pl.read_parquet(inputs["results"])
    enriched = {r.name: r for r in m6._gsea(results, lite.read_gmt(inputs["gene_sets"]))}
    if (hit := enriched.get(p.gene_set)) is None:
        raise NoFigureError(f"{p.gene_set}: {p.contrast} did not test it, so M6 draws no curve")
    said = every(spec.contrasts, spec.columns)[p.contrast]
    return m6._running(said, hit), running_text(said, hit, len(hit.running))


def background(inputs: Mapping[str, Path], p: OneContrast, spec: Dataset) -> tuple[Any, FigureText]:
    """The background lesson as built: the up calls against the tested genes and the union."""
    results = pl.read_parquet(inputs["results"])
    sets = lite.read_gmt(inputs["gene_sets"])
    tested = m6._symbols(results)
    # The build's view: the first background, the first calls, the first sets drawn.
    wide, direction = m6_background.BACKGROUNDS[0], m6_background.DIRECTIONS[0]
    how = m6_background.SHOWN[0]
    against = m6_background.backgrounds(tested, sets, None)[wide]
    found = m6_background.tested_against(m6._symbols(results, direction), tested, against, sets,
                                         direction, wide)  # fmt: skip
    said = every(spec.contrasts, spec.columns)[p.contrast]
    calls = m6._tested(results, direction).height
    bars = m6_background.bars(said, calls, len(tested), found, how)
    more_less, shown = found.moved(how)[:2], found.moved(how)[2]
    return bars, trap_text(len(tested), found.size, more_less, shown, found.passing())
