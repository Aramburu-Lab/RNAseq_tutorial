"""The three drawn panels of a card: GTEx tissues, HPA immune cells, GTEx exon use.

Each returns a figure contract and the numbers its words are filled with; the page writes
them. Every panel is descriptive: n, median, IQR, tau, CV; external data get no test.
"""

import math
from typing import Any

import numpy as np
import polars as pl

from ctx.genecard import annotation, gtex, hpa
from ctx.genecard.cards import Card
from ctx.genecard.summary import tau
from ctx.report import text
from ctx.viz.bar import BarData
from ctx.viz.heatmap import HeatmapData, diverging_levels
from ctx.viz.style import DATA_COLOURS, PROJECTOR_UNSAFE, category_colours

PAGE = "genecard"
SEXES = ("female", "male")
# One colour per HPA lineage: every data colour but the yellow a projector washes out.
LINEAGE_COLOURS = tuple(c for c in DATA_COLOURS if c not in PROJECTOR_UNSAFE)


def columns(key: str) -> dict[str, tuple[str, str]]:
    """(unit, meaning) of every column of one figure's table, from the content file."""
    return {c: (u, m) for c, (u, m) in text.raw(PAGE, f"figures.{key}.columns").items()}


def shown(value: float | None) -> str:
    """A tau as a caption prints it; an undefined one as the page's word for none."""
    return text.word(PAGE, "none") if value is None else f"{value:.2f}"


def _log(value: float) -> float:
    """log2 of a TPM plus one, the axis every tissue panel draws on."""
    return math.log2(value + 1.0)


def tissues(card: Card) -> tuple[BarData, dict[str, Any]]:
    """GTEx median per tissue and sex, as bars, tissues by their all-donor median.

    Bars, not points over a range: a horizontal bar chart grows with its 68 tissue names,
    which a strip's fixed height could not hold. Quartiles and n are on hover and in the
    table.
    """
    assert card.tissues is not None  # noqa: S101 - the page draws this panel only with data
    spread = card.tissues
    every = spread.filter(pl.col("sex") == gtex.ALL).sort("median", descending=True)
    order = every["tissue"].to_list()
    rank = {t: i for i, t in enumerate(order)}
    bars = (
        spread.filter(pl.col("sex").is_in(SEXES))
        .with_columns(value=pl.col("median").map_elements(_log, return_dtype=pl.Float64),
                      slot=pl.col("tissue").replace_strict(rank, return_dtype=pl.Int64))
        .sort("slot", "sex").select("tissue", "sex", "value", "n", "median", "q1", "q3")
    )  # fmt: skip
    numbers = {"symbol": card.symbol, "release": gtex.RELEASE, "top": order[0],
               "top_median": every["median"][0], "bottom": order[-1],
               "bottom_median": every["median"][-1], "tau": shown(tau(every["median"].to_numpy())),
               "tissues": len(order), "samples": f"{int(every['n'].sum()):,}",
               "least": every["n"].min(), "most": every["n"].max()}  # fmt: skip
    data = BarData(frame=bars, columns=columns("tissues"), category="tissue", value="value",
                   category_label=text.word(PAGE, "figures.tissues.x"),
                   value_label=text.word(PAGE, "figures.tissues.y"), series="sex",
                   colours=category_colours(list(SEXES)),
                   caption=text.fill(PAGE, "figures.tissues.caption", **numbers),
                   horizontal=True, hover=("n", "median", "q1", "q3"))  # fmt: skip
    return data, numbers


def immune_tau(xml: hpa.Xml) -> float | None:
    """Our tau over HPA's immune cell types, the pooled PBMC total left out."""
    return tau([v for cell, _, v in xml.immune if cell != hpa.POOLED])


def immune(card: Card, release: str) -> tuple[BarData, dict[str, Any]]:
    """HPA immune cell nTPM, one bar per cell type, coloured by lineage; our tau beside it."""
    assert card.xml is not None  # noqa: S101 - the page draws this panel only with data
    rows = [(c, lineage, v) for c, lineage, v in card.xml.immune if c != hpa.POOLED]
    frame = pl.DataFrame(rows, schema={"cell_type": pl.String, "lineage": pl.String,
                                       "ntpm": pl.Float64}, orient="row")  # fmt: skip
    lineages = list(dict.fromkeys(frame["lineage"].to_list()))
    if len(lineages) > len(LINEAGE_COLOURS):
        raise ValueError(f"{len(lineages)} HPA lineages; {len(LINEAGE_COLOURS)} colours")
    specificity = immune_tau(card.xml)
    top = frame.sort("ntpm", descending=True).row(0)
    numbers = {"symbol": card.symbol, "release": release, "top": top[0], "top_ntpm": top[2],
               "tau": shown(specificity), "cells": frame.height}  # fmt: skip
    bars = BarData(frame=frame, columns=columns("immune"), category="cell_type", value="ntpm",
                   category_label=text.word(PAGE, "figures.immune.x"),
                   value_label=text.word(PAGE, "figures.immune.y"), series="lineage",
                   colours=dict(zip(lineages, LINEAGE_COLOURS, strict=False)),
                   caption=text.fill(PAGE, "figures.immune.caption", **numbers),
                   horizontal=True, stacked=True)  # fmt: skip
    return bars, numbers


def exons(card: Card) -> tuple[HeatmapData, dict[str, Any], annotation.Pair | None]:
    """Each exon's median share per tissue against its median over tissues, with its CV."""
    assert card.exons is not None  # noqa: S101 - the page draws this panel only with data
    shares = card.exons
    names = list(dict.fromkeys(shares["exon"].to_list()))
    order = list(dict.fromkeys(shares["tissue"].to_list()))
    middle = shares.group_by("exon").agg(mid=pl.col("share").median())
    framed = shares.join(middle, on="exon").with_columns(
        log2_ratio=pl.when((pl.col("share") > 0) & (pl.col("mid") > 0))
        .then((pl.col("share") / pl.col("mid")).log(2)).otherwise(None)
    )  # fmt: skip
    levels, swatches = diverging_levels(framed["log2_ratio"].to_list())
    frame = framed.with_columns(level=pl.Series(levels, dtype=pl.String)).select(
        "exon", "tissue", "n", "share", "log2_ratio", "level"
    )
    cv = annotation.exon_cv(shares)
    pair = annotation.qpcr_pair(cv)
    finite = {e: v for e, v in cv.items() if np.isfinite(v)}
    uneven = max(finite, key=lambda e: finite[e]) if finite else names[0]
    numbers = {"symbol": card.symbol, "release": gtex.RELEASE, "exons": len(names),
               "tissues": len(order), "dropped": f"{card.exons_dropped:,}", "uneven": uneven,
               "uneven_cv": finite.get(uneven, float("nan"))}  # fmt: skip
    if pair is not None:
        numbers |= {"first": pair.first, "second": pair.second, "pair_cv": pair.cv}
    heat = HeatmapData(frame=frame, columns=columns("exons"), row="exon", column="tissue",
                       value="log2_ratio", level="level",
                       row_label=text.word(PAGE, "figures.exons.y"),
                       column_label=text.word(PAGE, "figures.exons.x"), swatches=swatches,
                       caption=text.fill(PAGE, "figures.exons.caption", **numbers),
                       row_order=tuple(names), column_order=tuple(order),
                       row_stats=cv,
                       row_stats_label=text.word(PAGE, "figures.exons.stat"),
                       row_stats_format=".2f")  # fmt: skip
    return heat, numbers, pair
