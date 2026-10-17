"""The drawn panels of a card: GTEx tissues, HPA immune cells, GTEx exon use, V10 junctions.

Each returns a figure contract and the numbers its words are filled with; the page writes
them. Every panel is descriptive: n, median, IQR, tau, CV; external data get no test.
"""

import math
from collections.abc import Mapping
from typing import Any

import polars as pl

from ctx.de.stats import ALPHA, LFC_MIN
from ctx.genecard import annotation, gtex, hpa
from ctx.genecard.cards import Card
from ctx.genecard.junctions import RELEASE as RELEASE_JUNCTIONS
from ctx.genecard.summary import DETECTED, tau
from ctx.report import text
from ctx.viz.axes import axis_title
from ctx.viz.bar import BarData
from ctx.viz.heatmap import HeatmapData, Swatch, diverging_levels
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


def tissues(card: Card, donors: tuple[str, ...] = SEXES) -> tuple[BarData, dict[str, Any]]:
    """GTEx median per tissue of ``donors``, as bars, tissues by their all-donor median.

    ``donors`` is both sexes side by side (the figure as built), ``(gtex.ALL,)``, or one sex:
    F2's views (D-056: summaries, so the donor control is three precomputed options). Every
    view keeps the all-donor order, so a bar that moves between views is the sex's own; its
    lead and caption are the view's own numbers. A sex with no sample in a tissue has no bar.

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
        spread.filter(pl.col("sex").is_in(donors))
        .with_columns(value=pl.col("median").map_elements(_log, return_dtype=pl.Float64),
                      slot=pl.col("tissue").replace_strict(rank, return_dtype=pl.Int64))
        .sort("slot", "sex").select("tissue", "sex", "value", "n", "median", "q1", "q3")
    )  # fmt: skip
    # Both sexes side by side are read against the all-donor medians that order them.
    own = every if donors == SEXES else bars.sort("median", descending=True)
    numbers = {"symbol": card.symbol, "release": gtex.RELEASE, "top": own["tissue"][0],
               "top_median": own["median"][0], "bottom": own["tissue"][-1],
               "bottom_median": own["median"][-1], "tau": shown(tau(own["median"].to_numpy())),
               "tissues": own.height, "samples": f"{int(own['n'].sum()):,}",
               "least": own["n"].min(), "most": own["n"].max(),
               "donors": text.word(PAGE, f"figures.tissues.donors.{'_'.join(donors)}")}  # fmt: skip
    palette = category_colours(list(SEXES))
    data = BarData(frame=bars, columns=columns("tissues"), category="tissue", value="value",
                   category_label=text.word(PAGE, "figures.tissues.x"),
                   value_label=text.word(PAGE, "figures.tissues.y"), series="sex",
                   colours={d: palette.get(d) for d in donors},
                   caption=text.fill(PAGE, "figures.tissues.caption", **numbers),
                   horizontal=True, hover=("n", "median", "q1", "q3"))  # fmt: skip
    return data, numbers


# F2's views, the figure as built first (D-056, figure controls plan).
DONOR_VIEWS = (SEXES, (gtex.ALL,), ("female",), ("male",))


def tissue_views(card: Card) -> list[tuple[str, BarData, dict[str, Any]]]:
    """(label, bars, numbers) of each donor view a reader may choose; one view per sex present."""
    assert card.tissues is not None  # noqa: S101 - the page draws this panel only with data
    present = set(card.tissues["sex"].to_list())
    return [(text.word(PAGE, f"figures.tissues.views.{'_'.join(d)}"), *tissues(card, d))
            for d in DONOR_VIEWS if set(d) & present]  # fmt: skip


def immune_tau(xml: hpa.Xml) -> float | None:
    """Our tau over HPA's immune cell types, the pooled PBMC total left out."""
    return tau([v for cell, _, v in xml.immune if cell != hpa.POOLED])


def immune_subset(rows: list[tuple[str, str, float]], lineages: set[str],
                  log: bool = True) -> dict[str, Any]:  # fmt: skip
    """F3's numbers over the cell types of ``lineages``: tau on the scale, the top cell type.

    The page's twin is ``ctxM9.immune``; ``rows`` are (cell type, lineage, nTPM), the pooled
    total already left out. With no cell type left, tau and the top are None.
    """
    kept = [(c, v) for c, lineage, v in rows if lineage in lineages]
    top = max(kept, key=lambda r: r[1]) if kept else None
    return {"tau": tau([v for _, v in kept], log), "cells": len(kept),
            "top": top[0] if top else None, "top_ntpm": top[1] if top else None}  # fmt: skip


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
                   horizontal=True, stacked=True,
                   live=immune_controls(lineages, frame.height))  # fmt: skip
    return bars, numbers


def immune_controls(lineages: list[str], cells: int) -> dict[str, Any]:
    """F3's checkbox set: the lineages kept and the scale tau is taken on (``ctxM9.immune``).

    One bar trace per lineage, in ``lineages``' order, as ``BarData`` draws a series; the page
    hides a lineage's trace and recomputes tau and the top cell type over what is left. The
    lead is the build's own template, filled in the page.
    """
    words = text.raw(PAGE, "figures.immune.live")
    start = {"cells": cells, "total": cells, "kept": len(lineages), "lineages": len(lineages),
             "scale": words["log"]}  # fmt: skip
    return {"kind": "subset", "twin": "immune", "options": list(lineages),
            "traces": list(lineages),
            "scales": [["log", words["log"]], ["linear", words["linear"]]],
            "template": text.raw(PAGE, "figures.immune.lead"), "none": text.word(PAGE, "none"),
            "words": words,
            "summary": text.fill(PAGE, "figures.immune.live.summary", **start)}  # fmt: skip


def least_even(cv: Mapping[str, float | None]) -> tuple[str, float] | None:
    """The exon with the largest defined CV and that CV; None where no CV is defined.

    The first of tied exons in ``cv``'s order; the page's twin is ``ctxM9.uneven``.
    """
    finite = {e: v for e, v in cv.items() if v is not None}
    if not finite:
        return None
    exon = max(finite, key=lambda e: finite[e])
    return exon, finite[exon]


def exons(card: Card) -> tuple[HeatmapData, dict[str, Any], annotation.Pair | None] | None:
    """Each exon's median share per tissue against its median over the tissues that express it.

    Only tissues at ``DETECTED`` TPM or more are drawn (``gtex.detected``, H.3); fewer than two
    leave no panel, and None. The CV is printed only where it is defined.
    """
    assert card.exons is not None and card.tissues is not None  # noqa: S101 - drawn with data
    use = annotation.exon_use(card.exons, gtex.detected(card.tissues))
    if len(use.tissues) < 2:  # noqa: PLR2004 - a reference across tissues needs two
        return None
    names = list(dict.fromkeys(use.frame["exon"].to_list()))
    levels, swatches = diverging_levels(list(use.steps))
    frame = use.frame.with_columns(level=pl.Series(levels, dtype=pl.String)).select(
        "exon", "tissue", "n", "share", "reference", "log2_ratio", "level"
    )
    pair = annotation.qpcr_pair(use.cv)
    finite = {e: v for e, v in use.cv.items() if v is not None}
    uneven = (least_even(use.cv) or (names[0], math.nan))[0]
    unreferenced = frame.filter(pl.col("reference") == 0)["exon"].n_unique()
    numbers = {"symbol": card.symbol, "release": gtex.RELEASE, "exons": len(names),
               "tissues": len(use.tissues), "left": len(use.left_out), "floor": DETECTED,
               "unreferenced": unreferenced, "uneven": uneven,
               "dropped": text.fill(PAGE, "figures.exons.dropped_one" if card.exons_dropped == 1
                                    else "figures.exons.dropped_many", n=card.exons_dropped),
               "uneven_cv": finite.get(uneven, float("nan"))}  # fmt: skip
    if pair is not None:
        numbers |= {"first": pair.first, "second": pair.second, "pair_cv": pair.cv}
    heat = HeatmapData(frame=frame, columns=columns("exons"), row="exon", column="tissue",
                       value="log2_ratio", level="level",
                       row_label=text.word(PAGE, "figures.exons.y"),
                       column_label=text.word(PAGE, "figures.exons.x"), swatches=swatches,
                       caption=text.fill(PAGE, "figures.exons.caption", **numbers),
                       row_order=tuple(names), column_order=use.tissues,
                       row_stats=finite,
                       row_stats_label=text.word(PAGE, "figures.exons.stat"),
                       row_stats_format=".2f")  # fmt: skip
    return heat, numbers, pair


def exon_controls(heat: HeatmapData) -> dict[str, Any]:
    """F4's checkbox set: the tissues kept, over which ``ctxM9.exonUse`` retakes the panel.

    Every drawn tissue starts kept. A tissue left out keeps its column, drawn empty, so the
    plot keeps the size it was laid out at. The page's heatmap has no PDF and SVG of a changed
    view, so the live line says the files draw every tissue (D-051).
    """
    unit = heat.columns[heat.value][0]
    tissues = len(heat.columns_in_order())
    return {"kind": "subset", "twin": "exons", "options": heat.columns_in_order(),
            "exons": heat.rows_in_order(), "template": text.raw(PAGE, "figures.exons.lead"),
            "template_none": text.raw(PAGE, "figures.exons.lead_none"),
            "words": text.raw(PAGE, "figures.exons.live"),
            "hover": [heat.row_label, heat.column_label, axis_title(heat.value, unit)],
            "levels": [s.level for s in heat.swatches],
            "summary": text.fill(PAGE, "figures.exons.live.summary", kept=tissues,
                                 tissues=tissues)}  # fmt: skip


# The evidence view draws at most this many genes; the rest are listed in a table (GS.9).
EVIDENCE_ROWS = 100
CALLS = ("up", "down")


def evidence_view(
    evidence: pl.DataFrame, symbols: Mapping[int, str], labels: Mapping[str, str]
) -> tuple[HeatmapData, dict[str, Any], list[int]]:
    """The listed genes by contrasts: a called cell in its log2 FC's step, an uncalled one grey.

    ``evidence`` is ``cards.evidence``'s long table, genes in list order; ``labels`` each
    contrast's words. A gene untested in a contrast draws nothing there. The first
    ``EVIDENCE_ROWS`` genes are drawn; the rest are returned for a table.

    Returns:
        The figure, its numbers, and the numbers of the genes left out of it.
    """
    order = list(dict.fromkeys(evidence["number"].to_list()))
    drawn, rest = order[:EVIDENCE_ROWS], order[EVIDENCE_ROWS:]
    frame = evidence.filter(pl.col("number").is_in(drawn)).with_columns(
        gene=pl.col("number").replace_strict(dict(symbols), return_dtype=pl.String),
        contrast=pl.col("contrast").replace_strict(dict(labels), return_dtype=pl.String),
    )
    called = frame["call"].is_in(CALLS).fill_null(value=False).to_list()
    steps, swatches = diverging_levels(frame["log2FoldChange"].to_list())
    none = text.word(PAGE, "figures.evidence.not_called")
    levels = [s if c else (none if p is not None else None)
              for s, c, p in zip(steps, called, frame["padj"].to_list(), strict=True)]  # fmt: skip
    used = set(levels)
    keys = (*(w for w in swatches if w.level in used), Swatch(none, None))
    frame = frame.with_columns(level=pl.Series(levels, dtype=pl.String)).select(
        "gene", "contrast", "log2FoldChange", "log2FoldChange_shrunk", "padj", "call", "level"
    )
    numbers = {"genes": len(drawn), "listed": len(order), "contrasts": len(labels),
               "called": sum(called), "cells": frame["padj"].is_not_null().sum(),
               "alpha": ALPHA, "lfc_min": LFC_MIN}  # fmt: skip
    heat = HeatmapData(frame=frame, columns=columns("evidence"), row="gene", column="contrast",
                       value="log2FoldChange", level="level",
                       row_label=text.word(PAGE, "figures.evidence.y"),
                       column_label=text.word(PAGE, "figures.evidence.x"), swatches=keys,
                       caption=text.fill(PAGE, "figures.evidence.caption", **numbers),
                       row_order=tuple(symbols[n] for n in drawn),
                       column_order=tuple(labels.values()))  # fmt: skip
    return heat, numbers, rest


def junctions(card: Card) -> tuple[HeatmapData, dict[str, Any]]:
    """GTEx V10 junction medians per tissue against each junction's median over tissues.

    A junction is labelled by its coordinates, with a mark where it is an intron of the GTF's
    MANE Select transcript base for base: the GTF's model, or Ensembl REST's where the GTF tags
    none (GDC's GENCODE v36 has no MANE Select for NDUFA4L2; release 116 has).
    """
    assert card.junctions is not None  # noqa: S101 - the page draws this panel only with data
    found = card.junctions
    model = card.gtf or card.rest
    introns = set(model.introns()) if model is not None else set()
    order = list(dict.fromkeys(found.sort("start", "end")["junction"].to_list()))
    tissues = list(dict.fromkeys(found["tissue"].to_list()))
    middle = found.group_by("junction").agg(mid=pl.col("median").median())
    framed = found.join(middle, on="junction").with_columns(
        mane=pl.struct("start", "end").map_elements(lambda r: (r["start"], r["end"]) in introns,
                                                    return_dtype=pl.Boolean),
        log2_ratio=pl.when((pl.col("median") > 0) & (pl.col("mid") > 0))
        .then((pl.col("median") / pl.col("mid")).log(2)).otherwise(None),
    )  # fmt: skip
    levels, swatches = diverging_levels(framed["log2_ratio"].to_list())
    frame = framed.with_columns(level=pl.Series(levels, dtype=pl.String)).select(
        "junction", "tissue", "median", "mane", "log2_ratio", "level"
    )
    mane = set(framed.filter(pl.col("mane"))["junction"].to_list())
    marker = text.word(PAGE, "figures.junctions.mane")
    # Short labels in genomic order: full coordinates left 54 tissue columns 80 px at 1280 px.
    display = {j: f"J{k} {marker}" if j in mane else f"J{k}" for k, j in enumerate(order, 1)}
    mids = dict(middle.iter_rows())
    numbers = {"symbol": card.symbol, "release": RELEASE_JUNCTIONS, "junctions": len(order),
               "tissues": len(tissues), "matched": len(mane), "introns": len(introns),
               "gencode_id": card.junction_id}  # fmt: skip
    heat = HeatmapData(frame=frame, columns=columns("junctions"), row="junction", column="tissue",
                       value="log2_ratio", level="level",
                       row_label=text.word(PAGE, "figures.junctions.y"),
                       column_label=text.word(PAGE, "figures.junctions.x"), swatches=swatches,
                       caption=text.fill(PAGE, "figures.junctions.caption", **numbers),
                       row_order=tuple(order), column_order=tuple(tissues),
                       row_stats={j: float(mids[j]) for j in order},
                       row_stats_label=text.word(PAGE, "figures.junctions.stat"),
                       row_stats_format=".3g", row_display=display)  # fmt: skip
    return heat, numbers
