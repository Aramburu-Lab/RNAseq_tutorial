"""The tables M9's figures are redrawn from that its page tables do not already hold (F.1).

``ctx figure`` redraws a card's figure from tables alone, never from the cache: the evidence
view needs every listed gene's M5 rows, cards left off the page included; the box its family's
tests, which span the listed genes; the exon panel how many samples had no read; the junction
panel the MANE exons its introns come from; the immune panel the HPA release it names. Each
is written before any figure, so a figure's spec finds it.
"""

from collections.abc import Sequence
from typing import Any

import polars as pl

from ctx.de.run import Writer
from ctx.genecard.cards import EVIDENCE, Card
from ctx.report import text

PAGE = "genecard"
# Each MANE model a card carries, by the source the table names it with.
SOURCES = ("gtf", "rest")


def _columns(raw: dict[str, Any]) -> dict[str, tuple[str, str]]:
    """(unit, meaning) of each column, from a content table."""
    return {c: (u, m) for c, (u, m) in raw.items()}


def _shared(*names: str) -> dict[str, tuple[str, str]]:
    """The page tables' shared column words, for these columns."""
    shared = _columns(text.raw(PAGE, "tables.columns"))
    return {n: shared[n] for n in names}


def evidence_listed(cards: Sequence[Card], evidence: pl.DataFrame) -> pl.DataFrame:
    """Every listed gene's M5 rows in every contrast, as the evidence view reads them."""
    symbols = {c.number: c.symbol for c in cards}
    ids = {c.number: c.gene_id for c in cards}
    return evidence.select(
        pl.col("number").replace_strict(symbols, return_dtype=pl.String).alias("symbol"),
        pl.col("number").replace_strict(ids, return_dtype=pl.String).alias("gene_id"),
        "contrast", *EVIDENCE,
    )  # fmt: skip


def own_family(cards: Sequence[Card], family: dict[str, Any]) -> pl.DataFrame:
    """The box's family in list order: each listed gene's test in the featured contrast."""
    symbols = {c.gene_id: c.symbol for c in cards}
    rows = [(symbols[g], g, t.n_reference, t.n_tested, t.statistic, t.p, t.r, t.p_adj)
            for g, t in family.items()]  # fmt: skip
    schema = {"symbol": pl.String, "gene_id": pl.String, "n_reference": pl.Int64,
              "n_tested": pl.Int64, "statistic": pl.Float64, "p": pl.Float64, "r": pl.Float64,
              "p_adj": pl.Float64}  # fmt: skip
    return pl.DataFrame(rows, schema=schema, orient="row")


def exons_dropped(cards: Sequence[Card]) -> pl.DataFrame:
    """Each kept card with exon shares: the GTEx samples left out for no read on the gene."""
    rows = [(c.symbol, c.gene_id, c.exons_dropped) for c in cards if c.exons is not None]
    schema = {"symbol": pl.String, "gene_id": pl.String, "exons_dropped": pl.Int64}
    return pl.DataFrame(rows, schema=schema, orient="row")


def mane(cards: Sequence[Card]) -> pl.DataFrame:
    """Each kept card's MANE Select exons from the GTF and from Ensembl REST, one row per exon."""
    rows = []
    for c in cards:
        for source, model in zip(SOURCES, (c.gtf, c.rest), strict=True):
            if model is not None:
                head = (c.symbol, c.gene_id, source, model.transcript, model.refseq, model.chrom,
                        model.strand)  # fmt: skip
                rows += [(*head, k, a, b) for k, (a, b) in enumerate(model.exons, 1)]
    schema = {"symbol": pl.String, "gene_id": pl.String, "source": pl.String,
              "transcript": pl.String, "refseq": pl.String, "chrom": pl.String,
              "strand": pl.String, "exon": pl.Int64, "start": pl.Int64,
              "end": pl.Int64}  # fmt: skip
    return pl.DataFrame(rows, schema=schema, orient="row")


def write(writer: Writer, listed: Sequence[Card], evidence: pl.DataFrame,
          kept: Sequence[Card], family: dict[str, Any], hpa_release: str) -> None:  # fmt: skip
    """Write the five tables: the list's evidence and box family, the kept cards' panels' rest.

    Args:
        writer: The run's writer.
        listed: Every card of the list, those the page could not hold included.
        evidence: ``cards.evidence`` of the whole list.
        kept: The cards the page holds, whose figures are written.
        family: The box's tests by gene ID, BH over the list (``samples.family``).
        hpa_release: The HPA release every HPA panel names.
    """
    shared = _shared("symbol", "gene_id")
    writer.table(evidence_listed(listed, evidence), "evidence_listed",
                 _shared("symbol", "gene_id", "contrast", *EVIDENCE))  # fmt: skip
    writer.table(own_family(listed, family), "own_family",
                 {**shared, **_columns(text.raw(PAGE, "figure_tables.own_family"))})  # fmt: skip
    writer.table(exons_dropped(kept), "exons_dropped",
                 {**shared, **_columns(text.raw(PAGE, "figure_tables.exons_dropped"))})  # fmt: skip
    writer.table(mane(kept), "mane", {**shared, **_columns(text.raw(PAGE, "figure_tables.mane"))})
    writer.table(pl.DataFrame({"hpa_release": [hpa_release]}), "figure_options",
                 _columns(text.raw(PAGE, "figure_tables.figure_options")))  # fmt: skip
