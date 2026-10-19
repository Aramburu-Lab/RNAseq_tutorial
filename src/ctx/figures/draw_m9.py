"""M9's figures redrawn from the tables ``ctx genecard`` wrote, by M9's own functions.

A card is rebuilt from its rows of M9's tables, never from the response cache or GTEx's
files: the tissue spread, exon shares, junction medians and immune nTPM as the page tables
hold them, the MANE exons, dropped samples and HPA release from the tables written for the
figures (``ctx.genecard.figure_tables``). The own-sample panels read M5's normalised values
as M9 wrote them, its contrasts named by their words, and the box its family's tests, which
span the whole card list. Each drawer calls the page's own function for its panel.
"""

from collections.abc import Mapping
from dataclasses import replace
from pathlib import Path
from typing import Any

import polars as pl

from ctx import runtime
from ctx.dataset import Dataset
from ctx.figures.kind import NoFigureError, Plain
from ctx.figures.kinds_m9 import OneCard
from ctx.figures.recorder import Recorder
from ctx.genecard import annotation, boxes, figures, hpa, page, samples
from ctx.genecard.cards import EVIDENCE, Card
from ctx.report import genes
from ctx.report.builder import FigureText

FOLDER = "11_genecard"
# A run's settings, which the page's functions for one panel never read.
SETTINGS = runtime.Settings("laptop", Path(), 1, None, 0.0)


def _read(inputs: Mapping[str, Path], name: str) -> pl.DataFrame:
    """One of M9's tables, read back."""
    return pl.read_parquet(inputs[name])


def _numbered(frame: pl.DataFrame) -> pl.DataFrame:
    """A per-gene table with each gene's unversioned Ensembl number, as M9's cards key them."""
    ids = frame["gene_id"].unique(maintain_order=True).to_list()
    numbers = {g: genes.number(g) for g in ids}
    return frame.with_columns(
        number=pl.col("gene_id").replace_strict(numbers, return_dtype=pl.Int64)
    )


def _card(frame: pl.DataFrame, p: OneCard) -> tuple[Card, pl.DataFrame]:
    """The card ``p`` names and its rows of a per-gene table, ``symbol`` and ``gene_id`` dropped.

    Raises:
        NoFigureError: When the table holds no card of that name.
    """
    number = int(p.gene.rsplit("_", 1)[1])
    rows = _numbered(frame).filter(pl.col("number") == number)
    if rows.is_empty():
        raise NoFigureError(f"{p.gene}: no card of this gene in M9's tables")
    card = Card(number, rows["gene_id"][0], rows["symbol"][0], "")
    if page._slug(card) != p.gene:
        raise NoFigureError(f"{p.gene}: M9's tables name this gene {page._slug(card)}")
    return card, rows.drop("symbol", "gene_id", "number")


def _plain(words: FigureText) -> FigureText:
    """Words as a redraw prints them: no link back to a module page, which is not beside it."""
    return replace(words, back="")


def evidence_view(inputs: Mapping[str, Path], _: Plain, spec: Dataset) -> tuple[Any, FigureText]:
    """Every listed gene in every contrast, cards left off the page included.

    Raises:
        NoFigureError: For a dataset of one contrast, where M9 draws no evidence view.
    """
    listed = _numbered(_read(inputs, "evidence_listed"))
    who = listed.select("number", "gene_id", "symbol").unique(maintain_order=True)
    cards = [Card(n, g, s, "") for n, g, s in who.rows()]
    evidence = listed.select("number", "contrast", *EVIDENCE)
    ctx = page.Context(spec=spec, options=None, chosen=cards, evidence=evidence, provenance={},
                       settings=SETTINGS, command="", cards=())  # fmt: skip
    if (found := page.evidence_view(Recorder(FOLDER, spec.slug), ctx)) is None:
        raise NoFigureError(f"{spec.slug}: one contrast, so M9 draws no evidence view")
    heat, words = found.figure[0], found.figure[2]
    return heat, _plain(words)


def _own(inputs: Mapping[str, Path], p: OneCard, spec: Dataset, kind: str,
         evidence: pl.DataFrame, family: dict[str, samples.Test]
         ) -> tuple[Any, FigureText]:  # fmt: skip
    """One own-sample panel, from M5's values as M9 wrote them, by the page's ``_own``.

    Raises:
        NoFigureError: When the page draws no such panel of the gene.
    """
    names = {label: name for name, label in boxes.names(spec).items()}
    table = _read(inputs, "own_values")
    card, _ = _card(table, p)
    frame = table.with_columns(contrast=pl.col("contrast").replace_strict(names)).select(
        "gene_id", "contrast", "sample", "donor", "group", "batch", "value"
    )
    paired = {c.name: c.paired for c in spec.contrasts}
    own = samples.Values(frame, {c: paired.get(c, False) for c in frame["contrast"].unique()})
    ctx = page.Context(spec=spec, options=None, chosen=[card], evidence=evidence, provenance={},
                       settings=SETTINGS, command="", cards=(), own=own,
                       family=family)  # fmt: skip
    drawn = {k: (data, words) for k, data, words, _ in page._own(card, ctx)[1]}
    if kind not in drawn:
        raise NoFigureError(f"{p.gene}: M9 draws no {kind} panel of this gene from these tables")
    data, words = drawn[kind]
    return data, _plain(words)


def box(inputs: Mapping[str, Path], p: OneCard, spec: Dataset) -> tuple[Any, FigureText]:
    """The gene in the featured contrast's samples, its p from the card list's family."""
    evidence = _numbered(_read(inputs, "evidence")).select("number", "contrast", *EVIDENCE)
    family = {row["gene_id"]: samples.Test(row["n_reference"], row["n_tested"], row["statistic"],
                                           row["p"], row["r"], row["p_adj"])
              for row in _read(inputs, "own_family").iter_rows(named=True)}  # fmt: skip
    return _own(inputs, p, spec, "box", evidence, family)


def contrasts(inputs: Mapping[str, Path], p: OneCard, spec: Dataset) -> tuple[Any, FigureText]:
    """F1: the gene in every contrast, every sample's first page of contrasts."""
    # The box beside it needs the family and M5's rows; with neither, the page's _own skips it.
    none = pl.DataFrame(schema={"number": pl.Int64, "contrast": pl.String, "padj": pl.Float64})
    return _own(inputs, p, spec, "contrasts", none, {})


def tissues(inputs: Mapping[str, Path], p: OneCard, _: Dataset) -> tuple[Any, FigureText]:
    """F2: GTEx median per tissue, both sexes side by side, as built."""
    card, spread = _card(_read(inputs, "tissues"), p)
    card.tissues = spread
    data, numbers = figures.tissue_views(card)[0][1:]
    return data, _plain(page._words("tissues", numbers))


def exons(inputs: Mapping[str, Path], p: OneCard, _: Dataset) -> tuple[Any, FigureText]:
    """F4: each exon's share per tissue that expresses the gene.

    Raises:
        NoFigureError: Where too few tissues express the gene, as the page then prints a line.
    """
    card, shares = _card(_read(inputs, "exons"), p)
    card.exons = shares
    card.tissues = _card(_read(inputs, "tissues"), p)[1]
    dropped = _card(_read(inputs, "exons_dropped"), p)[1]
    card.exons_dropped = int(dropped["exons_dropped"][0])
    if isinstance(drawn := page._exon_panel(card), str):
        raise NoFigureError(drawn)
    data, words = drawn
    return data, _plain(words)


def _mane(frame: pl.DataFrame, source: str) -> annotation.Mane | None:
    """A card's MANE model from one source, exons in transcription order; None without one."""
    rows = frame.filter(pl.col("source") == source).sort("exon")
    if rows.is_empty():
        return None
    first = rows.row(0, named=True)
    return annotation.Mane(first["transcript"], first["refseq"], first["chrom"], first["strand"],
                           tuple(zip(rows["start"], rows["end"], strict=True)))  # fmt: skip


def junctions(inputs: Mapping[str, Path], p: OneCard, _: Dataset) -> tuple[Any, FigureText]:
    """GTEx V10's junctions of the gene, MANE introns marked from the model the page read."""
    card, rows = _card(_read(inputs, "junctions"), p)
    card.junction_id, card.junctions = rows["gtex_id"][0], rows.drop("gtex_id")
    mane = _numbered(_read(inputs, "mane")).filter(pl.col("number") == card.number)
    card.gtf, card.rest = _mane(mane, "gtf"), _mane(mane, "rest")
    heat, numbers = figures.junctions(card)
    return heat, _plain(page._words("junctions", numbers))


def immune(inputs: Mapping[str, Path], p: OneCard, _: Dataset) -> tuple[Any, FigureText]:
    """F3: HPA immune cell nTPM of the gene, at the HPA release the run named."""
    card, rows = _card(_read(inputs, "immune"), p)
    card.xml = hpa.Xml("", (), tuple(rows.select("cell_type", "lineage", "ntpm").rows()), (), (),
                       "")  # fmt: skip
    release = _read(inputs, "genecard_options")["hpa_release"][0]
    bars, numbers = figures.immune(card, release)
    return bars, _plain(page._words("immune", numbers))
