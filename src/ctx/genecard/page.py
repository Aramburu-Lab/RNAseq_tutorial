"""The M9 page: one group per gene card, the per-gene tables, and the sources footer.

Words come from ``content/genecard.toml``; numbers from ``run``'s cards. A source that could
not be read leaves a line naming why in place of its panels (degraded, never dropped).
"""

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import polars as pl

from ctx import runtime
from ctx.dataset import Dataset
from ctx.de.run import Writer
from ctx.genecard import annotation, figures, gtex, hpa, junctions
from ctx.genecard.cards import Card
from ctx.genecard.summary import spread, tau
from ctx.io.naming import stem
from ctx.report import text
from ctx.report.builder import (
    FigureText,
    Group,
    KeyNumber,
    ModulePage,
    ParameterCard,
    Reproduction,
    Result,
    ResultTable,
    write,
)

PAGE = "genecard"
# HPA single cell types listed per card.
TOP_CELLS = 3
PG_PER_UG = 1e6


@dataclass(frozen=True)
class Context:
    """What the page is built from: the run's options and every card's data."""

    spec: Dataset
    options: Any
    chosen: list[Card]
    evidence: pl.DataFrame
    provenance: dict[str, dict[str, str]]
    settings: runtime.Settings
    command: str
    cards: tuple[ParameterCard, ...]


def fill(key: str, **values: Any) -> str:
    """A template of this page's content file, filled."""
    return text.fill(PAGE, key, **values)


def word(key: str) -> str:
    """A fixed text of this page's content file."""
    return text.word(PAGE, key)


def _slug(card: Card) -> str:
    """A figure name for the gene: its symbol in lower case, then its Ensembl number."""
    return f"{re.sub(r'[^a-z0-9]+', '_', card.symbol.lower()).strip('_')}_{card.number}"


def _words(key: str, numbers: dict[str, Any], lead: str = "lead") -> FigureText:
    """Title, lead and reading guide of one figure, filled from its numbers."""
    return FigureText(title=fill(f"figures.{key}.title", **numbers),
                      lead=fill(f"figures.{key}.{lead}", **numbers),
                      reading=text.words(PAGE, f"figures.{key}.reading"))  # fmt: skip


def _reason(reason: str) -> str:
    """Why a gene has a card, in page words."""
    kind, _, contrast = reason.partition(" in ")
    return fill(f"reasons.{kind}", contrast=contrast)


def _unavailable(card: Card, source: str) -> str:
    """The line standing in for a source's panels."""
    return fill("gene.unavailable", source=word(f"gene.sources.{source}"), symbol=card.symbol,
                why=card.notes[source])  # fmt: skip


def _number(value: float | None, spec: str) -> str:
    """A number as a cell prints it; a missing one as the page's word for none."""
    return word("none") if value is None else format(value, spec)


def _evidence(card: Card, evidence: pl.DataFrame) -> ResultTable:
    """The gene's M5 row in every contrast."""
    rows = []
    for r in evidence.filter(pl.col("number") == card.number).iter_rows(named=True):
        tested = r["padj"] is not None
        rows.append((r["contrast"], _number(r["log2FoldChange"], ".2f"),
                     _number(r["log2FoldChange_shrunk"], ".2f"),
                     _number(r["padj"], ".2g"),
                     r["call"] if tested and r["call"] else word("evidence.untested")))  # fmt: skip
    return ResultTable(fill("evidence.title", symbol=card.symbol),
                       text.words(PAGE, "evidence.columns"), tuple(rows),
                       numeric=(1, 2, 3))  # fmt: skip


def _structure(card: Card) -> list[Result]:
    """The MANE Select transcript from the GTF and from Ensembl REST, and whether they agree."""
    rows = []
    rest = fill("structure.rest", release=annotation.ENSEMBL_RELEASE)
    for name, mane in ((word("structure.gtf"), card.gtf), (rest, card.rest)):
        if mane is not None:
            rows.append((name, mane.transcript, mane.refseq or word("none"),
                         f"{len(mane.exons)}", mane.strand))  # fmt: skip
    if not rows:
        return []
    out: list[Result] = [ResultTable(fill("structure.title", symbol=card.symbol),
                                     text.words(PAGE, "structure.columns"), tuple(rows),
                                     numeric=(3,))]  # fmt: skip
    if card.gtf is not None and card.rest is not None:
        same = (card.gtf.transcript, card.gtf.exons) == (card.rest.transcript, card.rest.exons)
        out.append(word("structure.agree" if same else "structure.differ"))
    if card.rest is not None and card.rest.refseq:
        out.append(fill("structure.primers", refseq=card.rest.refseq,
                        url=annotation.PRIMER_BLAST))  # fmt: skip
    return out


def consensus_tau(xml: hpa.Xml) -> float | None:
    """Our tau over HPA's consensus tissues."""
    return tau([v for _, v in xml.tissues])


def _amount(value: float) -> str:
    """Three significant digits, never an exponent: 2.49, 0.22, 7,900."""
    return f"{value:,.0f}" if value >= 100 else f"{value:.3g}"  # noqa: PLR2004 - 3 digits


def _lines(lines: hpa.Lines, named: hpa.CellLines) -> tuple[str, str, str]:
    """The facts row of the dataset file's cell lines: each one's nTPM, then the group's."""
    values = dict(lines.values)
    shown = [fill("facts.line", line=n, ntpm=values[n]) if n in values
             else fill("facts.line_absent", line=n) for n in named.lines]  # fmt: skip
    group = spread([v for _, v in lines.values])
    return (fill("facts.lines", group=lines.group), "; ".join(shown),
            fill("facts.lines_note", n=group.n, median=group.median))  # fmt: skip


def _facts(card: Card, release: str, named_lines: hpa.CellLines | None) -> ResultTable:
    """HPA's categories beside our tau, protein, location, blood, interactions, cell lines."""
    assert card.xml is not None and card.facts is not None  # noqa: S101 - drawn only with data
    xml, facts = card.xml, card.facts
    named = text.raw(PAGE, "facts.categories")
    rows = [(named[k], v or word("none"), "") for k, v in facts.categories.items()]
    specificity = consensus_tau(xml)
    rows.append((word("facts.tau"), _number(specificity, ".4f"),
                 fill("facts.tau_note", n=len(xml.tissues))))  # fmt: skip
    if xml.tissues:
        name, top = max(xml.tissues, key=lambda t: t[1])
        rows.append((word("facts.consensus"), fill("facts.consensus_value", tissue=name,
                                                   ntpm=top), ""))  # fmt: skip
    cells = sorted(xml.single_cell, key=lambda c: -c[1])[:TOP_CELLS]
    rows.append((word("facts.single"), "; ".join(f"{n} {v:g} nCPM" for n, v in cells)
                 or word("none"), ""))  # fmt: skip
    detected = sum(level not in {"", "not detected"} for _, level in xml.ihc)
    rows.append((word("facts.ihc"), fill("facts.ihc_value", detected=detected,
                                         tissues=len(xml.ihc)), xml.reliability))  # fmt: skip
    rows.append((word("facts.locations"), "; ".join(facts.locations[0]) or word("none"),
                 "; ".join(facts.locations[1])))  # fmt: skip
    for key, value in (("blood_im", facts.blood_im), ("blood_ms", facts.blood_ms)):
        shown = word("none") if value is None else fill("facts.concentration",
                                                         ug=_amount(value / PG_PER_UG))  # fmt: skip
        rows.append((word(f"facts.{key}"), shown, ""))
    rows.append((word("facts.interactions"), _number(facts.interactions, "d"),
                 word("facts.interactions_note")))  # fmt: skip
    if card.lines is not None and named_lines is not None:
        rows.append(_lines(card.lines, named_lines))
    rows.append((word("facts.page"), hpa.page(release, f"ENSG{card.number:011d}"), ""))
    return ResultTable(fill("facts.title", symbol=card.symbol, release=release),
                       text.words(PAGE, "facts.columns"), tuple(rows))  # fmt: skip


def _card(writer: Writer, card: Card, ctx: Context) -> Group:
    """One gene's card: own evidence, structure, then each reference panel or why it is absent."""
    release = ctx.options.hpa_release
    items: list[Any] = [_evidence(card, ctx.evidence), *_structure(card)]
    if card.renamed:
        items.insert(0, fill("gene.renamed", typed=card.renamed, release=release,
                             gene=f"ENSG{card.number:011d}", symbol=card.symbol))  # fmt: skip
    if card.tissues is not None:
        strip, numbers = figures.tissues(card)
        items.append(writer.figure(strip, f"tissues_{_slug(card)}", _words("tissues", numbers)))
    if card.exons is not None:
        heat, numbers, pair = figures.exons(card)
        lead = "lead" if pair is not None else "lead_none"
        items.append(writer.figure(heat, f"exons_{_slug(card)}",
                                   _words("exons", numbers, lead)))  # fmt: skip
    if card.junctions is not None:
        heat, numbers = figures.junctions(card)
        items.append(writer.figure(heat, f"junctions_{_slug(card)}",
                                   _words("junctions", numbers)))  # fmt: skip
    if card.xml is not None and card.facts is not None:
        bars, numbers = figures.immune(card, release)
        items.append(writer.figure(bars, f"immune_{_slug(card)}", _words("immune", numbers)))
        items.append(_facts(card, release, hpa.cell_lines(ctx.spec.info, ctx.spec.path.name)))
    items += [_unavailable(card, source) for source in card.notes]
    lead = fill("gene.lead", reason=_reason(card.reason), gene_id=card.gene_id,
                gtex_id=card.gtex_id or word("none"))  # fmt: skip
    return Group(fill("gene.title", symbol=card.symbol, gene=f"ENSG{card.number:011d}"),
                 tuple(items), lead=lead)  # fmt: skip


def _named(frame: pl.DataFrame, card: Card) -> pl.DataFrame:
    """A per-gene frame headed by the gene's symbol and our ID."""
    return frame.select(pl.lit(card.symbol).alias("symbol"), pl.lit(card.gene_id).alias("gene_id"),
                        pl.all())  # fmt: skip


def _tables(writer: Writer, ctx: Context) -> None:
    """The per-gene tables with their READMEs: every number the cards print."""
    shared = {k: tuple(v) for k, v in text.raw(PAGE, "tables.columns").items()}
    by_number = {c.number: c for c in ctx.chosen}
    evidence = ctx.evidence.with_columns(
        symbol=pl.col("number").replace_strict({n: c.symbol for n, c in by_number.items()}),
        gene_id=pl.col("number").replace_strict({n: c.gene_id for n, c in by_number.items()}),
    ).select("symbol", "gene_id", "contrast", "log2FoldChange", "log2FoldChange_shrunk",
             "padj", "call")  # fmt: skip
    writer.table(evidence, "evidence", {c: shared[c] for c in evidence.columns})
    kinds = {
        "tissues": [_named(c.tissues, c) for c in ctx.chosen if c.tissues is not None],
        "exons": [_named(c.exons, c) for c in ctx.chosen if c.exons is not None],
        "junctions": [_named(c.junctions.with_columns(gtex_id=pl.lit(c.junction_id)), c)
                      for c in ctx.chosen if c.junctions is not None],
        "immune": [_named(pl.DataFrame(c.xml.immune, schema=["cell_type", "lineage", "ntpm"],
                                       orient="row"), c)
                   for c in ctx.chosen if c.xml is not None],
        "consensus": [_named(pl.DataFrame(c.xml.tissues, schema=["tissue", "ntpm"],
                                          orient="row"), c)
                      for c in ctx.chosen if c.xml is not None],
        "lines": [_named(pl.DataFrame(c.lines.values, schema=["cell_line", "ntpm"], orient="row")
                         .with_columns(group=pl.lit(c.lines.group)), c)
                  for c in ctx.chosen if c.lines is not None],
    }  # fmt: skip
    for kind, frames in kinds.items():
        if frames:
            frame = pl.concat(frames)
            own = figures.columns(kind) if kind in {"tissues", "exons", "immune"} else {}
            meanings = {**shared, **own, **text.raw(PAGE, f"tables.{kind}_columns")}
            writer.table(frame, kind, {c: meanings[c] for c in frame.columns})


def _sources(ctx: Context) -> Group:
    """Each source with its release, retrieval date and citation; GTEx's dated acknowledgement."""
    gtex_seen = ctx.provenance["gtex"].get("retrieved", "")[:10] or word("sources.never")
    remote = ctx.provenance["remote"]
    release = ctx.options.hpa_release
    v10 = remote.get("junctions", "")[:10] or word("sources.never")
    rows = (
        (word("sources.gtex"), f"{gtex.RELEASE} ({gtex.RELEASE_DATE}), GENCODE {gtex.GENCODE}",
         gtex_seen, fill("sources.gtex_cite", date=gtex_seen)),
        (word("sources.gtex_api"), fill("sources.gtex_api_release", release=junctions.RELEASE,
                                        dataset=junctions.DATASET, gencode=junctions.GENCODE),
         v10, fill("sources.gtex_cite", date=v10)),
        (word("sources.hpa"), release, remote.get("hpa", "")[:10] or word("sources.never"),
         fill("sources.hpa_cite", versioned=hpa.page(release, "<ENSG>"))),
        (word("sources.ensembl"), annotation.ENSEMBL_RELEASE,
         remote.get("ensembl", "")[:10] or word("sources.never"), word("sources.ensembl_cite")),
        (word("sources.gencode"), ctx.spec.path.name, word("sources.never")
         if not any(c.gtf for c in ctx.chosen) else "", word("sources.gencode_cite")),
    )  # fmt: skip
    table = ResultTable(word("sources.title"), text.words(PAGE, "sources.columns"), rows)
    return Group(word("sources.title"), (table, word("sources.gtex_current"),
                                         word("sources.gtex_rerelease")))  # fmt: skip


def build(writer: Writer, ctx: Context) -> list[Path]:
    """Write every table and figure, then the page; the page first in the returned list."""
    _tables(writer, ctx)
    results: list[Result] = [_card(writer, c, ctx) for c in ctx.chosen]
    results.append(_sources(ctx))
    gtex_numbers = ctx.provenance["gtex"]
    tissues = next((c.tissues for c in ctx.chosen if c.tissues is not None), None)
    keys = (KeyNumber(word("keys.genes"), f"{len(ctx.chosen):,}"),
            KeyNumber(word("keys.samples"), gtex_numbers.get("samples", word("none"))),
            KeyNumber(word("keys.tissues"), word("none") if tissues is None else
                      f"{tissues['tissue'].n_unique():,}"),
            KeyNumber(word("keys.hpa"), ctx.options.hpa_release))  # fmt: skip
    usage, budget = runtime.measure(ctx.settings.started), runtime.budget(PAGE)
    built = ModulePage(
        title=fill("title", dataset=ctx.spec.name), **text.skeleton(PAGE),
        parameters=ctx.cards, results=tuple(results), key_numbers=keys,
        reproduce=Reproduction(
            command=ctx.command, wall_seconds=usage.wall_seconds,
            peak_rss_self=usage.peak_rss_self, peak_rss_children=usage.peak_rss_children,
            profile=f"{ctx.settings.profile}, {ctx.settings.threads} threads",
            budget=fill("budget", minutes=budget.minutes, ram_gb=budget.ram_gb),
            versions=runtime.versions(),
        ),
    )  # fmt: skip
    page = writer.out / f"{stem(PAGE, f'{writer.dataset}_overview', 1, day=writer.day)}.html"
    return [write(built, page), *writer.files]
