"""The M9 page: one group per gene card, the per-gene tables, and the sources footer.

Words come from ``content/genecard.toml``; numbers from ``run``'s cards. A source that could
not be read leaves a line naming why in place of its panels (degraded, never dropped).
"""

import re
from collections import Counter
from collections.abc import Callable
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

import polars as pl

from ctx import runtime
from ctx.dataset import Dataset
from ctx.de.page import every
from ctx.de.run import Writer
from ctx.de.stats import ALPHA, LFC_MIN
from ctx.genecard import (
    annotation,
    boxes,
    figure_tables,
    figures,
    gtex,
    hpa,
    junctions,
    listing,
    samples,
)
from ctx.genecard.cards import Card
from ctx.genecard.summary import DETECTED, spread, tau
from ctx.io.naming import stem
from ctx.report import text
from ctx.report.builder import (
    MAX_BYTES,
    Code,
    FigurePage,
    FigureText,
    Group,
    KeyNumber,
    ModulePage,
    Panel,
    ParameterCard,
    Reproduction,
    Result,
    ResultTable,
    packed_bytes,
    rendered,
    write,
)
from ctx.report.entry import write_entry
from ctx.report.lookups import lookup
from ctx.viz.save import figure_block, write_figure

PAGE = "genecard"
# HPA single cell types listed per card.
TOP_CELLS = 3
PG_PER_UG = 1e6


@dataclass(frozen=True)
class Context:
    """What the page is built from: the run's options and every card's data.

    ``tested`` is every tested gene's M5 evidence for the list from a link
    (``listing.evidence``); None leaves the page without that panel. ``own`` is the listed
    genes' M5 normalised values (``samples.read``), None for none, and ``family`` the box's
    tests in the featured contrast, BH over the card list (``samples.family``).
    """

    spec: Dataset
    options: Any
    chosen: list[Card]
    evidence: pl.DataFrame
    provenance: dict[str, dict[str, str]]
    settings: runtime.Settings
    command: str
    cards: tuple[ParameterCard, ...]
    tested: dict[str, Any] | None = None
    own: samples.Values | None = None
    family: dict[str, samples.Test] = field(default_factory=dict)


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


@dataclass(frozen=True)
class Drawn:
    """One gene's card as the page prints it, and the figures whose files it links.

    The files are written only for a card the page keeps (``build``): 50 TCGA cards drawn and
    written whole took ten minutes for a page the cap then refused (C3.1).
    """

    card: Card
    items: tuple[Any, ...]
    figures: tuple[tuple[Any, str, FigureText, dict[str, Any] | None], ...]


def _card(writer: Writer, card: Card, ctx: Context) -> Drawn:
    """One gene's card: own evidence, structure, then each reference panel or why it is absent."""
    release = ctx.options.hpa_release
    items: list[Any] = [*_renamed(card, release), _evidence(card, ctx.evidence),
                        *_structure(card)]  # fmt: skip
    panels: list[tuple[Any, str, FigureText, dict[str, Any] | None]] = []

    def panel(data: Any, kind: str, words: FigureText,
              views: list[tuple[str, Any, FigureText]] | None = None,
              controls: dict[str, Any] | None = None) -> None:  # fmt: skip
        """A figure of the card: its page block now, its files once the card is kept.

        With ``views``, the first being the figure as built, the block carries a select
        between them (the D-051 view kind, each option with its own PDF and SVG through
        ``lookups.lookup``: ``views.viewed`` read a box's caption object as text and stopped the
        run); the files draw the first. ``controls`` are a
        figure's page controls where its contract carries none, as F4's heatmap.
        """
        what, words = (
            f"{kind}_{_slug(card)}",
            replace(words, back=f"{writer.stem('overview')}.html"),
        )
        panels.append((data, what, words, controls))
        block = figure_block(data, writer.stem(what), words)
        if views is not None and len(views) > 1:
            block = lookup(block, views, word(f"figures.{kind}.view"))
        if controls is not None:
            block = replace(block, controls=controls)
        items.append(block)

    lines, own = _own(card, ctx)
    items += lines
    for kind, data, words, choices in own:
        panel(data, kind, words, choices)
    if card.tissues is not None:
        views = [(label, data, _words("tissues", numbers))
                 for label, data, numbers in figures.tissue_views(card)]  # fmt: skip
        panel(views[0][1], "tissues", views[0][2], views)
    if card.exons is not None:
        if isinstance(exon := _exon_panel(card), str):
            items.append(exon)
        else:
            panel(exon[0], "exons", exon[1], controls=figures.exon_controls(exon[0]))
    if card.junctions is not None:
        heat, numbers = figures.junctions(card)
        panel(heat, "junctions", _words("junctions", numbers))
    if card.xml is not None and card.facts is not None:
        bars, numbers = figures.immune(card, release)
        panel(bars, "immune", _words("immune", numbers))
        items.append(_facts(card, release, hpa.cell_lines(ctx.spec.info, ctx.spec.path.name)))
    items += [_unavailable(card, source) for source in card.notes]
    return Drawn(card, tuple(items), tuple(panels))


def _renamed(card: Card, release: str) -> list[str]:
    """The line naming the symbol HPA gave a gene the table lacks, where it was renamed."""
    if not card.renamed:
        return []
    return [fill("gene.renamed", typed=card.renamed, release=release,
                 gene=f"ENSG{card.number:011d}", symbol=card.symbol)]  # fmt: skip


def _exon_panel(card: Card) -> tuple[Any, FigureText] | str:
    """F4's figure and words, or the line standing in where too few tissues express the gene."""
    if (drawn := figures.exons(card)) is None:
        return fill("figures.exons.too_few", symbol=card.symbol, floor=DETECTED)
    heat, numbers, pair = drawn
    return heat, _words("exons", numbers, "lead" if pair is not None else "lead_none")


def _own(card: Card, ctx: Context) -> tuple[list[str], list[tuple[str, Any, FigureText, Any]]]:
    """The card's own-sample panels: the featured contrast's box, then F1 with its batch views.

    Each panel is (kind, figure, words, views); a box the family holds no test for is a line
    printed in its place.
    """
    if ctx.own is None or ctx.spec.de is None:
        return [], []
    labels = every(ctx.spec.contrasts, ctx.spec.columns)
    featured = ctx.spec.de.featured
    row = ctx.evidence.filter((pl.col("number") == card.number)
                              & (pl.col("contrast") == featured))  # fmt: skip
    padj = row["padj"][0] if row.height else None
    out: list[tuple[str, Any, FigureText, Any]] = []
    lines = []
    drawn = boxes.box(card, ctx.own, featured, labels[featured], ctx.family, padj)
    if drawn is not None:
        out.append(("box", drawn[0], _words("box", drawn[1]), None))
    else:
        lines.append(fill("figures.box.none", symbol=card.symbol, contrast=labels[featured].label))
    if len(ctx.spec.contrasts) > 1:
        names = boxes.names(ctx.spec)
        found = boxes.views(card, ctx.own, names, samples.batches(ctx.own, boxes.MAX_BATCHES))
        choices = [(name, data, _words("contrasts", numbers)) for name, data, numbers in found]
        if choices:
            out.append(("contrasts", choices[0][1], choices[0][2], choices))
    return lines, out


def _group(drawn: Drawn, packed: bool) -> Group:
    """A card's group; packed, its figures travel gzipped and are drawn when reached."""
    card = drawn.card
    lead = fill("gene.lead", reason=_reason(card.reason), gene_id=card.gene_id,
                gtex_id=card.gtex_id or word("none"))  # fmt: skip
    return Group(fill("gene.title", symbol=card.symbol, gene=f"ENSG{card.number:011d}"),
                 drawn.items, lead=lead, packed=packed)  # fmt: skip


def _argv(ctx: Context) -> list[str]:
    """This run's ``ctx genecard`` arguments, its gene list apart."""
    argv = ["--root", str(ctx.settings.root), "genecard", "--dataset", ctx.spec.slug,
            "--dataset-file", str(ctx.options.dataset_file)]  # fmt: skip
    if ctx.options.hpa_release != hpa.RELEASE:
        argv += ["--hpa-release", ctx.options.hpa_release]
    if ctx.options.offline:
        argv.append("--offline")
    return argv


def left_out_command(ctx: Context, cards: list[Card]) -> str:
    """The command that builds the cards a page could not hold, on a page of their own.

    Each gene by its Ensembl ID: a symbol can name two genes, and an ID cannot.
    """
    return runtime.retypable([*_argv(ctx), "--genes",
                              ",".join(f"ENSG{c.number:011d}" for c in cards)])  # fmt: skip


def _list_panel(ctx: Context, drawn: list["Drawn"], first: int) -> Group | None:
    """The list from a link: every tested gene's evidence, and where each kept card sits.

    ``first`` is the group number of the first card, as the page numbers its groups.
    """
    if ctx.tested is None:
        return None
    labels = {n: w.label for n, w in every(ctx.spec.contrasts, ctx.spec.columns).items()}
    contrasts = [(c.name, labels[c.name]) for c in ctx.spec.contrasts]
    featured = ctx.spec.de.featured if ctx.spec.de is not None else contrasts[0][0]
    cards = [(d.card.number, f"group-{first + i}", d.card.symbol) for i, d in enumerate(drawn)]
    calls = {code: name for name, code in listing.CALL_CODES.items()}
    calls[listing.UNTESTED] = word("evidence.untested")
    words = {"list": text.raw(PAGE, "list"), "picker": text.raw("ui", "picker"),
             "errors": text.raw("ui", "select"), "calls": calls}  # fmt: skip
    words["list"] = {**words["list"], "heat_caption": fill(
        "list.heat_caption", alpha=ALPHA, lfc_min=LFC_MIN, contrasts=len(contrasts),
        genes="{genes}", listed="{listed}")}  # fmt: skip
    words["heat"] = figures.live_scale([label for _, label in contrasts])
    settings = listing.controls(ctx.tested, contrasts, featured, cards,
                                runtime.retypable(_argv(ctx)), words)  # fmt: skip
    return Group(word("list.group"), (Panel("genelist", settings),), lead=word("list.lead"))


def _left_out(ctx: Context, cards: list[Card]) -> Group:
    """The genes of the card list the page could not hold, and the command for their page."""
    names = ", ".join(c.symbol for c in cards)
    lead = fill("left_out.lead", n=len(cards), listed=len(ctx.chosen), cap=MAX_BYTES / 1e6)
    return Group(word("left_out.title"), (fill("left_out.genes", genes=names),
                                          Code(word("left_out.command"),
                                               left_out_command(ctx, cards))),
                 lead=lead)  # fmt: skip


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


@dataclass(frozen=True)
class Evidence:
    """The genes-by-contrasts view of the whole card list, and its figure for the files."""

    group: Group
    figure: tuple[Any, str, FigureText]


def _shown_names(cards: list[Card]) -> dict[int, str]:
    """Each card's row name in the evidence view: its symbol, with its ID where two share one."""
    seen = Counter(c.symbol for c in cards)
    return {c.number: c.symbol if seen[c.symbol] == 1 else f"{c.symbol} {c.gene_id}" for c in cards}


def evidence_view(writer: Writer, ctx: Context) -> Evidence | None:
    """Every listed gene in every contrast (GS.9's tier 1), before the cards; None for one contrast.

    It covers the whole list, cards left out included, since its numbers are M5's alone. A
    dataset with one contrast has its one column in each card's table already.
    """
    if len(ctx.spec.contrasts) < 2:  # noqa: PLR2004 - a view across contrasts needs two
        return None
    labels = {n: w.label for n, w in every(ctx.spec.contrasts, ctx.spec.columns).items()}
    names = _shown_names(ctx.chosen)
    heat, numbers, rest = figures.evidence_view(ctx.evidence, names, labels)
    words = replace(_words("evidence", numbers), back=f"{writer.stem('overview')}.html")
    items: list[Any] = [figure_block(heat, writer.stem("evidence_view"), words)]
    if rest:
        rows = []
        for n in rest:
            calls = ctx.evidence.filter(pl.col("number") == n)
            counts = [int((calls["call"] == k).sum()) for k in figures.CALLS]
            rows.append((names[n], *(f"{c:,}" for c in counts),
                         f"{calls['padj'].is_not_null().sum():,}"))  # fmt: skip
        items.append(ResultTable(fill("figures.evidence.rest", n=len(rest)),
                                 text.words(PAGE, "figures.evidence.rest_columns"), tuple(rows),
                                 numeric=(1, 2, 3)))  # fmt: skip
    group = Group(word("figures.evidence.group"), tuple(items))
    return Evidence(group, (heat, "evidence_view", words))


def _module_page(ctx: Context, drawn: list[Drawn], left: list[Card],
                 evidence: Evidence | None = None) -> ModulePage:  # fmt: skip
    """The page over the cards it keeps: the first drawn whole, the rest packed (C3.1).

    The list from a link comes first, then the evidence view, then the cards: a card's group
    number counts the groups before it.
    """
    results: list[Result] = [] if evidence is None else [evidence.group]
    if (listed := _list_panel(ctx, drawn, len(results) + 2)) is not None:
        results.insert(0, listed)
    results += [_group(d, packed=i > 0) for i, d in enumerate(drawn)]
    if left:
        results.append(_left_out(ctx, left))
    results.append(_sources(ctx))
    gtex_numbers = ctx.provenance["gtex"]
    tissues = next((d.card.tissues for d in drawn if d.card.tissues is not None), None)
    keys = (KeyNumber(word("keys.genes"), f"{len(drawn):,}",
                      fill("keys.genes_listed", listed=len(ctx.chosen)) if left else ""),
            KeyNumber(word("keys.samples"), gtex_numbers.get("samples", word("none"))),
            KeyNumber(word("keys.tissues"), word("none") if tissues is None else
                      f"{tissues['tissue'].n_unique():,}"),
            KeyNumber(word("keys.hpa"), ctx.options.hpa_release))  # fmt: skip
    usage, budget = runtime.measure(ctx.settings.started), runtime.budget(PAGE)
    return ModulePage(
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


def fit(pages: Callable[[int], tuple[str, list[int]]], cards: int, cap: int = MAX_BYTES) -> int:
    """How many cards, in list order, a page holds under ``cap`` bytes: measured, never guessed.

    ``pages(k)`` renders the page over the first ``k`` cards and gives its HTML and the bytes of
    each packed card, in order. Cards leave from the end, the most packed bytes first needed to
    cover the excess, and the page is measured again, since the lines naming the cards left out
    grow as they leave. The first card, drawn whole, always stays.
    """
    kept = cards
    while kept > 1:
        html, packed = pages(kept)
        over = len(html.encode()) - cap
        if over <= 0:
            break
        freed = 0
        while kept > 1 and freed < over:
            kept -= 1
            freed += packed[kept - 1] if kept - 1 < len(packed) else 0
    return kept


def build(writer: Writer, ctx: Context) -> list[Path]:
    """Write the kept cards' tables and figures, then the page; the page first in the list.

    The card list is packed to the page cap (M9.8): every card is laid out without its
    files, the page is measured, and the cards past what fits are named with the command
    that builds them.
    """
    page = writer.out / f"{stem(PAGE, f'{writer.dataset}_overview', 1, day=writer.day)}.html"
    drawn = [_card(writer, c, ctx) for c in ctx.chosen]
    view = evidence_view(writer, ctx)

    def pages(k: int) -> tuple[str, list[int]]:
        """The page over the first ``k`` cards, as written, and each packed card's bytes."""
        html = rendered(_module_page(ctx, drawn[:k], [d.card for d in drawn[k:]], view), page)
        return html, packed_bytes(html)

    kept = fit(pages, len(drawn), MAX_BYTES)
    numbers = [d.card.number for d in drawn[:kept]]
    shown = replace(
        ctx,
        chosen=[d.card for d in drawn[:kept]],
        evidence=ctx.evidence.filter(pl.col("number").is_in(numbers)),
    )
    _tables(writer, shown)
    if shown.own is not None:
        boxes.write_tables(writer, shown.own, shown.spec, shown.chosen)
    figure_tables.write(writer, ctx.chosen, ctx.evidence, shown.chosen, ctx.family,
                        ctx.options.hpa_release)  # fmt: skip
    figures_kept = [f for d in drawn[:kept] for f in d.figures]
    shown_view = [(*view.figure, None)] if view is not None else []
    for data, what, words, controls in shown_view + figures_kept:
        writer.files += write_figure(data, writer.out, writer.stem(what), title=words.title,
                                     data_releases=writer.releases, text=words)  # fmt: skip
        if controls is not None:
            # Opened alone, the figure keeps the control its contract does not carry (E.1d).
            block = replace(figure_block(data, writer.stem(what), words), controls=controls)
            write(FigurePage(words.title, block), writer.out / f"{writer.stem(what)}.html")
    # Laid out again once the files are written: a block finds its Reproduce panel by its stem,
    # and the panel exists only after write_figure found the figure's spec (F.1).
    again = [_card(writer, d.card, ctx) for d in drawn[:kept]]
    built = _module_page(ctx, again, [d.card for d in drawn[kept:]], evidence_view(writer, ctx))
    return [write(built, page), write_entry(writer.out, "genecard"), *writer.files]
