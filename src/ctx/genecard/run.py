"""``ctx genecard``: one page of gene cards for a dataset, its tables and figures.

Genes come from ``cards``; GTEx V11 from the files ``ctx setup`` registered (``gtex``); HPA,
Ensembl and GTEx V10's junctions (``junctions``) from the response cache (``cache``); MANE
exons from the dataset file's GTF (``annotation``). A source that cannot be read leaves its
panels unavailable, naming why, and the rest of the card stands. ``page`` draws the result.
"""

import json
from collections.abc import Callable
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any

import polars as pl

from ctx import runtime
from ctx.dataset import Dataset, DatasetError, load
from ctx.de.run import Writer
from ctx.genecard import (
    GenecardError,
    annotation,
    cards,
    gtex,
    hpa,
    junctions,
    listing,
    page,
    samples,
)
from ctx.genecard.cache import Cache
from ctx.io.naming import MODULES
from ctx.report import genes
from ctx.report.builder import ParameterCard
from ctx.setup import MANIFEST as SETUP_MANIFEST

MODULE = "genecard"
GTF = "annotation_gtf"


@dataclass(frozen=True)
class Options:
    """Every option of ``ctx genecard`` but the dataset.

    Attributes:
        dataset_file: The dataset file: contrasts, box genes, GTF and GTEx resources.
        genes: Gene symbols or Ensembl IDs typed on the command line.
        genes_from: An M5 results table or a saved gene list.
        contrast: The contrast whose top calls the default list takes; None, the featured.
        hpa_release: The HPA release every HPA request names, as ``v<NN>``.
        offline: Read the cache only.
        refresh: Move each gene's cached responses aside and fetch them again.
    """

    dataset_file: Path
    genes: tuple[str, ...] = ()
    genes_from: Path | None = None
    contrast: str | None = None
    hpa_release: str = hpa.RELEASE
    offline: bool = False
    refresh: bool = False


def _registered(root: Path, spec: Dataset, name: str) -> Path | str:
    """A resource's file, if ``ctx setup`` recorded it with its declared hash; else why not."""
    try:
        resource = spec.resource(name)
    except StopIteration:
        return f"{spec.path.name} registers no {name}"
    manifest_path = root / SETUP_MANIFEST
    manifest = json.loads(manifest_path.read_text()) if manifest_path.is_file() else {}
    entry, path = manifest.get(resource.dest), root / resource.dest
    if not (entry and entry["sha256"] == resource.sha256 and path.is_file()):
        return f"{resource.dest} is not registered: ctx setup --dataset-file {spec.path.name}"
    if path.stat().st_size != resource.bytes:
        return f"{resource.dest}: {path.stat().st_size} B, registered as {resource.bytes} B"
    return path


def _gtex(root: Path, spec: Dataset, chosen: list[cards.Card]) -> dict[str, str]:
    """Fill each card's GTEx tissue spread and exon shares; return GTEx's provenance."""
    found = {k: _registered(root, spec, v) for k, v in gtex.RESOURCES.items()}
    if missing := [str(v) for v in found.values() if isinstance(v, str)]:
        for card in chosen:
            card.notes["gtex"] = "; ".join(missing)
        return {}
    files = gtex.Files(*(Path(found[k]) for k in gtex.RESOURCES))
    who = gtex.donors(files)
    samples, rows = gtex.read_tpm(files.tpm, [c.number for c in chosen])
    exons = gtex.read_exons(files.exons, list(rows), samples)
    for card in chosen:
        if (row := rows.get(card.number)) is None:
            card.notes["gtex"] = f"not in GTEx {gtex.RELEASE} (GENCODE {gtex.GENCODE})"
            continue
        card.gtex_id, card.tissues = row.gene_id, gtex.tissue_spread(row.tpm, samples, who)
        if card.number in exons:
            card.exons, card.exons_dropped = gtex.exon_shares(exons[card.number], samples, who)
    manifest = json.loads((root / SETUP_MANIFEST).read_text())
    return {"retrieved": manifest[spec.resource(gtex.RESOURCES["tpm"]).dest]["retrieved"],
            "samples": f"{len(samples):,}", "donors": f"{who['sample'].len():,}"}  # fmt: skip


def _ensembl_release(cache: Cache) -> str:
    """Why Ensembl REST cannot be read for the pinned release, or "" when it can."""
    answer = cache.get("ensembl", annotation.ENSEMBL_RELEASE, "rest", "info_data",
                       annotation.release_url())  # fmt: skip
    if answer.body is None:
        return answer.error
    served = annotation.rest_release(answer.body)
    if annotation.ENSEMBL_RELEASE not in served.split(","):
        return f"Ensembl REST serves release {served}, not {annotation.ENSEMBL_RELEASE}"
    return ""


def _listed(root: Path, spec: Dataset, options: Options, table: genes.GeneTable,
            renamer: cards.Renamer) -> list[tuple[int, str]]:  # fmt: skip
    """The cards' genes and reasons; a name HPA could not be asked about says why.

    Raises:
        GenecardError: If a name is ambiguous or names no gene.
    """
    try:
        named = None
        if options.genes_from is not None:
            named = cards.from_file(options.genes_from, table, renamer)
        if options.genes:
            named = (named or []) + cards.named("\n".join(options.genes), table, "--genes",
                                                renamer)  # fmt: skip
        return cards.card_list(spec, table, named, options.contrast, root, renamer)
    except GenecardError as e:
        if not renamer.why:
            raise
        why = "; ".join(f"{word}: {reason}" for word, reason in renamer.why.items())
        raise GenecardError(f"{e}; {why}") from e


def _lines(cache: Cache, release: str, card: cards.Card, lines: hpa.CellLines) -> str:
    """Fill the card's cell-line values from HPA's search; return the retrieval, or ""."""
    gene = genes.ensembl(card.number)
    url = hpa.search_url(release, gene, ("g", "eg", lines.group))
    answer = cache.get("hpa", release, gene, f"lines_{lines.group}", url)
    if answer.body is None:
        card.notes["lines"] = answer.error
        return ""
    card.lines = hpa.parse_lines(answer.body, gene, lines.group)
    return answer.retrieved


def _junctions(cache: Cache, chosen: list[cards.Card]) -> str:
    """Fill each card's GTEx V10 junction medians from the API; return the newest retrieval."""
    newest = ""
    for card in chosen:
        gene = genes.ensembl(card.number)
        record = cache.get("gtex", junctions.DATASET, gene, f"reference_{junctions.GENCODE}",
                           junctions.gene_url(gene))  # fmt: skip
        if record.body is None:
            card.notes["junctions"] = record.error
            continue
        gencode_id = junctions.versioned(record.body, gene)
        if gencode_id is None:
            card.notes["junctions"] = (
                f"not in GENCODE {junctions.GENCODE}, GTEx {junctions.RELEASE}'s annotation"
            )
            continue
        answer = cache.get("gtex", junctions.DATASET, gene, "junctions",
                           junctions.junction_url(gencode_id))  # fmt: skip
        if answer.body is None:
            card.notes["junctions"] = answer.error
            continue
        newest = max(newest, record.retrieved, answer.retrieved)
        table = junctions.parse(answer.body, gencode_id)
        if table.is_empty():
            card.notes["junctions"] = f"no junction of {gencode_id} in {junctions.DATASET}"
            continue
        card.junction_id, card.junctions = gencode_id, table
    return newest


def _remote(cache: Cache, release: str, chosen: list[cards.Card],
            lines: hpa.CellLines | None) -> dict[str, str]:  # fmt: skip
    """Fill each card's HPA and Ensembl panels from the cache; return the newest retrieval."""
    newest: dict[str, str] = {}
    rest_note = _ensembl_release(cache)
    for card in chosen:
        gene = genes.ensembl(card.number)
        answers = {e: cache.get("hpa", release, gene, e, hpa.url(release, gene, e))
                   for e in hpa.ENDPOINTS}  # fmt: skip
        if failed := [f"{e}: {a.error}" for e, a in answers.items() if a.body is None]:
            card.notes["hpa"] = "; ".join(failed)
        else:
            card.xml = hpa.parse_xml(answers["xml"].body, release)  # type: ignore[arg-type]
            card.facts = hpa.parse_json(answers["json"].body)  # type: ignore[arg-type]
            seen = [a.retrieved for a in answers.values()]
            if lines is not None:
                seen.append(_lines(cache, release, card, lines))
            newest["hpa"] = max(newest.get("hpa", ""), *seen)
        if rest_note:
            card.notes["ensembl"] = rest_note
            continue
        lookup = cache.get("ensembl", annotation.ENSEMBL_RELEASE, gene, "lookup",
                           annotation.lookup_url(gene))  # fmt: skip
        if lookup.body is None:
            card.notes["ensembl"] = lookup.error
        else:
            card.rest = annotation.rest_mane(lookup.body)
            newest["ensembl"] = max(newest.get("ensembl", ""), lookup.retrieved)
    return newest


def _own(root: Path, spec: Dataset, chosen: list[cards.Card], evidence: pl.DataFrame
         ) -> tuple[samples.Values | None, dict[str, samples.Test]]:  # fmt: skip
    """The listed genes' M5 normalised values and the box's family, BH over the card list.

    The family is the listed genes M5 tested in the featured contrast that a test can rank.
    Without M5's size factors every card names why, and draws no own-sample panel.
    """
    own = samples.read(root, spec, [c.gene_id for c in chosen])
    if isinstance(own, str):
        for card in chosen:
            card.notes["own"] = own
        return None, {}
    if spec.de is None:
        return own, {}
    featured = spec.de.featured
    tested = set(evidence.filter((pl.col("contrast") == featured)
                                 & pl.col("padj").is_not_null())["number"].to_list())  # fmt: skip
    return own, samples.family(own, featured, [c.gene_id for c in chosen if c.number in tested])


def run(settings: runtime.Settings, cards_: tuple[ParameterCard, ...], command: str,
        dataset: str, options: Options, opener: Callable[[str], Any] | None = None
        ) -> list[Path]:  # fmt: skip
    """Build the gene cards page, its tables and figures under ``outputs/11_genecard``.

    Args:
        settings: The resolved runtime.
        cards_: Parameter cards for every option of this command.
        command: The command line, as the reader should retype it.
        dataset: The dataset's slug.
        options: The command's options.
        opener: How the cache opens a URL; None is urllib. Tests pass their own.

    Returns:
        The page, then every table and figure file.

    Raises:
        DatasetError: If the dataset file names another slug, or an M5 table is missing.
        GenecardError: If a gene is named ambiguously or not at all, a file breaks its
            format, or a cached response changed.
    """
    spec = load(options.dataset_file)
    if spec.slug != dataset:
        raise DatasetError(f"{options.dataset_file} is for {spec.slug!r}, not {dataset!r}")
    if not hpa.RELEASE_PATTERN.fullmatch(options.hpa_release):
        raise GenecardError(f"--hpa-release {options.hpa_release!r}: a major release, as v25")
    root = settings.root
    lines = hpa.cell_lines(spec.info, spec.path.name)
    cache = Cache(root, offline=options.offline, refresh=options.refresh)
    if opener is not None:
        cache.opener = opener
    frame, table = cards.gene_table(root, dataset)
    renamer = cards.Renamer(cache, options.hpa_release, table)
    listed = _listed(root, spec, options, table, renamer)
    ids = frame["gene_id"].to_list()
    chosen = [cards.Card(n, ids[table.row[n]], table.symbols[table.row[n]] or genes.ensembl(n), why,
                         renamed=renamer.renamed.get(n, ""))
              for n, why in listed]  # fmt: skip
    evidence = cards.evidence(root, spec, [c.number for c in chosen])
    tested = listing.evidence({c.name: cards.results(root, dataset, c.name)
                               for c in spec.contrasts}, table)  # fmt: skip
    own, family = _own(root, spec, chosen, evidence)
    provenance = {"gtex": _gtex(root, spec, chosen)}
    gtf = _registered(root, spec, GTF)
    mane = annotation.gtf_mane(gtf, [c.number for c in chosen]) if isinstance(gtf, Path) else {}
    for card in chosen:
        card.gtf = mane.get(card.number)
    provenance["remote"] = _remote(cache, options.hpa_release, chosen, lines)
    provenance["remote"]["junctions"] = _junctions(cache, chosen)
    out = root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}"
    releases = {"GTEx": f"{gtex.RELEASE} ({gtex.RELEASE_DATE}), GENCODE {gtex.GENCODE}",
                "HPA": options.hpa_release, "Ensembl": annotation.ENSEMBL_RELEASE,
                "dataset file": str(spec.path)}  # fmt: skip
    writer = Writer(out, dataset, date.today(), releases)
    context = page.Context(spec=spec, options=options, chosen=chosen, evidence=evidence,
                           provenance=provenance, settings=settings, command=command,
                           cards=cards_, tested=tested, own=own, family=family)  # fmt: skip
    return page.build(writer, context)
