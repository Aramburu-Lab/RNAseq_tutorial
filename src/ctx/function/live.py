"""M6's live over-representation panel (GS.5): what the page embeds for ``ctxOra``.

The page reruns ``lite.ora`` on the reader's list (``templates/_ora.js.j2``, tested against it
under node). It embeds one index of gene symbols, each contrast's tested symbols as a bitset
over it, the gene sets as index lists, the featured contrast's calls as the default list, and
the Ensembl numbers a link's list names, each with its symbol. With no list from a link, the
page's ORA is the build's: the featured up and down calls against its tested symbols.
"""

import base64
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import polars as pl

from ctx.de.stats import ALPHA
from ctx.function import lite
from ctx.report import genes, text

PAGE = "function"
# The sets drawn per family when the reader asks for the smallest FDR rather than the passing
# ones (figure controls plan: "7 smallest FDR per direction").
TOP = 7


def bitset(members: Sequence[int], size: int) -> str:
    """Index positions as a bitset, least significant bit first, base64: 17,830 genes in 3 kB."""
    bits = bytearray((size + 7) // 8)
    for i in members:
        bits[i // 8] |= 1 << (i % 8)
    return base64.b64encode(bytes(bits)).decode("ascii")


def members(encoded: str, size: int) -> list[int]:
    """The positions a ``bitset`` holds: its inverse, as the page reads it."""
    bits = base64.b64decode(encoded)
    return [i for i in range(size) if bits[i // 8] >> (i % 8) & 1]


def symbols(results: pl.DataFrame, call: str | None = None) -> list[str]:
    """Symbols of the tested genes, or of the genes with one call, as ``run._symbols`` takes."""
    tested = set(lite.background(results["gene_id"].to_list(), results["padj"].to_list()))
    rows = results.filter(pl.col("gene_id").is_in(list(tested)))
    rows = rows.filter(pl.col("call") == call) if call else rows
    return sorted({s for s in rows["gene_name"].drop_nulls() if s})


def by_contrast(results: pl.DataFrame, sets: Mapping[str, Sequence[str]]) -> pl.DataFrame:
    """ORA of one contrast's own calls, the live panel's contrast select (WO7, V2).

    The one definition the page's choice follows, ``ctxOra.contrast`` its twin: the contrast's
    up and down calls, each its own family, against its own tested symbols. It does not move
    with the GSEA bars' contrast view, a lookup of its own (D-073); each select is its own
    state.

    Returns:
        One row per kept set and direction, in ``run.ORA_COLUMNS``'s order.
    """
    tested = symbols(results)
    rows = [(d, r.name, r.size, r.hits, r.expected, r.p, r.fdr, ";".join(r.genes))
            for d in ("up", "down")
            for r in lite.ora(symbols(results, d), tested, sets)]  # fmt: skip
    schema = ["direction", "gene_set", "size", "hits", "expected", "p", "fdr", "genes"]
    return pl.DataFrame(rows, schema=schema, orient="row")


def listed(path: Path, results: pl.DataFrame, sets: Mapping[str, Sequence[str]]) -> pl.DataFrame:
    """ORA of a saved list, as the page's live panel tests one from a link (GS.5).

    Its genes are signed by their raw log2 FC in ``results``, the featured contrast, as a
    link's are, and each direction is its own BH family against the tested symbols.

    Returns:
        One row per kept set and direction, in ``run.ORA_COLUMNS``'s order.
    """
    tested = results.filter(pl.col("padj").is_not_null() & pl.col("gene_name").is_not_null())
    table = genes.table(tested["gene_id"].to_list(), tested["gene_name"].to_list(),
                        tested["log2FoldChange"].to_list())  # fmt: skip
    found = genes.read_list(path, table).genes.found
    names = tested["gene_name"].to_list()
    hits: dict[str, list[str]] = {"up": [], "down": []}
    for n in found:
        at = table.row[n]
        hits["down" if (table.lfc or ())[at] < 0 else "up"].append(names[at])
    rows = [(d, r.name, r.size, r.hits, r.expected, r.p, r.fdr, ";".join(r.genes))
            for d, h in hits.items() for r in lite.ora(h, symbols(results), sets)]  # fmt: skip
    schema = ["direction", "gene_set", "size", "hits", "expected", "p", "fdr", "genes"]
    return pl.DataFrame(rows, schema=schema, orient="row")


def panel(
    results: Mapping[str, pl.DataFrame],
    featured: str,
    sets: Mapping[str, Sequence[str]],
    labels: Mapping[str, str],
    words: Mapping[str, str],
) -> dict[str, Any]:
    """The live ORA's settings and data, the build's default list first.

    Args:
        results: Each contrast's M5 results, in page order.
        featured: The contrast whose calls are the default list and whose tested genes are the
            default background.
        sets: The gene sets, in file order.
        labels: Each contrast's label, and ``union`` for the trap background.
        words: Each gene set's name in words, as the bars name it.

    Returns:
        JSON-ready settings for the page script.
    """
    tested = {c: symbols(r) for c, r in results.items()}
    union = sorted({g for gs in sets.values() for g in gs} | set(tested[featured]))
    every = sorted(set().union(*tested.values(), union))
    at = {s: i for i, s in enumerate(every)}
    backgrounds = [[c, labels[c], bitset([at[s] for s in t], len(every)), len(t)]
                   for c, t in tested.items()]  # fmt: skip
    backgrounds.append(["union", labels["union"], bitset([at[s] for s in union], len(every)),
                        len(union)])  # fmt: skip
    numbers: dict[int, int] = {}
    for r in results.values():
        for gene_id, name in r.select("gene_id", "gene_name").iter_rows():
            if not (name and name in at):
                continue
            try:
                numbers.setdefault(genes.number(gene_id), at[name])
            except ValueError:
                continue  # an ID no link can carry (not Ensembl) is reached by its symbol
    calls = {k: [at[s] for s in symbols(results[featured], k)] for k in ("up", "down")}
    # Every contrast's own calls, each tested against its own tested genes (WO7): the build's
    # first, as the select lists them.
    contrasts = [[c, labels[c], *([at[s] for s in symbols(results[c], k)] for k in ("up", "down"))]
                 for c in [featured, *(n for n in results if n != featured)]]  # fmt: skip
    return {
        "kind": "ora", "featured": featured, "symbols": every, "labels": dict(words),
        "backgrounds": backgrounds, "contrasts": contrasts, "top": TOP,
        "sets": [[name, [at[g] for g in gs if g in at]] for name, gs in sets.items()],
        "numbers": sorted(numbers.items()), "calls": calls,
        "min_size": lite.MIN_SIZE, "max_size": lite.MAX_SIZE, "alpha": ALPHA,
        "named": genes.PAGE_ROWS, "words": text.raw(PAGE, "live"),
        "errors": text.raw("ui", "select"),
    }  # fmt: skip
