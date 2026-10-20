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
# The sets a family draws when the reader asks for the smallest p rather than the passing ones,
# and the count a family with fewer passing sets is filled to in grey (figure controls plan:
# "7 smallest FDR per direction"; PK.3 took this cap rather than a new number).
TOP = 7
HOW = ("passing", "top")


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


def drawn(family: pl.DataFrame, how: str = "passing", alpha: float = ALPHA,
          top: int = TOP) -> pl.DataFrame:  # fmt: skip
    """The sets one family's bars draw, in order, each marked by whether it passes (PK.3).

    One definition, ``ctxOra.drawn`` its twin. The order: drop the sets holding none of the
    list (p = 1 says nothing of them), sort by p, a tie by set name, mark FDR < ``alpha``;
    then ``passing`` keeps every set that passes and, where fewer than ``top`` do, the
    failing sets of smallest p up to ``top`` in all; ``top`` keeps the first ``top``.
    Rytis's 6 genes on Step-1 passed none, and the panel drew an empty chart.

    Args:
        family: One family's ORA rows, with ``gene_set``, ``hits``, ``p`` and ``fdr``.
        how: ``passing`` or ``top``, the panel's "Sets drawn".
        alpha: The FDR a set must stay under to pass.
        top: The count a short family is filled to, and the count ``top`` draws.

    Returns:
        The rows drawn, with ``passes``, passing sets first.

    Raises:
        ValueError: If ``how`` is neither.
    """
    if how not in HOW:
        raise ValueError(f"sets drawn: {how!r}, not one of {HOW}")
    held = family.filter(pl.col("hits") > 0).sort("p", "gene_set")
    held = held.with_columns(passes=pl.col("fdr") < alpha)
    if how == "top":
        return held.head(top)
    passing = held.filter(pl.col("passes"))
    if passing.height >= top:
        return passing
    return pl.concat([passing, held.filter(~pl.col("passes")).head(top - passing.height)])


def bars(frame: pl.DataFrame, how: str = "passing", alpha: float = ALPHA,
         top: int = TOP) -> pl.DataFrame:  # fmt: skip
    """Every family's drawn sets as the bars stand, ``ctxOra.blocks`` its twin (PK.3).

    Each family by ``drawn``; the family whose first set has the smallest p first, a tie in
    the frame's family order; ``series`` is the family for a set that passes and
    ``<family>_failing`` for one that does not, the grey bars.

    Returns:
        The rows drawn, with ``passes`` and ``series``; empty when no set holds a gene.
    """
    families = [drawn(frame.filter(pl.col("direction") == d), how, alpha, top)
                for d in frame["direction"].unique(maintain_order=True)]  # fmt: skip
    families = sorted((f for f in families if f.height), key=lambda f: f["p"][0])
    if not families:
        return frame.clear().with_columns(passes=pl.lit(True), series=pl.lit(""))
    return pl.concat(families).with_columns(
        series=pl.when("passes").then("direction").otherwise(pl.col("direction") + "_failing")
    )


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
