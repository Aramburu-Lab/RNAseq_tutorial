"""Which genes get a card, who each is, and its rows in the dataset's own M5 results.

Genes are named through ``ctx.report.genes``, the reader the paste box and ``--genes-from``
already share, against gather's genes table: a symbol naming two genes, or none, stops the
run with every such name listed. With no gene named, the list is the dataset file's
registered box genes, then the featured contrast's top up calls and top down calls, each
ordered by padj, then absolute log2 FC, then Ensembl number.

A symbol absent from the table is looked up by name in HPA's search (``Renamer``): HPA lists a
gene's symbol and its synonyms, previous symbols among them, so NDUFA4L2, renamed COXFA4L2 by
GENCODE 50, still names ENSG00000185633, and the card states the rename.
"""

import re
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

from ctx.dataset import Dataset, DatasetError
from ctx.de.results import newest
from ctx.genecard import GenecardError, annotation, hpa
from ctx.genecard.cache import Cache
from ctx.io.naming import MODULES
from ctx.report import genes

# Calls of each sign the default list takes from the featured contrast.
TOP = 25
EVIDENCE = ("log2FoldChange", "log2FoldChange_shrunk", "padj", "call")


@dataclass
class Card:
    """Everything one gene's card draws; a source's value is None with a note when unread."""

    number: int
    gene_id: str
    symbol: str
    reason: str
    gtex_id: str = ""
    tissues: pl.DataFrame | None = None
    exons: pl.DataFrame | None = None
    exons_dropped: int = 0
    xml: hpa.Xml | None = None
    facts: hpa.Facts | None = None
    gtf: annotation.Mane | None = None
    rest: annotation.Mane | None = None
    lines: hpa.Lines | None = None
    renamed: str = ""
    junctions: pl.DataFrame | None = None
    junction_id: str = ""
    junction_tpm: dict[str, float] | None = None
    notes: dict[str, str] = field(default_factory=dict)


def _out(root: Path, module: str) -> Path:
    """A module's output folder under the root."""
    return root / "outputs" / f"{dict(MODULES)[module]:02d}_{module}"


def gene_table(root: Path, dataset: str) -> tuple[pl.DataFrame, genes.GeneTable]:
    """Gather's genes table and the lookup built on it."""
    frame = pl.read_parquet(newest(_out(root, "gather"), dataset, "genes"))
    return frame, genes.table(frame["gene_id"].to_list(), frame["gene_name"].to_list())


@dataclass
class Renamer:
    """Finds the gene of a symbol the table lacks by HPA's symbol and synonym fields.

    Attributes:
        cache: The response cache the search goes through.
        release: The HPA release.
        table: Gather's genes; a gene HPA names must be one of them.
        renamed: Each resolved gene's number to the name typed for it.
        why: Each unresolved name to why its search gave nothing, where the search failed.
    """

    cache: Cache
    release: str
    table: genes.GeneTable
    renamed: dict[int, str] = field(default_factory=dict)
    why: dict[str, str] = field(default_factory=dict)

    def __call__(self, word: str) -> list[int]:
        """The table's genes HPA gives ``word`` to as its symbol or a synonym, exactly."""
        try:
            genes.number(word)
            return []  # an Ensembl ID the table lacks has no other name to try
        except ValueError:
            pass
        key = re.sub(r"[^A-Z0-9._-]", "_", word.upper())
        url = hpa.search_url(self.release, word, ("g", "gs", "eg"))
        answer = self.cache.get("hpa", self.release, "search", f"name_{key}", url)
        if answer.body is None:
            self.why[word] = f"HPA search: {answer.error}"
            return []
        found = []
        for gene, symbol, synonyms in hpa.parse_names(answer.body):
            # HPA's search matches any field; only the symbol or a synonym names the gene.
            if word.upper() in {symbol.upper(), *(s.upper() for s in synonyms)}:
                number = genes.number(gene)
                if number in self.table.row:
                    found.append(number)
        if len(found) == 1:
            self.renamed[found[0]] = word
        return found


def named(text: str, table: genes.GeneTable, source: str,
          resolve: Callable[[str], list[int]] | None = None) -> list[int]:  # fmt: skip
    """The genes a text names, refused whole if any name is ambiguous or unknown.

    ``resolve`` is asked for each name the table lacks; its single gene follows the genes
    found in the table, several make the name ambiguous, none leaves it unknown.
    """
    found = genes.parse(text, table)
    problems = [f"{word} names {', '.join(ids)}" for word, ids in found.ambiguous]
    numbers = list(found.found)
    for word in found.missing:
        other = resolve(word) if resolve is not None else []
        if len(other) == 1:
            numbers.append(other[0])
        elif other:
            problems.append(f"{word} names {', '.join(genes.ensembl(n) for n in other)}")
        else:
            problems.append(f"{word} names no gene")
    if problems:
        raise GenecardError(f"{source}: " + "; ".join(problems))
    return list(dict.fromkeys(numbers))


def from_file(path: Path, table: genes.GeneTable,
              resolve: Callable[[str], list[int]] | None = None) -> list[int]:  # fmt: skip
    """The genes of ``--genes-from``: an M5 results table by ``gene_name``, or a GS.1 list."""
    if path.suffix in {".parquet", ".csv"}:
        frame = pl.read_parquet(path) if path.suffix == ".parquet" else pl.read_csv(path)
        if "gene_id" in frame.columns:
            return named("\n".join(frame["gene_id"].to_list()), table, path.name)
        if "gene_name" not in frame.columns:
            raise GenecardError(f"{path.name}: neither gene_id nor gene_name, as M5 writes")
        names = "\n".join(frame["gene_name"].drop_nulls().to_list())
        return named(names, table, path.name, resolve)
    listed = genes.read_list(path, table).genes
    if listed.ambiguous or listed.missing:
        return named(path.read_text(encoding="utf-8-sig"), table, path.name, resolve)
    return list(listed.found)


def results(root: Path, dataset: str, contrast: str) -> pl.DataFrame:
    """One contrast's M5 results, keyed by unversioned number."""
    frame = pl.read_parquet(newest(_out(root, "de"), dataset, f"{contrast}_results"))
    numbers = [genes.number(g) for g in frame["gene_id"].to_list()]
    return frame.select(*EVIDENCE).with_columns(number=pl.Series(numbers, dtype=pl.Int64))


def ranked(table: pl.DataFrame, call: str, top: int = TOP) -> list[int]:
    """The top ``call`` genes: lowest padj, then largest absolute log2 FC, then number."""
    chosen = table.filter(pl.col("call") == call).sort(
        [pl.col("padj"), pl.col("log2FoldChange").abs(), pl.col("number")],
        descending=[False, True, False],
    )
    return chosen["number"].head(top).to_list()


def card_list(spec: Dataset, table: genes.GeneTable, requested: list[int] | None,
              contrast: str | None, root: Path,
              resolve: Callable[[str], list[int]] | None = None
              ) -> list[tuple[int, str]]:  # fmt: skip
    """Each card's gene and why it has one, in card order, each gene once.

    Raises:
        DatasetError: With no gene named, if the dataset file has no ``[de]`` featured
            contrast to take the top calls from.
    """
    if requested is not None:
        return [(n, "named") for n in dict.fromkeys(requested)]
    if spec.de is None:
        raise DatasetError(f"{spec.path.name}: no [de]; name genes with --genes or --genes-from")
    featured = contrast or spec.de.featured
    boxed = named("\n".join(spec.de.box_genes), table, f"{spec.path.name} box_genes", resolve)
    listed = [(n, "registered") for n in boxed]
    top = results(root, spec.slug, featured)
    listed += [(n, f"up in {featured}") for n in ranked(top, "up")]
    listed += [(n, f"down in {featured}") for n in ranked(top, "down")]
    first: dict[int, str] = {}
    for n, why in listed:
        first.setdefault(n, why)
    return list(first.items())


def evidence(root: Path, spec: Dataset, numbers: list[int]) -> pl.DataFrame:
    """Every listed gene's M5 row in every contrast; a gene a contrast never tested has nulls.

    Raises:
        DatasetError: If a contrast of the dataset file has no results table.
    """
    frames = []
    for contrast in spec.contrasts:
        table = results(root, spec.slug, contrast.name).filter(pl.col("number").is_in(numbers))
        frames.append(
            pl.DataFrame({"number": numbers}, schema={"number": pl.Int64})
            .join(table, on="number", how="left", maintain_order="left")
            .with_columns(contrast=pl.lit(contrast.name))
        )
    if not frames:
        raise DatasetError(f"{spec.path.name}: no contrasts, so no own evidence")
    return pl.concat(frames).select("number", "contrast", *EVIDENCE)
