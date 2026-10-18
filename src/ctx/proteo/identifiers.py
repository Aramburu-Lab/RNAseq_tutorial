"""The one way a protein table's identifiers become the RNA side's genes, with what did not map.

A row's identifier cell may name a protein group: several identifiers joined by ``;``. Each is
first reduced to a gene key of its kind (a symbol as written, a UniProt accession through the
dataset's map with its isoform suffix and FASTA wrapping removed, an Ensembl gene ID without its
version), then the keys meet gather's genes table. A row is mapped only when its group names one
gene and that gene is one row of the RNA side; a group naming several genes is never resolved
to one of them, and is listed instead (CONVENTIONS.md, proteomics; the R3 brief).
"""

import re
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import polars as pl

ID_KINDS = ("symbol", "uniprot", "ensembl")
# Why a row is not mapped, in the order the checks run.
NO_ID, NO_MAP, SEVERAL, ABSENT, SHARED = "no_id", "no_map", "several", "absent", "shared"
REASONS = (NO_ID, NO_MAP, SEVERAL, ABSENT, SHARED)
# MaxQuant and DIA-NN join a group's members with a semicolon.
GROUP_SEPARATOR = ";"
# A FASTA header kept whole by the search: sp|P04637|P53_HUMAN or tr|...; the accession is the
# second field.
FASTA_ID = re.compile(r"^(?:sp|tr)\|([^|]+)\|")
# A UniProt isoform: P04637-2 is P04637's gene.
ISOFORM = re.compile(r"-\d+$")
# An Ensembl gene ID's version: ENSG00000141510.18 is ENSG00000141510. GENCODE's _PAR_Y copy
# keeps its suffix, so it stays apart from the X copy it would otherwise collide with.
VERSION = re.compile(r"\.\d+$")


@dataclass(frozen=True)
class Mapped:
    """Each row's RNA gene, or why it has none.

    Attributes:
        genes: Per row, the RNA gene's name as pages show it, or None.
        reasons: Per row, ``mapped`` or one of ``REASONS``.
        several: (identifier cell, the genes it names) of every row naming several genes.
    """

    genes: tuple[str | None, ...]
    reasons: tuple[str, ...]
    several: tuple[tuple[str, str], ...]

    def counts(self) -> dict[str, int]:
        """Rows per reason, every reason present, ``mapped`` included."""
        found = Counter(self.reasons)
        return {r: found.get(r, 0) for r in ("mapped", *REASONS)}


def rna_names(genes: pl.DataFrame) -> list[str]:
    """Each RNA row's name as pages show it, so two rows never share a name.

    The name is the row's symbol where no other row shares it, else its unversioned gene ID.
    """
    ids = [VERSION.sub("", g) for g in genes["gene_id"].to_list()]
    symbols = [n or "" for n in genes["gene_name"].to_list()]
    shared = {n for n, k in Counter(symbols).items() if k > 1}
    return [n if n and n not in shared else i for i, n in zip(ids, symbols, strict=True)]


def key(identifier: str, kind: str, id_map: Mapping[str, str]) -> str | None:
    """One identifier's gene key: the symbol, the mapped accession's symbol, or the bare ID."""
    text = identifier.strip()
    if kind == "uniprot":
        fasta = FASTA_ID.match(text)
        accession = ISOFORM.sub("", fasta.group(1) if fasta else text)
        return id_map.get(accession) or None
    return VERSION.sub("", text) if kind == "ensembl" else text


def _gene(found: str, index: Mapping[str, list[int]], names: list[str]) -> tuple[str | None, str]:
    """One gene key's RNA gene and ``mapped``, or None and why not."""
    rows = index.get(found, [])
    if len(rows) != 1:
        return None, ABSENT if not rows else SHARED
    return names[rows[0]], "mapped"


def to_genes(
    cells: Sequence[str | None],
    kind: str,
    genes: pl.DataFrame,
    id_map: Mapping[str, str] | None = None,
) -> Mapped:
    """Map each row's identifier cell to one RNA gene, or say why not.

    The order is fixed: split the group, key each member, refuse a member with no key, refuse
    several genes, then meet the RNA side.

    Args:
        cells: Each row's identifier cell, a group joined by ``GROUP_SEPARATOR``.
        kind: One of ``ID_KINDS``.
        genes: Gather's genes table, ``gene_id`` and ``gene_name``.
        id_map: UniProt accession to gene symbol; needed for ``uniprot`` alone.

    Returns:
        Each row's gene or reason, and the several-gene groups.

    Raises:
        ValueError: If the kind is unknown, or UniProt identifiers come without a map.
    """
    if kind not in ID_KINDS:
        raise ValueError(f"identifier kind {kind!r}: one of {ID_KINDS}")
    if kind == "uniprot" and not id_map:
        raise ValueError("UniProt identifiers need a map from accession to gene symbol")
    names = rna_names(genes)
    column = genes["gene_id"].to_list() if kind == "ensembl" else genes["gene_name"].to_list()
    index: dict[str, list[int]] = {}
    for i, value in enumerate(column):
        index.setdefault(VERSION.sub("", value) if kind == "ensembl" else value or "", []).append(i)
    out: list[str | None] = []
    reasons, several = [], []
    for cell in cells:
        members = [m for m in (cell or "").split(GROUP_SEPARATOR) if m.strip()]
        keys = [key(m, kind, id_map or {}) for m in members]
        found = sorted({k for k in keys if k})
        gene = None
        if not members:
            reason = NO_ID
        elif None in keys:
            reason = NO_MAP
        elif len(found) > 1:
            reason = SEVERAL
            several.append((cell or "", GROUP_SEPARATOR.join(found)))
        else:
            gene, reason = _gene(found[0], index, names)
        out.append(gene)
        reasons.append(reason)
    return Mapped(tuple(out), tuple(reasons), tuple(several))
