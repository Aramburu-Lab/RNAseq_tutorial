"""The door of a count matrix from any repository: read it, check it, map its gene IDs.

A genes-by-samples table of raw counts, as GEO, ArrayExpress, recount or a student's own
pipeline gives it, enters here. Every problem is collected, so one run lists them all: values
that are not raw counts, gene IDs and samples that repeat, and gene IDs the annotation cannot
place. Nothing is guessed: an ID that maps to no gene, or to several, is left out and listed
with its counts, never matched by a nearest name. ``ctx.gather.counts_run`` writes the tables.
"""

import csv
import gzip
import io
import re
from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl

from ctx.gather.curation import Problem

# The column only this route writes, so a module can tell the sheet from its columns alone.
SOURCE_COLUMN = "counts_file"
ENSEMBL = re.compile(r"(ENS[A-Z]*G\d{11})(?:\.(\d+))?")
ENTREZ = re.compile(r"\d+")
# featureCounts writes these five after its Geneid column (subread's documented output).
FEATURECOUNTS = ("Chr", "Start", "End", "Strand", "Length")
# Rows that count reads, not genes: HTSeq's five start with "__"; STAR's ReadsPerGene has four.
HTSEQ_PREFIX = "__"
STAR_ROWS = ("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
# TPM and CPM put every sample's sum at one million; raw counts put it at the library size.
PER_MILLION = 1_000_000
# A sum within 1 % of a million, in every sample, reads as TPM or CPM. A course convention:
# raw libraries of a whole study would all have to land within 1 % of one exact number.
SCALED_TOLERANCE = 0.01
# Log2 expression rarely tops 20; raw counts of a whole matrix reach thousands at least.
LOG_CEILING = 30
INT32_MAX = 2**31 - 1
ENTREZ_RESOURCE = "gencode_entrez"
TRANSCRIPT = re.compile(r'transcript_id "([^"]*)"')
GENE = re.compile(r'gene_id "([^"]*)"')


class CountsError(ValueError):
    """A matrix that cannot enter; ``problems`` holds (check, key, detail) for each."""

    def __init__(self, problems: list[Problem]) -> None:
        """Keep the problems and summarise them in the message."""
        self.problems = problems
        lines = [f"{check}: {key}: {detail}" for check, key, detail in problems]
        super().__init__(f"{len(problems)} count matrix problems:\n  " + "\n  ".join(lines))


@dataclass(frozen=True)
class Raw:
    """A matrix as read: its row IDs, sample names and values, before any check.

    Attributes:
        ids: The first column, one per gene row, summary rows removed.
        samples: The sample columns' names, in file order.
        values: Genes by samples, float64; NaN where a cell was empty.
        summary: (row name, total over samples) of each read-summary row removed.
        dropped: The featureCounts annotation columns removed, or none.
        format: "featureCounts" or "table".
        problems: Cells that are not numbers, found while reading.
    """

    ids: tuple[str, ...]
    samples: tuple[str, ...]
    values: np.ndarray
    summary: tuple[tuple[str, float], ...]
    dropped: tuple[str, ...]
    format: str
    problems: tuple[Problem, ...]


@dataclass(frozen=True)
class Mapped:
    """The matrix's rows placed on the annotation's genes.

    Attributes:
        counts: ``gene_id`` then one integer column per sample, the mapped rows only.
        genes: ``gene_id``, ``gene_name`` and ``source_id``, the matrix's own ID.
        left_out: ``source_id``, ``kind``, ``reason`` (one of ``REASONS``), ``detail`` and
            ``total``, the counts the row held over every sample.
        kinds: Rows per ID kind: ensembl, entrez, symbol.
        versions: Ensembl rows whose version equals, differs from, or is absent in the GTF's.
    """

    counts: pl.DataFrame
    genes: pl.DataFrame
    left_out: pl.DataFrame
    kinds: Mapping[str, int]
    versions: Mapping[str, int]


def is_sheet(sheet: pl.DataFrame) -> bool:
    """Whether a gather sheet came from a count matrix: it carries ``SOURCE_COLUMN``."""
    return SOURCE_COLUMN in sheet.columns


def refuse_sheet(sheet: pl.DataFrame, path: Path, command: str) -> list[str]:
    """Why a module that reads FASTQ files cannot read this sheet, or nothing."""
    if not is_sheet(sheet):
        return []
    return [f"{path.name} is a count matrix's sheet: ctx {command} reads FASTQ files, and a "
            "count matrix holds none; ctx qc, de, function and report read its counts"]  # fmt: skip


def _text(path: Path) -> str:
    """The file as text, gunzipped when its first bytes say gzip, whatever its name."""
    data = path.read_bytes()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return data.decode("utf-8-sig")


def _delimiter(path: Path) -> str:
    """Comma for a ``.csv`` or ``.csv.gz``, tab otherwise."""
    name = path.name.lower().removesuffix(".gz")
    return "," if name.endswith(".csv") else "\t"


def _body(text: str) -> tuple[str, str]:
    """The header line and the data lines, after leading ``#`` comment lines (featureCounts)."""
    lines = text.splitlines(keepends=True)
    start = next((i for i, line in enumerate(lines) if not line.startswith("#")), len(lines))
    header = lines[start] if start < len(lines) else ""
    return header, "".join(lines[start + 1 :])


def _columns(header: list[str]) -> tuple[list[str], tuple[str, ...], str]:
    """The header's names, the featureCounts columns to drop, and the format read."""
    names = [n.strip() for n in header]
    if tuple(names[1:6]) == FEATURECOUNTS:
        return names, FEATURECOUNTS, "featureCounts"
    return names, (), "table"


def _numeric(frame: pl.DataFrame) -> tuple[np.ndarray, list[Problem]]:
    """Every sample column as float64, and a problem per column holding a non-number."""
    values, problems = [], []
    for column in frame.columns:
        text = frame[column].str.strip_chars()
        number = text.cast(pl.Float64, strict=False)
        bad = text.is_not_null() & (text != "") & number.is_null()
        if bad.any():
            first = int(bad.arg_true()[0])
            problems.append(("values", column, f"{int(bad.sum()):,} cells are not numbers, "
                             f"the first {text[first]!r}"))  # fmt: skip
        values.append(number.fill_null(np.nan).to_numpy())
    matrix = np.column_stack(values) if values else np.empty((frame.height, 0))
    return matrix, problems


def read(path: Path) -> Raw:
    """Read a genes-by-samples matrix: gene IDs first, then one column per sample.

    TSV, CSV (by the ``.csv`` ending) and gzip (by the file's first bytes) are read, and so
    is featureCounts' own output, whose comment line and five annotation columns are dropped.
    HTSeq's and STAR's read-summary rows are removed and kept apart: they are not genes.

    Raises:
        CountsError: If the file holds no header or no sample column.
    """
    delimiter = _delimiter(path)
    line, body = _body(_text(path))
    names, dropped, kind = _columns(next(csv.reader([line], delimiter=delimiter), []))
    if len(names) < 2:  # noqa: PLR2004 - the gene IDs, then at least one sample
        read_as = "commas" if delimiter == "," else "tabs, as its name does not end in .csv"
        raise CountsError([("header", path.name, f"no sample column after the gene IDs, "
                            f"read with {read_as}")])  # fmt: skip
    positional = [f"c{i}" for i in range(len(names))]
    frame = pl.read_csv(io.StringIO(body), separator=delimiter, has_header=False,
                        new_columns=positional, infer_schema=False)  # fmt: skip
    ids = frame["c0"].fill_null("").str.strip_chars()
    keep = [c for c, n in zip(positional[1:], names[1:], strict=True) if n not in dropped]
    values, found = _numeric(frame.select(keep))
    problems = tuple((c, names[int(k[1:])], d) for c, k, d in found)
    summary = ids.str.starts_with(HTSEQ_PREFIX) | ids.is_in(list(STAR_ROWS))
    totals = np.nansum(values[summary.to_numpy()], axis=1)
    removed = tuple(zip(ids.filter(summary).to_list(), totals.tolist(), strict=True))
    samples = tuple(names[int(c[1:])] for c in keep)
    return Raw(tuple(ids.filter(~summary).to_list()), samples, values[~summary.to_numpy()],
               removed, dropped, kind, problems)  # fmt: skip


def _scaled(sums: np.ndarray) -> bool:
    """Whether every sample sums to a million within ``SCALED_TOLERANCE``: TPM or CPM."""
    return bool(sums.size) and bool(
        (np.abs(sums - PER_MILLION) <= SCALED_TOLERANCE * PER_MILLION).all()
    )


def value_problems(raw: Raw) -> list[Problem]:
    """Why the values are not raw counts, each with how it was told, or nothing.

    Raw counts are whole numbers of reads, none below 0. A fraction, a negative, a sum of a
    million in every sample or a ceiling of log scale each name the scale they point to.
    """
    v = raw.values
    finite = np.isfinite(v)
    problems: list[Problem] = []
    if (missing := int((~finite).sum())) > 0:
        said = f"{missing:,} cells are empty or not finite; a raw count is never missing"
        problems.append(("values", "missing", said))
    seen = v[finite]
    if (negative := int((seen < 0).sum())) > 0:
        said = f"{negative:,} values below 0, the lowest {seen.min():g}: counts never are"
        problems.append(("values", "negative", f"{said}; log2 or centred values are"))
    fractional = seen != np.rint(seen)
    sums = np.nansum(v, axis=0)
    if fractional.any():
        scale = "FPKM, RPKM or estimated counts"
        if _scaled(sums):
            scale = "TPM or CPM"
        elif seen.max() <= LOG_CEILING:
            scale = "log2 values"
        said = f"{int(fractional.sum()):,} of {seen.size:,} values are not whole numbers"
        shown = f"e.g. {seen[fractional][0]:g}; the largest is {seen.max():g}"
        problems.append(("values", "fractional", f"{said}, {shown}: {scale}, not raw counts"))
    elif _scaled(sums):
        said = f"every sample sums to {PER_MILLION:,} within {SCALED_TOLERANCE:.0%}"
        problems.append(("values", "scaled", f"{said}: rounded TPM or CPM, not raw counts"))
    return problems


def duplicate_problems(raw: Raw) -> list[Problem]:
    """Gene IDs and sample names that appear more than once, an Ensembl ID by its stable part."""
    problems: list[Problem] = []
    stable = [m[1] if (m := ENSEMBL.fullmatch(i)) else i for i in raw.ids]
    for kind, names in (("gene", stable), ("sample", list(raw.samples))):
        seen: dict[str, int] = defaultdict(int)
        for name in names:
            seen[name] += 1
        problems += [(f"duplicate {kind}", name, f"appears {n} times; a {kind} is one "
                      f"{'row' if kind == 'gene' else 'column'}")
                     for name, n in seen.items() if n > 1]  # fmt: skip
    if blank := sum(1 for i in raw.ids if not i):
        problems.append(("duplicate gene", "", f"{blank} rows have no gene ID"))
    return problems


def kind_of(source_id: str) -> str:
    """ensembl, entrez or symbol: what a matrix row's ID is, by its form alone."""
    if ENSEMBL.fullmatch(source_id):
        return "ensembl"
    return "entrez" if ENTREZ.fullmatch(source_id) else "symbol"


def entrez_genes(path: Path, gtf: Path) -> dict[str, set[str]]:
    """Each Entrez gene ID's GENCODE genes: its transcripts in the metadata, through the GTF.

    GENCODE's ``metadata.EntrezGene`` pairs transcripts with Entrez IDs; the GTF gives each
    transcript's gene. A transcript the GTF does not hold (another assembly's) is skipped.
    """
    gene_of = {}
    with gzip.open(gtf, "rt") if gtf.suffix == ".gz" else gtf.open() as handle:
        for line in handle:
            if "\ttranscript\t" in line and (t := TRANSCRIPT.search(line)):
                gene_of[t[1]] = GENE.search(line)[1]  # type: ignore[index]
    found: dict[str, set[str]] = defaultdict(set)
    with gzip.open(path, "rt") as handle:
        for line in handle:
            transcript, _, entrez = line.rstrip("\n").partition("\t")
            if transcript in gene_of:
                found[entrez.strip()].add(gene_of[transcript])
    return dict(found)


def _targets(annotation: pl.DataFrame, entrez: Mapping[str, set[str]] | None
             ) -> dict[str, dict[str, set[str]]]:  # fmt: skip
    """Per kind, each key's set of annotation genes: stable Ensembl ID, symbol, Entrez ID."""
    by: dict[str, dict[str, set[str]]] = {"ensembl": defaultdict(set), "symbol": defaultdict(set)}
    for gene_id, name in annotation.select("gene_id", "gene_name").rows():
        stable = ENSEMBL.fullmatch(gene_id)
        by["ensembl"][stable[1] if stable else gene_id].add(gene_id)
        if name:
            by["symbol"][name].add(gene_id)
    by["entrez"] = dict(entrez or {})
    return by


# (kind, gene_id or "", reason, detail) of one matrix row; the reason is one of REASONS.
Placed = tuple[str, str, str, str]
REASONS = ("not_in_annotation", "ambiguous", "shared_gene", "no_entrez_table")


def _place(ids: Iterable[str], targets: Mapping[str, Mapping[str, set[str]]],
           has_entrez: bool) -> list[Placed]:  # fmt: skip
    """Each row's one gene, or why it has none, with the genes or table it concerns."""
    placed: list[Placed] = []
    for source_id in ids:
        kind = kind_of(source_id)
        key = m[1] if (m := ENSEMBL.fullmatch(source_id)) else source_id
        genes = sorted(targets[kind].get(key, set()))
        if kind == "entrez" and not has_entrez:
            placed.append((kind, "", "no_entrez_table", f"no {ENTREZ_RESOURCE} resource"))
        elif len(genes) == 1:
            placed.append((kind, genes[0], "", ""))
        elif genes:
            placed.append((kind, "", "ambiguous", f"{len(genes)} genes: {', '.join(genes[:3])}"))
        else:
            placed.append((kind, "", "not_in_annotation", ""))
    return placed


def _versions(ids: Iterable[str], placed: Iterable[Placed]) -> dict[str, int]:
    """Ensembl rows mapped whose version equals the GTF's, differs, or is not given."""
    tally = {"equal": 0, "differ": 0, "unversioned": 0}
    for source_id, (kind, gene_id, _, _) in zip(ids, placed, strict=True):
        if kind != "ensembl" or not gene_id:
            continue
        given = ENSEMBL.fullmatch(source_id)
        ours = ENSEMBL.fullmatch(gene_id)
        if given is None or given[2] is None:
            tally["unversioned"] += 1
        else:
            tally["equal" if ours is not None and ours[2] == given[2] else "differ"] += 1
    return tally


def to_integers(values: np.ndarray) -> np.ndarray:
    """Whole counts as int32 where every one fits, int64 otherwise.

    A value that is not finite becomes 0 here: such a matrix has already failed
    ``value_problems``, and nothing built from it is written.
    """
    whole = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
    kind = np.int32 if whole.max(initial=0) <= INT32_MAX else np.int64
    integers: np.ndarray = whole.astype(kind)
    return integers


def _shared(placed: list[Placed], ids: tuple[str, ...]) -> None:
    """Rows placed on one gene leave together: summing them would be a guess."""
    hits: dict[str, list[int]] = defaultdict(list)
    for i, (_, gene_id, _, _) in enumerate(placed):
        if gene_id:
            hits[gene_id].append(i)
    for gene_id, rows in hits.items():
        if len(rows) > 1:
            detail = f"{gene_id}, with {', '.join(ids[r] for r in rows)}"
            for r in rows:
                placed[r] = (placed[r][0], "", "shared_gene", detail)


def place(raw: Raw, annotation: pl.DataFrame, entrez: Mapping[str, set[str]] | None) -> Mapped:
    """Each row on one annotation gene, or left out with why; two rows on one gene both leave.

    Order: each ID by its own kind (Ensembl by its stable part, Entrez through GENCODE's
    table, a symbol by the annotation's gene name); a key matching several genes is ambiguous;
    then rows sharing a gene are left out together.
    """
    names = dict(annotation.select("gene_id", "gene_name").rows())
    placed = _place(raw.ids, _targets(annotation, entrez), entrez is not None)
    _shared(placed, raw.ids)
    kept = [i for i, (_, g, _, _) in enumerate(placed) if g]
    totals = np.nansum(raw.values, axis=1)
    left = [(raw.ids[i], k, why, detail, float(totals[i]))
            for i, (k, g, why, detail) in enumerate(placed) if not g]  # fmt: skip
    gene_ids = [placed[i][1] for i in kept]
    matrix = to_integers(raw.values[kept]) if kept else np.empty((0, len(raw.samples)), np.int32)
    data = {s: pl.Series(s, matrix[:, j]) for j, s in enumerate(raw.samples)}
    counts = pl.DataFrame({"gene_id": pl.Series("gene_id", gene_ids, pl.String)} | data)
    genes = pl.DataFrame({"gene_id": gene_ids, "gene_name": [names.get(g) for g in gene_ids],
                          "source_id": [raw.ids[i] for i in kept]},
                         schema=dict.fromkeys(("gene_id", "gene_name", "source_id"),
                                              pl.String))  # fmt: skip
    schema = {**dict.fromkeys(("source_id", "kind", "reason", "detail"), pl.String),
              "total": pl.Float64}  # fmt: skip
    left_out = pl.DataFrame(left, schema=schema, orient="row")
    kinds: dict[str, int] = defaultdict(int)
    for kind, *_ in placed:
        kinds[kind] += 1
    return Mapped(counts, genes, left_out, dict(kinds), _versions(raw.ids, placed))
