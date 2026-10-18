"""RSEM's expected counts at gather's door: per-sample ``.genes.results`` files and the rounding.

A dataset file declares ``[counts] kind = "rsem_expected"`` for a matrix of RSEM's
``expected_count``, or for a folder of RSEM's per-sample gene files, combined here by gene ID,
never by row order. Only the ``expected_count`` column is read: TPM and FPKM never become
counts. The values are checked as they are, then rounded to the nearest integer, as is common
before DESeq2, and the rounding is counted for the page and every README.
"""

import gzip
import io
import re
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import polars as pl

from ctx.dataset import CountsSpec
from ctx.gather import counts
from ctx.gather.curation import Problem

# RSEM's genes.results starts with these two columns (rsem-calculate-expression's output).
RSEM_FIRST = ("gene_id", "transcript_id(s)")
EXPECTED = "expected_count"
SUFFIXES = (".genes.results", ".genes.results.gz")
FORMAT = "RSEM genes.results"


@dataclass(frozen=True)
class Rounding:
    """What rounding did to the matrix as read, every row included.

    Attributes:
        values: The finite values rounded.
        fractional: How many of them were not whole numbers.
        largest: The largest absolute change, at most 0.5.
    """

    values: int
    fractional: int
    largest: float


def _first_line(path: Path) -> str:
    """The file's first line, gunzipped when its first bytes say gzip."""
    with path.open("rb") as handle:
        gzipped = handle.read(2) == b"\x1f\x8b"
    with gzip.open(path, "rt") if gzipped else path.open(encoding="utf-8-sig") as handle:
        return handle.readline().rstrip("\r\n")


def is_rsem(path: Path) -> bool:
    """Whether a file is RSEM's genes.results: its header starts ``gene_id, transcript_id(s)``."""
    return path.is_file() and tuple(_first_line(path).split("\t")[:2]) == RSEM_FIRST


def sample_name(path: Path, pattern: str) -> str | None:
    """The sample a file holds, or None when a pattern is given and does not match.

    The pattern's first group in the file's name, else the name less its RSEM ending.
    """
    if not pattern:
        name = path.name.removesuffix(".gz")
        return name.removesuffix(SUFFIXES[0])
    found = re.search(pattern, path.name)
    return found.group(1) if found and found.groups() else None


def _one(path: Path) -> tuple[pl.DataFrame | None, list[Problem]]:
    """One file's ``gene_id`` and ``expected_count``, read by name, or why it cannot be."""
    frame = pl.read_csv(io.StringIO(counts.file_text(path)), separator="\t", infer_schema=False)
    if missing := [c for c in (RSEM_FIRST[0], EXPECTED) if c not in frame.columns]:
        return None, [("rsem", path.name, f"no {missing} column; not RSEM's genes.results")]
    frame = frame.select(gene_id=pl.col(RSEM_FIRST[0]).str.strip_chars(),
                         value=pl.col(EXPECTED).str.strip_chars())  # fmt: skip
    problems: list[Problem] = []
    if (repeated := frame.height - frame["gene_id"].n_unique()) > 0:
        problems.append(("duplicate gene", path.name, f"{repeated:,} gene IDs repeat in it"))
    return frame, problems


def _combined(files: list[tuple[str, Path]]) -> counts.Raw:
    """The files' expected counts as one matrix, joined on gene ID; the first file's order.

    A gene in one file and not another is a problem, never a zero: RSEM writes every gene of
    its reference, so differing gene sets mean different references.
    """
    problems: list[Problem] = []
    read = []
    for sample, path in files:
        frame, found = _one(path)
        problems += found
        if frame is not None:
            read.append((sample, path, frame))
    if not read:
        raise counts.CountsError(problems)
    ids = read[0][2]["gene_id"].to_list()
    first = set(ids)
    columns = []
    for sample, path, frame in read:
        given = set(frame["gene_id"])
        if given != first:
            problems.append(("rsem", path.name, f"{len(first - given):,} genes of "
                             f"{read[0][1].name} are absent, {len(given - first):,} added; "
                             "the files come from different references"))  # fmt: skip
        by_id = dict(frame.unique("gene_id", keep="first").rows())
        columns.append(pl.Series(sample, [by_id.get(i) for i in ids], pl.String))
    values, found = counts.numeric_columns(pl.DataFrame(columns))
    dropped = ("transcript_id(s)", "length", "effective_length", "TPM", "FPKM")
    return counts.Raw(tuple(ids), tuple(s for s, _, _ in read), values, (), dropped, FORMAT,
                      tuple(problems + found))  # fmt: skip


def read(matrix: Path, spec: CountsSpec) -> counts.Raw:
    """A matrix, one RSEM genes.results file, or a folder of them, as ``counts.Raw``.

    Raises:
        counts.CountsError: If a folder holds no RSEM file, a file name misses the sample
            pattern, two files name one sample, or no file can be read.
    """
    if not matrix.is_dir() and not is_rsem(matrix):
        return counts.read(matrix)
    paths = sorted(p for p in matrix.iterdir() if p.name.endswith(SUFFIXES)) \
        if matrix.is_dir() else [matrix]  # fmt: skip
    if not paths:
        raise counts.CountsError([("rsem", matrix.name, f"no file ending {' or '.join(SUFFIXES)}")])
    named = [(sample_name(p, spec.sample_pattern), p) for p in paths]
    said = f"sample_pattern {spec.sample_pattern!r} finds no group in the name"
    problems: list[Problem] = [("rsem", p.name, said) for s, p in named if s is None]
    samples = [s for s, _ in named if s is not None]
    problems += [("duplicate sample", s, f"{samples.count(s)} files name it")
                 for s in sorted(set(samples)) if samples.count(s) > 1]  # fmt: skip
    if problems:
        raise counts.CountsError(problems)
    return _combined([(s, p) for s, p in named if s is not None])


def kind_problems(raw: counts.Raw, kind: str) -> list[Problem]:
    """An RSEM file read without ``kind = "rsem_expected"`` declared: the dataset file decides."""
    if raw.format == FORMAT and kind != counts.RSEM_EXPECTED:
        return [("values", "kind", f'RSEM\'s expected_count needs [counts] kind = '
                 f'"{counts.RSEM_EXPECTED}" in the dataset file')]  # fmt: skip
    return []


def rounded(raw: counts.Raw) -> tuple[counts.Raw, Rounding]:
    """Every value at its nearest integer, halves to even as R's ``round`` takes them.

    Called after ``counts.value_problems`` has passed the values as they were read.
    """
    whole = np.rint(raw.values)
    finite = np.isfinite(raw.values)
    change = np.abs(whole - raw.values)[finite]
    found = Rounding(int(finite.sum()), int((change > 0).sum()), float(change.max(initial=0.0)))
    return replace(raw, values=whole), found
