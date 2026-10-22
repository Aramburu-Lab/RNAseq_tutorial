"""One BAM's featureCounts matrix, compacted where it was counted (ED.6 pilot, 2026-10-22).

On GENCODE 50 one sample's matrix at one setting is 133 MB, nearly all of it the annotation
columns (Chr, Start, End, Strand, Length), the same in every file counted against one annotation.
``compact`` keeps the command line, Geneid and the count column, gzipped, and the annotation
columns once beside them; it reads the compact file back against the original, row by row, and
only then deletes the original. ``ctx.process.cluster`` reads the compact form and rebuilds the
full matrix from the one annotation copy.

Standard library only, and no syntax newer than Python 3.6: blocks B1 and E run this file on
Kebnekaise with the system ``python3``, outside pixi.
"""

from __future__ import annotations

import gzip
import hashlib
import os
import sys

SUFFIX = ".featureCounts.tsv.gz"
ANNOTATION = "featurecounts_annotation.tsv.gz"
MATRIX = ".featureCounts.txt"
META = 6


def compact_name(matrix: str) -> str:
    """The compact file of a featureCounts matrix: ``<setting>.featureCounts.tsv.gz``."""
    return matrix[: -len(MATRIX)] + SUFFIX


def _digest(text: str) -> str:
    """SHA-256 of a text, as UTF-8."""
    return hashlib.sha256(text.encode()).hexdigest()


def _annotation(meta: str, path: str) -> None:
    """Write the annotation columns once; every later file's must be the same.

    Raises:
        ValueError: If a file's annotation columns differ from the copy already kept.
    """
    if os.path.exists(path):
        with gzip.open(path, "rt") as handle:
            kept = handle.read()
        if _digest(kept) != _digest(meta):
            raise ValueError(f"the annotation columns differ from {path}: another annotation")
        return
    # Written beside and renamed, so array tasks racing for the first copy never half-write it.
    part = f"{path}.{os.getpid()}.part"
    with gzip.open(part, "wt") as handle:
        handle.write(meta)
    os.replace(part, path)


def compact(matrix: str, annotation: str) -> str:
    """Compact one featureCounts matrix of one BAM, check it, and delete the original.

    Args:
        matrix: ``<setting>.featureCounts.txt``: featureCounts' command line, its header
            (Geneid, Chr, Start, End, Strand, Length, the BAM), one row per gene.
        annotation: The one copy of the annotation columns for every file of this annotation.

    Returns:
        The compact file: the command line, ``Geneid`` and the BAM's header, one row per gene.

    Raises:
        ValueError: If the matrix is not one BAM's, its annotation columns differ from the copy,
            or the compact file read back differs from the original in a row or in the sum.
    """
    with open(matrix) as handle:
        command, header, *rows = handle.read().splitlines()
    fields = [row.split("\t") for row in [header, *rows]]
    if any(len(f) != META + 1 for f in fields):
        raise ValueError(f"{matrix}: not a matrix of one BAM ({len(fields[0])} header columns)")
    _annotation("".join("\t".join(f[:META]) + "\n" for f in fields), annotation)
    out = compact_name(matrix)
    with gzip.open(out, "wt") as handle:
        handle.write(command + "\n")
        handle.writelines(f"{f[0]}\t{f[META]}\n" for f in fields)
    _, back = read(out)
    wanted = [(f[0], f[META]) for f in fields]
    if back != wanted or sum(map(float, (c for _, c in back[1:]))) != sum(
        float(c) for _, c in wanted[1:]
    ):
        raise ValueError(f"{out} read back differs from {matrix}")
    os.remove(matrix)
    return out


def read(path: str) -> tuple[str, list[tuple[str, str]]]:
    """A compact file's command line and its rows, the header first, each (Geneid, count)."""
    with gzip.open(path, "rt") as handle:
        command, *rows = handle.read().splitlines()
    return command, [(row.split("\t")[0], row.split("\t")[1]) for row in rows]


def annotation_rows(path: str) -> list[list[str]]:
    """The kept annotation columns, the header first."""
    with gzip.open(path, "rt") as handle:
        return [line.split("\t") for line in handle.read().splitlines()]


def main(argv: list[str]) -> int:
    """``python3 fc_compact.py <annotation copy> <matrix> ...``: compact each matrix in turn."""
    annotation, *matrices = argv
    for matrix in matrices:
        out = compact(matrix, annotation)
        print(f"{out}: {os.path.getsize(out)} bytes; {matrix} checked and removed")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
