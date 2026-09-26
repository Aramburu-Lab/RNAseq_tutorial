#!/usr/bin/env python3
# /// codebase-meta
# exempt: instructor-side builder of tests/fixtures/chr22; needs the network
# ///
"""Build tests/fixtures/chr22: a GRCh38 chr22 slice, its annotation, and reads with known truth.

Instructor-side and online. It downloads the pinned GENCODE annotation and the Ensembl
chromosome 22 sequence of the same release, checks both against the publishers' checksums
and a pinned SHA-256, and writes an uncompressed fixture that the tests read offline.
Standard library only, so it runs before anything else is installed.

Usage:
    pixi run python scripts/python/make_chr22_fixture.py --out tests/fixtures/chr22
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import random
import re
import shutil
import sys
import tempfile
import urllib.request
from bisect import bisect_right
from collections.abc import Iterable, Iterator
from itertools import accumulate
from pathlib import Path
from typing import NamedTuple

GENCODE_RELEASE = 50
# GENCODE 50 is Ensembl 116 on GRCh38.p14 (gencodegenes.org/human/releases.html, and the
# GTF's own header line), so Ensembl's chromosome 22 carries the annotation's coordinates.
ENSEMBL_RELEASE = 116
GENCODE_DIR = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50"
ENSEMBL_DIR = "https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/dna"


class Source(NamedTuple):
    """A pinned download and the checksums that prove it is the published file."""

    url: str
    algorithm: str  # of the publisher's checksum: "md5" (GENCODE MD5SUMS), "bsd_sum" (CHECKSUMS)
    publisher: str
    sha256: str  # first-download value, 2026-09-26; the publishers give none


GTF = Source(
    f"{GENCODE_DIR}/gencode.v50.basic.annotation.gtf.gz",
    "md5",
    "11e77cf1d4b78b43f57d608e3f7f2202",
    "5b9ad9aa71f6ab5d55123c51f74723644a895bc0cbaa9cee8fbc11bb9dec0bf7",
)
FASTA = Source(
    f"{ENSEMBL_DIR}/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz",
    "bsd_sum",
    "31309 11123",
    "974f97ac8ef7ffae971b63b47608feda327403be40c27e391ee4a1a78b800df5",
)
CHROM = "chr22"
# 1-based closed. Both edges fall between genes, so no transcript is cut, and genes inside share
# exons on opposite strands (the README lists them): strandedness decides those reads.
WINDOW_START, WINDOW_END = 30_330_001, 30_431_000
OFFSET = WINDOW_START - 1
CONTIG = "chr22_slice"
READ_LENGTH = 100
# Uniform fragment lengths; the minimum exceeds the read length, so no read runs into adapter.
FRAGMENT_MIN, FRAGMENT_MAX = 150, 350
N_PAIRS = 2500
SEED = 42
QUALITY = "I"
FASTA_WIDTH = 60
MANIFEST = "MANIFEST.sha256"
_ATTR = re.compile(r'(\S+) "([^"]*)"')
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


class Transcript(NamedTuple):
    """One mature transcript, 5' to 3'."""

    transcript_id: str
    gene_id: str
    gene_name: str
    sequence: str


def revcomp(sequence: str) -> str:
    """Reverse complement of an ACGT sequence."""
    return sequence.translate(_COMPLEMENT)[::-1]


def digest(path: Path, algorithm: str) -> str:
    """Hex digest of a file, read in chunks."""
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            hasher.update(chunk)
    return hasher.hexdigest()


def bsd_sum(path: Path) -> str:
    """BSD ``sum``: 16-bit rotating checksum and 1 KiB block count, as Ensembl's CHECKSUMS."""
    checksum = size = 0
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            size += len(chunk)
            for byte in chunk:
                checksum = (((checksum >> 1) | ((checksum & 1) << 15)) + byte) & 0xFFFF
    return f"{checksum:05d} {-(-size // 1024)}"


def fetch(source: Source, cache: Path) -> Path:
    """Download a pinned file once into the cache, then prove it is the published one.

    Args:
        source: The pinned URL and checksums.
        cache: Directory outside the repository for downloads.

    Returns:
        Path of the verified file.

    Raises:
        ValueError: If the file fails the publisher's checksum or the pinned SHA-256.
    """
    path = cache / source.url.rsplit("/", 1)[1]
    if not path.exists():
        cache.mkdir(parents=True, exist_ok=True)
        partial = path.with_name(path.name + ".part")
        urllib.request.urlretrieve(source.url, partial)
        partial.rename(path)
    found = digest(path, "md5") if source.algorithm == "md5" else bsd_sum(path)
    if found != source.publisher:
        raise ValueError(f"{path}: {source.algorithm} {found}, expected {source.publisher}")
    if (sha := digest(path, "sha256")) != source.sha256:
        raise ValueError(f"{path}: sha256 {sha}, expected {source.sha256}")
    return path


def read_chromosome(path: Path) -> str:
    """The one sequence of a gzipped single-chromosome FASTA, upper case as Ensembl ships it."""
    with gzip.open(path, "rt") as handle:
        handle.readline()
        return "".join(line.rstrip("\n") for line in handle)


def gtf_rows(path: Path, chrom: str) -> Iterator[list[str]]:
    """The tab-split rows of one chromosome in a gzipped GTF."""
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(chrom + "\t"):
                yield line.rstrip("\n").split("\t")


def window_features(rows: Iterable[list[str]], start: int, end: int) -> list[list[str]]:
    """Every row of every gene lying wholly inside [start, end], at genomic coordinates.

    Args:
        rows: GTF rows of one chromosome, each gene's row before its features, as GENCODE writes.
        start: First base of the window, 1-based.
        end: Last base of the window, inclusive.

    Returns:
        The kept rows in file order.

    Raises:
        ValueError: If a gene crosses a window edge, which would leave half a gene unannotated.
    """
    inside: dict[str, bool] = {}
    kept = []
    for row in rows:
        gene = dict(_ATTR.findall(row[8]))["gene_id"]
        if row[2] == "gene":
            low, high = int(row[3]), int(row[4])
            inside[gene] = start <= low and high <= end
            if not inside[gene] and low <= end and high >= start:
                raise ValueError(f"{gene} ({low}-{high}) crosses the window {start}-{end}")
        if inside[gene]:
            kept.append(row)
    return kept


def rebase(rows: Iterable[list[str]], offset: int, contig: str) -> list[list[str]]:
    """Move rows onto a contig that starts at genomic position offset + 1."""
    return [
        [contig, *r[1:3], str(int(r[3]) - offset), str(int(r[4]) - offset), *r[5:]] for r in rows
    ]


def splice(sequence: str, rows: Iterable[list[str]]) -> list[Transcript]:
    """Join each transcript's exons in transcript order, reverse-complementing minus-strand exons.

    Args:
        sequence: The contig the rows' 1-based closed coordinates refer to.
        rows: GTF rows; only exon rows are read.

    Returns:
        Transcripts in the order their first exon appears.

    Raises:
        ValueError: If an exon lies outside the sequence, where slicing would silently truncate.
    """
    exons: dict[tuple[str, str, str], list[tuple[int, int, str]]] = {}
    for row in rows:
        if row[2] != "exon":
            continue
        low, high = int(row[3]), int(row[4])
        if not 1 <= low <= high <= len(sequence):
            raise ValueError(f"exon {low}-{high} outside a sequence of {len(sequence)} bases")
        attrs = dict(_ATTR.findall(row[8]))
        key = (attrs["transcript_id"], attrs["gene_id"], attrs["gene_name"])
        exons.setdefault(key, []).append((low, high, row[6]))
    transcripts = []
    for key, parts in exons.items():
        minus = parts[0][2] == "-"
        ordered = sorted(parts, reverse=minus)
        pieces = [sequence[low - 1 : high] for low, high, _ in ordered]
        transcripts.append(Transcript(*key, "".join(revcomp(p) if minus else p for p in pieces)))
    return transcripts


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    """Write (header, sequence) pairs, 60 bases per line, with Unix line ends on every OS."""
    with path.open("w", newline="\n") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for i in range(0, len(sequence), FASTA_WIDTH):
                handle.write(sequence[i : i + FASTA_WIDTH] + "\n")


def read_transcripts(path: Path) -> list[Transcript]:
    """Parse a transcripts.fa whose headers read ``<id> gene_id=<id> gene_name=<name>``."""
    records: list[tuple[str, list[str]]] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            records.append((line[1:], []))
        else:
            records[-1][1].append(line)
    out = []
    for header, lines in records:
        name, *fields = header.split()
        attrs = dict(field.split("=", 1) for field in fields)
        out.append(Transcript(name, attrs["gene_id"], attrs["gene_name"], "".join(lines)))
    return out


def require_empty(directory: Path) -> None:
    """Refuse a directory that holds anything, so a rerun never mixes old files with new."""
    if directory.exists() and any(directory.iterdir()):
        raise FileExistsError(f"{directory} is not empty")


def reads(transcripts: Path, n_pairs: int, seed: int = SEED, *, out_dir: Path) -> dict[str, int]:
    """Simulate error-free, reverse-stranded (dUTP) read pairs whose names carry their origin.

    Only ``random()`` of a seeded ``random.Random`` is called, with plain arithmetic on top:
    Python guarantees that method's sequence for a given seed across versions, and not
    ``choice``, ``randrange`` or ``choices``. Names are ``<pair>:<transcript>:<start>:<length>``
    with a 0-based fragment start on the transcript, then ``/1`` or ``/2``.

    Args:
        transcripts: transcripts.fa as ``write_fasta`` writes it.
        n_pairs: Read pairs to simulate.
        seed: Seed for ``random.Random``.
        out_dir: Empty or missing directory for reads_R1.fq, reads_R2.fq and truth.tsv.

    Returns:
        Pairs per transcript, 0 for transcripts shorter than the longest fragment.
    """
    require_empty(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    records = read_transcripts(transcripts)
    rng = random.Random(seed)
    eligible = [r for r in records if len(r.sequence) >= FRAGMENT_MAX]
    # Fragment shares from 1 to 64 (a squared uniform) give low and high expressers. Products,
    # not pow(): IEEE multiplication is exactly rounded on every platform, libm pow() is not.
    shares = [(1 + 7 * u) * (1 + 7 * u) for u in (rng.random() for _ in eligible)]
    cumulative = list(accumulate(shares))
    counts = dict.fromkeys((r.transcript_id for r in records), 0)
    with (
        (out_dir / "reads_R1.fq").open("w", newline="\n") as read1,
        (out_dir / "reads_R2.fq").open("w", newline="\n") as read2,
    ):
        for pair in range(1, n_pairs + 1):
            origin = eligible[bisect_right(cumulative, rng.random() * cumulative[-1])]
            length = FRAGMENT_MIN + int(rng.random() * (FRAGMENT_MAX - FRAGMENT_MIN + 1))
            start = int(rng.random() * (len(origin.sequence) - length + 1))
            fragment = origin.sequence[start : start + length]
            name = f"{pair:06d}:{origin.transcript_id}:{start}:{length}"
            # dUTP: read 1 is antisense, from the fragment's 3' end; read 2 is sense, from its 5'.
            read1.write(
                f"@{name}/1\n{revcomp(fragment[-READ_LENGTH:])}\n+\n{QUALITY * READ_LENGTH}\n"
            )
            read2.write(f"@{name}/2\n{fragment[:READ_LENGTH]}\n+\n{QUALITY * READ_LENGTH}\n")
            counts[origin.transcript_id] += 1
    with (out_dir / "truth.tsv").open("w", newline="\n") as truth:
        truth.write("transcript_id\tgene_id\tn_pairs\n")
        truth.writelines(
            f"{r.transcript_id}\t{r.gene_id}\t{counts[r.transcript_id]}\n" for r in records
        )
    return counts


def shared_exons(rows: Iterable[list[str]]) -> list[str]:
    """Exonic bases each pair of genes shares: the reads a gene-level count has to arbitrate."""
    bases: dict[tuple[str, str], set[int]] = {}
    for row in rows:
        if row[2] == "exon":
            gene = (dict(_ATTR.findall(row[8]))["gene_name"], row[6])
            bases.setdefault(gene, set()).update(range(int(row[3]), int(row[4]) + 1))
    genes = list(bases)
    return [
        f"- {a[0]} ({a[1]}) and {b[0]} ({b[1]}): {shared} bp"
        for i, a in enumerate(genes)
        for b in genes[i + 1 :]
        if (shared := len(bases[a] & bases[b]))
    ]


def readme(
    annotation: list[list[str]],
    transcripts: list[Transcript],
    counts: dict[str, int],
    out_dir: Path,
) -> str:
    """The fixture's README. It carries no timestamp, so a rebuild reproduces it byte for byte.

    Args:
        annotation: The rows of slice.gtf.
        transcripts: What transcripts.fa holds.
        counts: Pairs per transcript, from ``reads``.
        out_dir: The fixture directory, holding every data file and nothing else yet.

    Returns:
        The README text.
    """
    per_gene: dict[str, list[int]] = {}
    for t in transcripts:
        per_gene.setdefault(t.gene_id, []).append(counts[t.transcript_id])
    gene_rows = []
    for g in (row for row in annotation if row[2] == "gene"):
        attrs = dict(_ATTR.findall(g[8]))
        n = per_gene.get(attrs["gene_id"], [])
        gene_rows.append(
            f"| {attrs['gene_name']} | {attrs['gene_id']} | {g[6]} | {attrs['gene_type']} "
            f"| {g[3]} | {g[4]} | {len(n)} | {sum(n)} |"
        )
    files = sorted(out_dir.iterdir())
    return README.format(
        **globals(),
        WIDTH=WINDOW_END - WINDOW_START + 1,
        GENE_ROWS="\n".join(gene_rows),
        SHARED="\n".join(shared_exons(annotation)),
        N_TRANSCRIPTS=len(transcripts),
        UNSAMPLED=sum(1 for t in transcripts if len(t.sequence) < FRAGMENT_MAX),
        FILES="\n".join(
            f"| {p.name} | {p.stat().st_size} | {digest(p, 'sha256')} |" for p in files
        ),
    )


README = """# chr22 fixture: a GRCh38 slice and simulated reads of known origin

Written by `scripts/python/make_chr22_fixture.py`; never edit it by hand. To rebuild, run
`pixi run python scripts/python/make_chr22_fixture.py --out <empty directory>` and compare
`MANIFEST.sha256`, or check this copy with `shasum -a 256 -c MANIFEST.sha256`. The sequence is
the reference genome and the reads are simulated: no sample or study is involved.

## Source

| | |
|---|---|
| Annotation | GENCODE {GENCODE_RELEASE} basic, GRCh38.p14 |
| URL | {GTF.url} |
| Publisher checksum | MD5 `{GTF.publisher}`, from `MD5SUMS` in the same directory |
| SHA-256 | `{GTF.sha256}` (first download, 2026-09-26; GENCODE publishes none) |
| Sequence | Ensembl {ENSEMBL_RELEASE}, chromosome 22, unmasked (`dna`) |
| URL | {FASTA.url} |
| Publisher checksum | BSD `sum` `{FASTA.publisher}`, from `CHECKSUMS` in the same directory |
| SHA-256 | `{FASTA.sha256}` (first download, 2026-09-26; Ensembl publishes none) |

GENCODE {GENCODE_RELEASE} is Ensembl {ENSEMBL_RELEASE}
(https://www.gencodegenes.org/human/releases.html; the GTF header reads
`version {GENCODE_RELEASE} (Ensembl {ENSEMBL_RELEASE})`). Ensembl names the chromosome `22` and
GENCODE `chr22`; the coordinates are the same.

## Window

`{CHROM}:{WINDOW_START}-{WINDOW_END}`, 1-based and closed, {WIDTH} bp, written as the contig
`{CONTIG}`. Offset {OFFSET}: {CHROM} position = {CONTIG} position + {OFFSET}. Both edges fall
between genes, and every gene lying wholly inside is kept, whatever its biotype. It was picked
by a scan of chr22 for the narrowest run of whole genes, 2 kb clear of other genes at both
ends, holding at least two multi-exon protein-coding genes on each strand in at most 45
transcripts; the edges were then rounded.

| gene | gene_id | strand | type | slice start | slice end | transcripts | pairs |
|---|---|---|---|---|---|---|---|
{GENE_ROWS}

`transcripts.fa` holds the {N_TRANSCRIPTS} transcripts spliced from `slice.fa` at the
coordinates in `slice.gtf`, minus-strand ones reverse-complemented.

Exonic bases shared by two genes, counted from `slice.gtf`. On opposite strands a stranded
count separates the genes and an unstranded one cannot; on the same strand neither can.

{SHARED}

## Reads

- Seed {SEED}, {N_PAIRS} pairs, 2 x {READ_LENGTH} nt, fragment length uniform on
  {FRAGMENT_MIN}-{FRAGMENT_MAX} nt.
- Reverse-stranded (dUTP): read 1 is antisense, from the fragment's 3' end; read 2 is sense.
  Salmon library type `ISR`; featureCounts `-s 2`.
- No sequencing errors and no adapter (every fragment is longer than a read); quality `I`
  (Phred 40) on every base.
- Read name `<pair>:<transcript_id>:<start>:<length>/1` (or `/2`); start is 0-based on the
  transcript, 5' to 3'.
- Each transcript of at least {FRAGMENT_MAX} nt gets a fragment share from 1 to 64 (a squared
  uniform draw); the {UNSAMPLED} shorter ones get no reads and stay in `truth.tsv` at 0.
- `truth.tsv`: `transcript_id`, `gene_id`, `n_pairs`.

## What it can test

- Parsing and pairing of FASTQ, FASTA and GTF, and their agreement with each other.
- Spliced alignment on both strands: every read's transcript and offset are in its name.
- Strandedness: the wrong strand setting loses reads, and the opposite-strand pairs above
  show what an unstranded count cannot separate.
- Transcript and gene quantification against exact counts, including zeros.
- Byte-identical regeneration from the seed.

## What it cannot test

- Sequencing errors, adapters, low-quality tails or trimming: none are planted.
- Intronic, intergenic, pre-mRNA, rRNA or genomic DNA reads: every read is from a mature
  transcript.
- Multi-mapping to paralogues or pseudogenes elsewhere: the reference is {WIDTH} bp.
- GC, positional or 3' bias, PCR duplicates, a realistic fragment-length distribution.
- Anything between samples, such as differential expression: there is one sample.
- Novel transcripts, or transcripts shorter than {FRAGMENT_MAX} nt, which get no reads.

## Files

| file | bytes | SHA-256 |
|---|---|---|
{FILES}
"""


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Read the command line.

    Args:
        argv: Argument list without the program name; ``sys.argv[1:]`` when None.

    Returns:
        The parsed arguments.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", required=True, type=Path, help="Empty fixture directory.")
    parser.add_argument(
        "--cache",
        type=Path,
        default=Path(tempfile.gettempdir()) / "chr22_fixture_cache",
        help="Download cache, outside the repository.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Download, verify, slice, annotate, simulate, and write the fixture with its manifest.

    Args:
        argv: Argument list without the program name; ``sys.argv[1:]`` when None.

    Returns:
        0 on success.
    """
    args = parse_args(argv)
    require_empty(args.out)
    genomic = window_features(gtf_rows(fetch(GTF, args.cache), CHROM), WINDOW_START, WINDOW_END)
    annotation = rebase(genomic, OFFSET, CONTIG)
    contig = read_chromosome(fetch(FASTA, args.cache))[WINDOW_START - 1 : WINDOW_END]
    if set(contig) - set("ACGT"):
        raise ValueError(f"the window holds {sorted(set(contig) - set('ACGT'))}, not only ACGT")
    # Extracted from the slice at the rebased coordinates, so all three files share any rebase
    # error and agree with each other. Offline, the codon anchor in test_gtf_fasta_consistent
    # catches a wrong offset; test_slice_rebuilds checks against the genome itself.
    transcripts = splice(contig, annotation)
    with tempfile.TemporaryDirectory() as staging_dir:
        staging = Path(staging_dir)
        write_fasta(
            staging / "slice.fa", [(f"{CONTIG} {CHROM}:{WINDOW_START}-{WINDOW_END}", contig)]
        )
        write_fasta(
            staging / "transcripts.fa",
            [
                (f"{t.transcript_id} gene_id={t.gene_id} gene_name={t.gene_name}", t.sequence)
                for t in transcripts
            ],
        )
        with (staging / "slice.gtf").open("w", newline="\n") as gtf:
            gtf.write(
                f"##description: GENCODE {GENCODE_RELEASE} basic (Ensembl {ENSEMBL_RELEASE}), "
                f"genes wholly inside {CHROM}:{WINDOW_START}-{WINDOW_END}\n"
                f"##offset: {OFFSET} ({CHROM} position = {CONTIG} position + offset)\n"
            )
            gtf.writelines("\t".join(row) + "\n" for row in annotation)
        counts = reads(staging / "transcripts.fa", N_PAIRS, SEED, out_dir=args.out)
        for staged in staging.iterdir():
            shutil.move(staged, args.out / staged.name)
    text = readme(annotation, transcripts, counts, args.out)
    (args.out / "README.md").write_text(text, newline="\n")
    manifest = "".join(f"{digest(p, 'sha256')}  {p.name}\n" for p in sorted(args.out.iterdir()))
    (args.out / MANIFEST).write_text(manifest, newline="\n")
    print(f"{len(transcripts)} transcripts, {sum(counts.values())} pairs", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
