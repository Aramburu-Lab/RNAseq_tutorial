"""The reference M3 aligns to: the genome, the annotation's lines on its contigs, the indices.

The genome is a FASTA given by path (plain or gzipped) and recorded by its SHA-256; its contig
names choose the annotation's lines, matched exactly. The annotation is the dataset's declared
GTF. Both are written into one content-keyed reference folder (``ctx.process.cache``): the
genome uncompressed, because STAR refuses gzipped FASTA, and the GTF lines on its contigs. With
``half``, each contig keeps its first half and the annotation only the genes that end inside it,
the "incomplete reference" the page teaches. The start codons check that genome and annotation
share one assembly: on a shifted or foreign genome they stop reading ATG.
"""

import gzip
import hashlib
import json
import math
import re
from collections import Counter
from collections.abc import Collection, Iterator, Mapping
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import IO, Any

from ctx import tools
from ctx.process import ProcessError, cache

GENOME = "genome.fa"
ANNOTATION = "annotation.gtf"
FACTS = "reference.json"
HISAT2_PREFIX = "genome"
FASTA_WIDTH = 60
# STAR's default sjdbOverhang, "ideally mate_length - 1" (manual 2.1): Step-1's longest reads
# are 101 nt, so the ideal is the default.
OVERHANG = 100
# STAR manual 2.2.5: genomeSAindexNbases = min(14, log2(GenomeLength)/2 - 1).
SA_INDEX_CAP = 14
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
GZIP_MAGIC = b"\x1f\x8b"
GENE_ID = re.compile(r'gene_id "([^"]+)"')
# GENCODE's "...version 50 (Ensembl 116)"; the test fixture's "GENCODE 50 basic (Ensembl 116)".
RELEASE = re.compile(r"(?:version|GENCODE) (\d+)\b[^(\n]*\(Ensembl (\d+)\)")
SHOWN = 10


@dataclass(frozen=True)
class Genome:
    """A FASTA file's facts: its SHA-256 as given, its size, and each contig's name and length."""

    sha256: str
    bytes: int
    contigs: tuple[str, ...]
    lengths: tuple[int, ...]
    n_bases: int

    @property
    def length(self) -> int:
        """The genome's length in bases, N included."""
        return sum(self.lengths)


@dataclass(frozen=True)
class Annotation:
    """The annotation's lines on the genome's contigs, as written."""

    sha256: str
    release: str | None
    genes: int
    lines: int


@dataclass(frozen=True)
class StartCodons:
    """The annotation's 3-nt start codons read from the genome, and the split ones skipped."""

    records: int
    atg: int
    split: int
    codons: dict[str, int]

    @property
    def share(self) -> float | None:
        """ATG over all 3-nt records; None without any."""
        return self.atg / self.records if self.records else None


@dataclass(frozen=True)
class Reference:
    """A prepared reference folder and what was measured while it was built."""

    folder: Path
    built: bool
    key: str
    half: bool
    source: Genome
    annotation_source_sha256: str
    genome: Genome
    annotation: Annotation
    start_codons: StartCodons

    @property
    def fasta(self) -> Path:
        """The uncompressed genome the indices are built from."""
        return self.folder / GENOME

    @property
    def gtf(self) -> Path:
        """The annotation's lines on the genome's contigs."""
        return self.folder / ANNOTATION


def sha256(path: Path) -> str:
    """SHA-256 of a file as it is on disk."""
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _open(path: Path) -> IO[str]:
    """A text handle on a file, gzipped or not, told by its first bytes."""
    with path.open("rb") as handle:
        gzipped = handle.read(2) == GZIP_MAGIC
    return gzip.open(path, "rt") if gzipped else path.open()


def records(path: Path) -> Iterator[tuple[str, str]]:
    """Each FASTA record: its name, the header's first word, and its sequence.

    Raises:
        ProcessError: If a sequence comes before any header, a header has no name, or a name
            repeats.
    """
    name: str | None = None
    seen: set[str] = set()
    chunks: list[str] = []
    with _open(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                words = line[1:].split()
                if not words or words[0] in seen:
                    raise ProcessError(f"{path}:{number}: a FASTA header with no name or a "
                                       f"repeated one: {line.strip()!r}")  # fmt: skip
                name, chunks = words[0], []
                seen.add(name)
            elif line.strip():
                if name is None:
                    raise ProcessError(f"{path}:{number}: sequence before the first header")
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def genome_facts(path: Path) -> Genome:
    """Contig names and lengths, N bases, SHA-256 and size of a FASTA, plain or gzipped.

    Raises:
        ProcessError: If the file holds no FASTA record.
    """
    names, lengths, n = [], [], 0
    for name, sequence in records(path):
        names.append(name)
        lengths.append(len(sequence))
        n += sequence.count("N") + sequence.count("n")
    if not names:
        raise ProcessError(f"{path} holds no FASTA record")
    return Genome(sha256(path), path.stat().st_size, tuple(names), tuple(lengths), n)


def write_genome(source: Path, dest: Path, *, half: bool) -> None:
    """The genome uncompressed, each contig whole or its first floor(length / 2) bases."""
    with dest.open("w") as out:
        for name, sequence in records(source):
            kept = sequence[: len(sequence) // 2] if half else sequence
            out.write(f">{name}\n")
            for start in range(0, len(kept), FASTA_WIDTH):
                out.write(kept[start : start + FASTA_WIDTH] + "\n")


def sa_index_bases(length: int) -> int:
    """STAR's --genomeSAindexNbases for a genome: min(14, floor(log2(length) / 2 - 1)).

    7 for the 101,000-base fixture slice, 11 for chr22 (50,818,468 bases).
    """
    return min(SA_INDEX_CAP, math.floor(math.log2(length) / 2 - 1))


def annotation_subset(
    gtf: Path, contigs: Collection[str], dest: Path, cut: Mapping[str, int] | None = None
) -> Annotation:
    """Write the GTF's lines whose seqid is one of the contigs, exactly, with its ## header.

    Args:
        gtf: The annotation, plain or gzipped.
        contigs: The genome's contig names.
        dest: The subset GTF to write.
        cut: Each contig's kept length; a gene stays only when its gene record ends inside it,
            and every line of a gene goes with it.

    Returns:
        The subset's SHA-256, release, gene count and line count.

    Raises:
        ProcessError: If no line falls on the contigs (the first names of each listed), or, with
            a cut, a feature comes before its gene record.
    """
    wanted, keep = set(contigs), dict[str, bool]()
    seqids: dict[str, None] = {}
    release, genes, lines = None, 0, 0
    with _open(gtf) as source, dest.open("w") as out:
        for line in source:
            if line.startswith("#"):
                if line.startswith("##"):
                    out.write(line)
                    found = RELEASE.search(line)
                    release = release or (f"GENCODE {found[1]} (Ensembl {found[2]})" if found
                                          else None)  # fmt: skip
                continue
            seqid = line[: line.find("\t")]
            if len(seqids) < SHOWN:
                seqids.setdefault(seqid)
            if seqid not in wanted:
                continue
            fields = line.split("\t", 8)
            if cut is not None and not _inside(fields, cut[seqid], keep, gtf):
                continue
            out.write(line)
            lines += 1
            genes += fields[2] == "gene"
    if not lines:
        raise ProcessError(
            f"the genome and {gtf.name} share no contig: the genome has "
            f"{sorted(wanted)[:SHOWN]}, the GTF {list(seqids)}"
        )
    return Annotation(sha256(dest), release, genes, lines)


def _inside(fields: list[str], length: int, keep: dict[str, bool], gtf: Path) -> bool:
    """Whether a GTF line's gene ends within the kept length; the gene record decides."""
    found = GENE_ID.search(fields[8])
    if found is None:
        raise ProcessError(f"{gtf.name}: a line with no gene_id: {fields[:5]}")
    if fields[2] == "gene":
        keep[found[1]] = int(fields[4]) <= length
    elif found[1] not in keep:
        raise ProcessError(f"{gtf.name}: {found[1]} has a {fields[2]} before its gene record")
    return keep[found[1]]


def start_codon_share(fasta: Path, gtf: Path) -> StartCodons:
    """Read each 3-nt start_codon record from the genome, on its strand.

    A minus-strand codon is reverse-complemented. Records split across an intron (two
    records of 2 and 1 nt) are counted as split and not read.
    """
    wanted: dict[str, list[tuple[int, int, str]]] = {}
    split = 0
    with _open(gtf) as handle:
        for line in handle:
            fields = line.split("\t", 8)
            if len(fields) < 9 or fields[2] != "start_codon":  # noqa: PLR2004 - GTF's 9 columns
                continue
            start, end = int(fields[3]), int(fields[4])
            if end - start != 2:  # noqa: PLR2004 - a whole codon spans 3 bases
                split += 1
                continue
            wanted.setdefault(fields[0], []).append((start, end, fields[6]))
    codons: Counter[str] = Counter()
    for name, sequence in records(fasta):
        for start, end, strand in wanted.get(name, []):
            codon = sequence[start - 1 : end].upper()
            codons[codon.translate(COMPLEMENT)[::-1] if strand == "-" else codon] += 1
    return StartCodons(sum(codons.values()), codons["ATG"], split, dict(codons.most_common()))


def prepare(genome: Path, gtf: Path, parent: Path, *, half: bool = False) -> Reference:
    """The reference folder for a genome and annotation: reused when intact, else built.

    Args:
        genome: The genome FASTA, plain or gzipped.
        gtf: The annotation, plain or gzipped.
        parent: Where keyed folders sit (``tools/index/``).
        half: Keep each contig's first half and the genes that end inside it.

    Returns:
        The folder and its facts.
    """
    annotation_sha256 = sha256(gtf)
    inputs = {"genome_sha256": sha256(genome), "annotation_sha256": annotation_sha256, "half": half}

    def build(dest: Path) -> None:
        """Write the genome copy, the annotation subset and what was measured on them."""
        source = genome_facts(genome)
        write_genome(genome, dest / GENOME, half=half)
        annotation = annotation_subset(gtf, source.contigs, dest / ANNOTATION,
                                       halves(source) if half else None)  # fmt: skip
        codons = start_codon_share(dest / GENOME, dest / ANNOTATION)
        facts = {"source": asdict(source), "genome": asdict(genome_facts(dest / GENOME)),
                 "annotation": asdict(annotation), "start_codons": asdict(codons)}  # fmt: skip
        (dest / FACTS).write_text(json.dumps(facts, indent=1) + "\n")

    folder, built = cache.keyed(parent, "reference", inputs, build)
    facts = json.loads((folder / FACTS).read_text())
    return Reference(
        folder=folder, built=built, key=cache.key(inputs), half=half,
        source=_genome(facts["source"]), annotation_source_sha256=annotation_sha256,
        genome=_genome(facts["genome"]), annotation=Annotation(**facts["annotation"]),
        start_codons=StartCodons(**facts["start_codons"]),
    )  # fmt: skip


def halves(source: Genome) -> dict[str, int]:
    """Each contig's kept length on the half reference: floor(length / 2)."""
    return {c: n // 2 for c, n in zip(source.contigs, source.lengths, strict=True)}


def from_facts(prepared: Path, gtf: Path, parent: Path) -> Reference:
    """A reference ``prepare`` built elsewhere, from its facts, its annotation rebuilt here.

    The cluster keeps the genome; its folder's ``reference.json`` and ``ctx_manifest.json``
    are copied (ED.6 block D). The annotation lines are written again from ``gtf`` by
    ``annotation_subset``, the code that wrote them there, and must match their SHA-256.

    Raises:
        ProcessError: If the facts or the manifest are absent, ``gtf`` is not the annotation
            the cluster read, or the rebuilt lines differ from the cluster's.
    """
    try:
        facts = json.loads((prepared / FACTS).read_text())
        given = json.loads((prepared / cache.MANIFEST).read_text())["inputs"]
    except (OSError, json.JSONDecodeError, KeyError) as error:
        raise ProcessError(f"{prepared} holds no readable {FACTS} and {cache.MANIFEST}: "
                           f"{error}") from error  # fmt: skip
    annotation_sha256 = sha256(gtf)
    if given.get("annotation_sha256") != annotation_sha256:
        raise ProcessError(
            f"{prepared}: the cluster read an annotation of SHA-256 "
            f"{given.get('annotation_sha256')}, {gtf.name} has {annotation_sha256}"
        )
    source, half = _genome(facts["source"]), bool(given.get("half"))

    def build(dest: Path) -> None:
        """Write the annotation lines on the cluster genome's contigs, and its facts."""
        annotation_subset(gtf, source.contigs, dest / ANNOTATION, halves(source) if half else None)
        (dest / FACTS).write_text(json.dumps(facts, indent=1) + "\n")

    folder, built = cache.keyed(parent, "reference_facts", given, build)
    if (rebuilt := sha256(folder / ANNOTATION)) != facts["annotation"]["sha256"]:
        raise ProcessError(f"{folder / ANNOTATION} has SHA-256 {rebuilt}, the cluster's lines "
                           f"{facts['annotation']['sha256']}")  # fmt: skip
    return Reference(
        folder=folder, built=built, key=cache.key(given), half=half, source=source,
        annotation_source_sha256=annotation_sha256, genome=_genome(facts["genome"]),
        annotation=Annotation(**facts["annotation"]),
        start_codons=StartCodons(**facts["start_codons"]),
    )  # fmt: skip


def _genome(raw: dict[str, Any]) -> Genome:
    """A Genome from its JSON form, where tuples came back as lists."""
    return Genome(raw["sha256"], raw["bytes"], tuple(raw["contigs"]), tuple(raw["lengths"]),
                  raw["n_bases"])  # fmt: skip


def star_index(reference: Reference, parent: Path, threads: int) -> tuple[Path, bool]:
    """STAR's genome index of a prepared reference, with its annotation's splice junctions.

    Returns:
        The index folder, and whether it was built now.
    """
    bases = sa_index_bases(reference.genome.length)
    inputs = {
        "tool": "STAR", "version": tools.TOOLS["STAR"].version, "half": reference.half,
        "genome_sha256": reference.genome.sha256, "annotation_sha256": reference.annotation.sha256,
        "sa_index_bases": bases, "overhang": OVERHANG,
    }  # fmt: skip

    def build(dest: Path) -> None:
        """Run STAR's genomeGenerate into the folder."""
        tools.run("STAR", [
            "--runMode", "genomeGenerate", "--runThreadN", str(threads), "--genomeDir", str(dest),
            "--genomeFastaFiles", str(reference.fasta), "--sjdbGTFfile", str(reference.gtf),
            "--sjdbOverhang", str(OVERHANG), "--genomeSAindexNbases", str(bases),
            "--outFileNamePrefix", f"{dest}/",
        ], cwd=dest, log=dest / "genomeGenerate.log")  # fmt: skip

    return cache.keyed(parent, "star", inputs, build)


def hisat2_index(reference: Reference, parent: Path, threads: int) -> tuple[Path, bool]:
    """HISAT2's plain genome index (no splice sites or exons) of a prepared reference.

    Returns:
        The index folder, whose index prefix is ``<folder>/genome``, and whether it was built
        now.
    """
    inputs = {"tool": "hisat2-build", "version": tools.TOOLS["hisat2-build"].version,
              "half": reference.half, "genome_sha256": reference.genome.sha256}  # fmt: skip

    def build(dest: Path) -> None:
        """Run hisat2-build into the folder."""
        tools.run("hisat2-build", ["-p", str(threads), str(reference.fasta),
                                   str(dest / HISAT2_PREFIX)],
                  cwd=dest, log=dest / "hisat2-build.log")  # fmt: skip

    return cache.keyed(parent, "hisat2", inputs, build)
