"""M4's laptop path: a transcriptome from M3's reference, Salmon's index, Salmon on the reads.

The transcriptome is spliced from the reference folder ``ctx.process.reference.prepare``
writes (the genome copy and the annotation's lines on its contigs, whose start codons check
that both share one assembly), one record per transcript, exons in transcript order and a
minus-strand transcript reverse-complemented; its transcript-to-gene table comes from the
same lines. The index and each quantification are content-keyed folders
(``ctx.process.cache``, D-037) under ``outputs/04_quant/tools/``. The quantification folder
holds Salmon's per-sample folders as nf-core lays them out, so MultiQC and ``ctx.quant.outdir``
read it as they read copy Q1, beside ctx's transcript-to-gene table and the samples'
transcript TPM joined into one table: quant.sf's TPM is tximport's transcript TPM unchanged.

Salmon is used only after it has mapped reads here: ``salmon_probe`` quantifies 10 reads
copied from a random two-transcript index and refuses a Salmon that does not map all 10, as
``star_probe`` refuses a STAR that passes ``--version`` and maps nothing (D-036).
"""

import functools
import gzip
import random
import re
import shutil
import tempfile
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from ctx import tools
from ctx.process import cache
from ctx.process.reference import COMPLEMENT, FASTA_WIDTH, Reference, records, sha256
from ctx.quant import QuantError, SalmonOptions, outdir

TRANSCRIPTS, TX2GENE, DECOYS, GENTROME, INDEX = (
    "transcripts.fa", "tx2gene.tsv", "decoys.txt", "gentrome.fa", "index",
)  # fmt: skip
# nf-core/rnaseq 3.27.0's salmon index k-mer, Salmon's own default.
KMER = 31
ATTRIBUTE = re.compile(r'(\w+) "([^"]*)"')
SEED = 42
# The probe: two random transcripts and reads copied from them, all of which must map.
PROBE_LENGTH, PROBE_TRANSCRIPTS, PROBE_READS, PROBE_READ_LENGTH, PROBE_KMER = 1000, 2, 10, 50, 19
PROBE_FIX = (
    "Salmon is used only once it maps the probe's 10 reads (D-036's rule for STAR). Set "
    "CTX_TOOLS to a folder holding a Salmon 2.7.0 that maps, or run ctx quant on an nf-core "
    "outdir with --rnaseq-dir."
)


@dataclass(frozen=True)
class Quantified:
    """A keyed quantification folder: Salmon's folder per sample, tx2gene and transcript TPM."""

    index: Path
    index_built: bool
    folder: Path
    built: bool


def _attributes(field: str) -> dict[str, str]:
    """A GTF line's attributes, first value of each key."""
    found: dict[str, str] = {}
    for key, value in ATTRIBUTE.findall(field):
        found.setdefault(key, value)
    return found


def exons(gtf: Path) -> dict[str, tuple[str, str, str, str, list[tuple[int, int]]]]:
    """Each transcript's contig, strand, gene, gene name and exons (1-based, closed), in order.

    Raises:
        QuantError: If an exon carries no transcript_id or gene_id.
    """
    found: dict[str, tuple[str, str, str, str, list[tuple[int, int]]]] = {}
    with gtf.open() as handle:
        for line in handle:
            fields = line.split("\t")
            if line.startswith("#") or len(fields) < 9 or fields[2] != "exon":  # noqa: PLR2004
                continue
            attrs = _attributes(fields[8])
            if "transcript_id" not in attrs or "gene_id" not in attrs:
                raise QuantError(f"{gtf.name}: an exon with no transcript_id or gene_id")
            entry = found.setdefault(attrs["transcript_id"], (
                fields[0], fields[6], attrs["gene_id"], attrs.get("gene_name", attrs["gene_id"]),
                []))  # fmt: skip
            entry[4].append((int(fields[3]), int(fields[4])))
    return found


def spliced(genome: Path, gtf: Path) -> Iterator[tuple[str, str]]:
    """Each transcript's name and sequence, in annotation order, its exons joined 5' to 3'."""
    models = exons(gtf)
    by_contig: dict[str, list[str]] = {}
    for name, (contig, *_) in models.items():
        by_contig.setdefault(contig, []).append(name)
    sequences = {}
    for contig, sequence in records(genome):
        for name in by_contig.get(contig, []):
            _, strand, _, _, parts = models[name]
            joined = "".join(sequence[a - 1 : b] for a, b in sorted(parts))
            sequences[name] = joined.translate(COMPLEMENT)[::-1] if strand == "-" else joined
    for name in models:
        if name in sequences:
            yield name, sequences[name]


def write_transcriptome(reference: Reference, dest: Path) -> int:
    """Write ``transcripts.fa`` and ``tx2gene.tsv`` (nf-core's columns); return the records."""
    models = exons(reference.gtf)
    n = 0
    with (dest / TRANSCRIPTS).open("w") as fasta:
        for name, sequence in spliced(reference.fasta, reference.gtf):
            fasta.write(f">{name}\n")
            fasta.writelines(sequence[i : i + FASTA_WIDTH] + "\n"
                             for i in range(0, len(sequence), FASTA_WIDTH))  # fmt: skip
            n += 1
    rows = [(t, m[2], m[3]) for t, m in models.items()]
    pl.DataFrame(rows, schema=list(outdir.TX2GENE_COLUMNS), orient="row").write_csv(
        dest / TX2GENE, separator="\t")  # fmt: skip
    return n


def index(reference: Reference, parent: Path, threads: int, *, decoys: bool) -> tuple[Path, bool]:
    """Salmon's index of the reference's transcriptome, the genome's contigs as decoys or not.

    Returns:
        The keyed folder, holding ``index/``, ``transcripts.fa`` and ``tx2gene.tsv``, and
        whether it was built now.
    """
    inputs = {"tool": "salmon", "version": tools.TOOLS["salmon"].version, "kmer": KMER,
              "genome_sha256": reference.genome.sha256, "decoys": decoys,
              "annotation_sha256": reference.annotation.sha256}  # fmt: skip

    def build(dest: Path) -> None:
        """Splice the transcripts, then index them, with the genome appended where decoyed."""
        write_transcriptome(reference, dest)
        target = TRANSCRIPTS
        extra: list[str] = []
        if decoys:
            (dest / DECOYS).write_text("".join(f"{c}\n" for c in reference.genome.contigs))
            with (dest / GENTROME).open("wb") as out:
                for part in (dest / TRANSCRIPTS, reference.fasta):
                    with part.open("rb") as source:
                        shutil.copyfileobj(source, out)
            target, extra = GENTROME, ["-d", DECOYS]
        tools.run("salmon", ["index", "-t", target, *extra, "-i", INDEX, "-k", str(KMER),
                             "-p", str(threads)], cwd=dest, log=dest / "index.log")  # fmt: skip
        # The transcripts and genome both stay on disk; their concatenation is scratch.
        (dest / GENTROME).unlink(missing_ok=True)

    return cache.keyed(parent, "salmon_index", inputs, build)


def quant_args(options: SalmonOptions) -> list[str]:
    """Salmon quant's options besides index, reads, threads and output."""
    return ["-l", options.libtype, *(["--gcBias"] if options.gc_bias else []),
            *(["--seqBias"] if options.seq_bias else [])]  # fmt: skip


def quantify(reads: Mapping[str, Sequence[Path]], index_folder: tuple[Path, bool], parent: Path,
             options: SalmonOptions, threads: int) -> Quantified:  # fmt: skip
    """Salmon quant on each sample's reads, into one keyed folder, then the joined TPM table.

    Args:
        reads: Each sample's read files, single-end; a sample's runs are given together.
        index_folder: ``index``'s folder and whether it was built now.
        parent: Where the keyed folders sit (``tools/<dataset>/``).
        options: The setting.
        threads: Salmon's ``-p``.
    """
    folder_in, index_built = index_folder
    inputs = {"tool": "salmon", "version": tools.TOOLS["salmon"].version,
              "index": folder_in.name, "args": quant_args(options),
              "samples": {s: [sha256(p) for p in paths] for s, paths in reads.items()}}  # fmt: skip

    def build(dest: Path) -> None:
        """Quantify each sample into its own folder, then write tx2gene and transcript TPM."""
        for sample, paths in reads.items():
            tools.run("salmon", ["quant", "-i", str(folder_in / INDEX), *quant_args(options),
                                 "-r", *(str(p.resolve()) for p in paths), "-p", str(threads),
                                 "-o", sample], cwd=dest, log=dest / f"{sample}.log")  # fmt: skip
        shutil.copy(folder_in / TX2GENE, dest / outdir.TX2GENE)
        transcript_tpm(dest, list(reads))

    folder, built = cache.keyed(parent, "salmon_quant", inputs, build)
    return Quantified(folder_in, index_built, folder, built)


def transcript_tpm(folder: Path, samples: Sequence[str]) -> None:
    """``salmon.merged.transcript_tpm.tsv``: tx, gene_id and each sample's quant.sf TPM."""
    genes = outdir.tx2gene(folder).select(tx="transcript_id", gene_id="gene_id")
    frame = outdir.quant(folder / samples[0] / "quant.sf").select(tx="Name")
    for sample in samples:
        tpm = outdir.quant(folder / sample / "quant.sf").select("TPM")["TPM"]
        frame = frame.with_columns(tpm.alias(sample))
    frame.join(genes, on="tx", how="left", maintain_order="left").select(
        "tx", "gene_id", *samples).write_csv(folder / outdir.MERGED.format("transcript_tpm"),
                                            separator="\t")  # fmt: skip


def _probe_inputs(work: Path) -> None:
    """The probe's transcripts and gzipped reads, from seed 42 and the standard library."""
    rng = random.Random(SEED)
    transcripts = ["".join(rng.choices("ACGT", k=PROBE_LENGTH)) for _ in range(PROBE_TRANSCRIPTS)]
    (work / TRANSCRIPTS).write_text("".join(f">probe{i}\n{s}\n" for i, s in enumerate(transcripts)))
    lines = []
    for n in range(PROBE_READS):
        source = transcripts[n % PROBE_TRANSCRIPTS]
        at = rng.randrange(PROBE_LENGTH - PROBE_READ_LENGTH)
        lines.append(f"@r{n}\n{source[at : at + PROBE_READ_LENGTH]}\n+\n"
                     f"{'I' * PROBE_READ_LENGTH}\n")  # fmt: skip
    with gzip.open(work / "probe.fastq.gz", "wt") as handle:
        handle.write("".join(lines))


@functools.cache
def _probe(salmon: Path) -> str | None:
    """Why the Salmon at ``salmon`` failed the probe, None when it mapped; once per executable.

    A failing probe's files are kept, and the message names their folder.
    """
    work = Path(tempfile.mkdtemp(prefix="ctx_salmon_probe_")).resolve()
    kept = f"the probe's files are kept in {work}. {PROBE_FIX}"
    try:
        _probe_inputs(work)
        tools.run("salmon", ["index", "-t", TRANSCRIPTS, "-i", INDEX, "-k", str(PROBE_KMER),
                             "-p", "1"], cwd=work, log=work / "index.log")  # fmt: skip
        tools.run("salmon", ["quant", "-i", INDEX, "-l", "A", "-r", "probe.fastq.gz", "-p", "1",
                             "-o", "out"], cwd=work, log=work / "quant.log")  # fmt: skip
        mapped = round(outdir.quant(work / "out" / "quant.sf")["NumReads"].sum())
    except (tools.ToolError, OSError, pl.exceptions.PolarsError) as error:
        return f"{salmon} cannot quantify the probe's reads; {kept}\n{error}"
    if mapped != PROBE_READS:
        return (f"{salmon} exited 0 but assigned {mapped} of the probe's {PROBE_READS} reads; "
                f"{kept}")  # fmt: skip
    # The probe's own scratch, never anyone's data.
    shutil.rmtree(work)
    return None


def salmon_probe() -> None:
    """Refuse a Salmon that cannot map: 10 exact reads must all reach the two transcripts.

    Raises:
        ToolError: If Salmon is missing, reports another version or does not assign the 10
            reads; the message names the routes left.
    """
    try:
        salmon = tools.require("salmon")
    except tools.ToolError as error:
        raise tools.ToolError(f"{error}. {PROBE_FIX}") from error
    if (failure := _probe(salmon)) is not None:
        raise tools.ToolError(failure)
