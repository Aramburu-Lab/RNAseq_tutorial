"""M3's reads to alignments: the cap, fastp, the aligner, samtools (D-036, D-038).

Reads enter per run, as the gather sheet names them (D-031). The cap and fastp act per run;
the aligner acts per sample, its runs comma-joined, so a sample split over two runs is one
alignment. Each step writes one content-keyed folder (``ctx.process.cache``, D-037), keyed by
the tool, its pinned version, every option that can change its output and the SHA-256 of
every file it reads; threads and paths are not keyed.

STAR is used only after it has mapped reads here: the pinned osx-arm64 build passes
``--version`` and exits 0 while reading no reads (D-036), so ``star_probe`` maps 10 reads
through the production options before any alignment, and nothing falls back to HISAT2 on its
own.
"""

import ast
import functools
import gzip
import math
import platform
import random
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from ctx import tools
from ctx.process import STRAND, Aligner, ProcessError, cache
from ctx.process.reference import HISAT2_PREFIX, OVERHANG, sa_index_bases, sha256
from ctx.qc import Fastp
from ctx.qc import multiqc_io as mqc_io
from ctx.qc.reads import run_fastp, trimmed_names

SEED = 42
# seqkit sample -p keeps each read with probability p, so it returns a binomial number of reads.
# p sits 5 standard deviations above the cap, so `seqkit head -n cap` finds the cap but for a
# chance of about 3e-7 (seed 42 at p = 100 / 1,000 keeps 90 of se_a's reads, measured).
CAP_SD = 5
GZIP_MAGIC, CHUNK = b"\x1f\x8b", 1 << 24
CAPPED, REPORTS, READS, LOGS = "capped", "fastp", "reads", "logs"
# Production STAR options: gzip -dcf reads plain FASTQ too (-f); unmapped reads stay in the BAM,
# so flagstat counts every read (I6). MultiQC names a STAR sample by its file prefix.
STAR_OPTIONS = (
    "--readFilesCommand", "gzip", "-dcf", "--quantMode", "GeneCounts",
    "--outSAMtype", "BAM", "SortedByCoordinate", "--outSAMunmapped", "Within",
)  # fmt: skip
# HISAT2's options besides threads, paths and the unspliced switch; keyed like STAR's. Unaligned
# reads stay in the SAM by default (no --no-unal), as with STAR's Within.
HISAT2_OPTIONS = ("--new-summary",)
# HISAT2's --new-summary for pairs: pairs aligned concordantly once or more, or discordantly
# once, hold two mapped mates; mates of its other pairs aligned alone count one each.
HISAT2_PAIRS = ("paired_aligned_one", "paired_aligned_multi", "paired_aligned_discord_one")
HISAT2_ALONE = ("unpaired_aligned_one", "unpaired_aligned_multi")
STAR_FILES = ("Log.final.out", "ReadsPerGene.out.tab", "Aligned.sortedByCoord.out.bam")
STAR_BAM = STAR_FILES[2]
# ReadsPerGene.out.tab: gene id, then STAR columns 2 to 4, named by the strand each counts.
RPG_COLUMNS = ["gene_id", *(s.value for s, _ in sorted(STRAND.items(), key=lambda i: i[1][0]))]
# The probe: a random genome of two genes on opposite strands, and reads copied from their exons.
# Two exons each: with one-exon genes linux STAR 2.7.11b aborted building this genome ("double
# free or corruption", SAindexNbases 4, linux-64 and aarch64), while two-exon genes built on
# every size tried; the cause is not traced (D-036).
PROBE_LENGTH, PROBE_READS, PROBE_READ_LENGTH = 2000, 10, 50
PROBE_GENES = (("probe_plus", "+", ((101, 450), (551, 900))),
               ("probe_minus", "-", ((1101, 1450), (1551, 1900))))  # fmt: skip
D036 = (
    "D-036: STAR is used only once it maps the probe's 10 reads. Run with --aligner hisat2, or "
    "set CTX_TOOLS to a folder holding a STAR 2.7.11b that maps."
)
D036_MAC = (
    "On darwin-arm64 this is the known fault of the pinned osx-arm64 build: it passes --version "
    "yet exits 0 reading no reads (STAR issue #2632), exits 102 spawning --readFilesCommand, "
    "and with --quantMode GeneCounts exits 109, having lost the path of geneInfo.tab."
)


@dataclass(frozen=True)
class Trimmed:
    """A keyed trim folder: each run's capped and trimmed reads, and fastp's reports."""

    folder: Path
    built: bool
    given: Mapping[str, tuple[Path, ...]]

    @property
    def reports(self) -> Path:
        """The folder of fastp's JSON reports, one per run, named by the run."""
        return self.folder / REPORTS

    def reads(self, run: str) -> tuple[Path, ...]:
        """The run's trimmed reads, read 1 first."""
        names = trimmed_names(run, len(self.given[run]))
        return tuple(self.folder / READS / name for name in names)

    def raw(self, run: str) -> tuple[Path, ...]:
        """What fastp read for the run: the capped copy when the cap cut it, else its FASTQs."""
        capped = tuple(self.folder / CAPPED / n for n in trimmed_names(run, len(self.given[run])))
        return capped if all(p.is_file() for p in capped) else self.given[run]


@dataclass(frozen=True)
class AlignOptions:
    """One alignment setting.

    ``multimap`` is STAR's --outFilterMultimapNmax (None leaves STAR's 10); ``spliced`` False
    adds HISAT2's --no-spliced-alignment. ``index`` is the index folder. A setting the aligner
    has no switch for is refused (``align``), never keyed as if it took effect.
    """

    aligner: Aligner
    index: Path
    threads: int
    multimap: int | None = None
    spliced: bool = True


@dataclass(frozen=True)
class Alignment:
    """A keyed alignment folder: one sorted, indexed BAM per sample and the logs MultiQC reads."""

    folder: Path
    built: bool
    aligner: Aligner
    samples: tuple[str, ...]
    paired: bool

    def bam(self, sample: str) -> Path:
        """The sample's coordinate-sorted BAM, unmapped reads included."""
        return self.folder / bam_name(sample)

    def logs(self) -> list[Path]:
        """What MultiQC reads: the aligner's summaries and samtools' flagstat and idxstats.

        Never the aligner's own log: hisat2 prints its summary to stderr too, and MultiQC
        would read the sample twice (measured in M3.2).
        """
        own = STAR_FILES[:2] if self.aligner is Aligner.STAR else ("hisat2.txt",)
        return [self.folder / f"{s}.{name}" for s in self.samples
                for name in (*own, "flagstat", "idxstats")]  # fmt: skip


def bam_name(sample: str) -> str:
    """A sample's BAM in its alignment folder; featureCounts heads its matrix column with it."""
    return f"{sample}.bam"


def runs(files: pl.DataFrame) -> dict[str, tuple[Path, ...]]:
    """Each run's FASTQ files, read 1 first, in sheet order (``reads.fastq_files``' table)."""
    found: dict[str, list[tuple[int, Path]]] = {}
    for run, mate, path in files.select("run_accession", "mate", "path").rows():
        found.setdefault(run, []).append((mate, Path(path)))
    return {run: tuple(p for _, p in sorted(mates)) for run, mates in found.items()}


def by_sample(files: pl.DataFrame, reads: Mapping[str, Sequence[Path]]
              ) -> dict[str, list[tuple[Path, ...]]]:  # fmt: skip
    """Each sample's runs in sheet order, each run's files from ``reads``: the aligner's input."""
    out: dict[str, list[tuple[Path, ...]]] = {}
    pairs = files.unique("run_accession", keep="first", maintain_order=True)
    for run, sample in pairs.select("run_accession", "sample").rows():
        out.setdefault(sample, []).append(tuple(reads[run]))
    return out


def proportion(cap: int, total: int) -> float:
    """The `seqkit sample -p` proportion that yields at least ``cap`` of ``total`` reads."""
    return min(1.0, (cap + CAP_SD * math.sqrt(cap)) / total)


def count_reads(path: Path) -> int:
    """Reads in a FASTQ, plain or gzipped, as its lines over 4.

    The cap's proportion needs the count before fastp runs. The FASTQ is read as data, so no
    tool's log is parsed outside MultiQC (D-038); the realised count comes from fastp's report
    through MultiQC (``fastp_reads``).

    Raises:
        ProcessError: If the lines are not a multiple of 4: not a four-line FASTQ.
    """
    with path.open("rb") as raw:
        packed = raw.read(len(GZIP_MAGIC)) == GZIP_MAGIC
    lines, last = 0, b"\n"
    with gzip.open(path, "rb") if packed else path.open("rb") as handle:
        while chunk := handle.read(CHUNK):
            lines, last = lines + chunk.count(b"\n"), chunk[-1:]
    lines += last != b"\n"
    if lines % 4:
        raise ProcessError(f"{path} holds {lines} lines, not 4 per read")
    return lines // 4


def cap_run(paths: Sequence[Path], cap: int, dest: Path, title: str) -> tuple[Path, ...]:
    """At most ``cap`` reads of a run, the same from each mate: a seeded sample cut at the cap.

    `seqkit sample -p` keeps each read with probability p and seed 42; `seqkit head -n cap`
    then keeps the first ``cap`` of them in file order, as seqkit's help advises for large
    files. p sits CAP_SD standard deviations above cap / reads, so the cut drops sampled reads
    only near the end of the file, where reads are kept less often. When p reaches 1 every
    read is sampled and the cap is the run's first ``cap`` reads. Both mates have one read
    count and one seed, so they keep the same reads.

    Returns:
        The capped files, or the run's own files when it holds no more than ``cap`` reads.
    """
    # Absolute: seqkit runs in the log folder, where a relative -o lands elsewhere, exit 0.
    dest = dest.resolve()
    logs = dest / LOGS
    logs.mkdir(parents=True, exist_ok=True)
    total = count_reads(Path(paths[0]))
    if total <= cap:
        return tuple(paths)
    out = []
    for path, name in zip(paths, trimmed_names(title, len(paths)), strict=True):
        sampled = logs / f"{name}.sampled.fastq.gz"
        tools.run("seqkit", ["sample", "-p", repr(proportion(cap, total)), "-s", str(SEED),
                             str(Path(path).resolve()), "-o", str(sampled)],
                  cwd=logs, log=logs / f"{name}.sample.log")  # fmt: skip
        tools.run("seqkit", ["head", "-n", str(cap), str(sampled), "-o", str(dest / name)],
                  cwd=logs, log=logs / f"{name}.head.log")  # fmt: skip
        # A scratch file of this build, never anyone's data: the capped copy replaces it.
        sampled.unlink()
        out.append(dest / name)
    return tuple(out)


def trim(files: pl.DataFrame, parent: Path, *, fastp: Fastp, cap: int | None,
         threads: int) -> Trimmed:  # fmt: skip
    """Every run capped (when ``cap`` is set) and trimmed by fastp, in one keyed folder.

    Args:
        files: ``reads.fastq_files``' first table.
        parent: Where keyed folders sit (``tools/<dataset>/``).
        fastp: The filter, M2's own (D-034), so M2 and M3 share one fastp setting.
        cap: Reads kept per run before fastp; None keeps all.
        threads: fastp's worker threads.

    Returns:
        The folder, and where each run's reads and report are.
    """
    given = runs(files)
    inputs = {"fastp": tools.TOOLS["fastp"].version, "seqkit": tools.TOOLS["seqkit"].version,
              "fastp_args": fastp.args(), "cap": cap, "seed": SEED,
              "runs": {r: [sha256(p) for p in paths] for r, paths in given.items()}}  # fmt: skip

    def build(dest: Path) -> None:
        """Cap and trim each run."""
        (dest / REPORTS).mkdir()
        for run, paths in given.items():
            read = cap_run(paths, cap, dest / CAPPED, run) if cap is not None else paths
            run_fastp([str(p) for p in read], dest / REPORTS, run, fastp, threads,
                      trimmed=dest / READS)  # fmt: skip

    folder, built = cache.keyed(parent, "trim", inputs, build)
    return Trimmed(folder, built, given)


def _mates(reads: Mapping[str, Sequence[Sequence[Path]]]) -> int:
    """The one number of files per run across every sample.

    Raises:
        ProcessError: If single-end and paired-end runs are mixed; each kind's samples listed.
    """
    kinds: dict[int, list[str]] = {}
    for sample, sample_runs in reads.items():
        for files in sample_runs:
            kinds.setdefault(len(files), []).append(sample)
    if len(kinds) != 1 or not set(kinds) <= {1, 2}:
        raise ProcessError(f"runs of one alignment must be all single-end or all paired-end; "
                           f"samples by files per run: {kinds}")  # fmt: skip
    return next(iter(kinds))


def _joined(sample_runs: Sequence[Sequence[Path]]) -> list[str]:
    """Each mate's files across a sample's runs, comma-joined and absolute, read 1 first."""
    return [",".join(str(Path(p).resolve()) for p in files)
            for files in zip(*sample_runs, strict=True)]  # fmt: skip


def star_align(sample: str, sample_runs: Sequence[Sequence[Path]], index: Path, out: Path,
               multimap: int | None, threads: int) -> Path:  # fmt: skip
    """One sample through STAR with the production options; its runs comma-joined.

    Returns:
        The BAM, renamed from STAR's name to ``<sample>.bam``.

    Raises:
        ToolError: If STAR fails, or exits 0 without writing one of its outputs.
    """
    out = out.resolve()
    (out / LOGS).mkdir(parents=True, exist_ok=True)
    mates = _joined(sample_runs)
    args = ["--runThreadN", str(threads), "--genomeDir", str(index.resolve()), "--readFilesIn",
            *mates, *STAR_OPTIONS, "--outFileNamePrefix", f"{out}/{sample}."]  # fmt: skip
    args += ["--outFilterMultimapNmax", str(multimap)] if multimap is not None else []
    tools.run("STAR", args, cwd=out, log=out / LOGS / f"{sample}.star.log")
    if missing := [n for n in STAR_FILES if not (out / f"{sample}.{n}").is_file()]:
        raise tools.ToolError(f"STAR exited 0 and wrote no {missing} for {sample} in {out}")
    return (out / f"{sample}.{STAR_BAM}").rename(out / bam_name(sample))


def hisat2_align(sample: str, sample_runs: Sequence[Sequence[Path]], options: AlignOptions,
                 out: Path) -> Path:  # fmt: skip
    """One sample through HISAT2, its runs comma-joined, sorted into ``<sample>.bam``.

    No --rna-strandness: the strand is called from the counts downstream (D-039), and the
    alignment must not assume it. Unaligned reads stay in the BAM, as with STAR's ``Within``.
    """
    out = out.resolve()
    (out / LOGS).mkdir(parents=True, exist_ok=True)
    mates = _joined(sample_runs)
    given = ["-U", mates[0]] if len(mates) == 1 else ["-1", mates[0], "-2", mates[1]]
    sam = out / f"{sample}.sam"
    tools.run("hisat2", [
        "-p", str(options.threads), "-x", str(options.index.resolve() / HISAT2_PREFIX), *given,
        *HISAT2_OPTIONS, "--summary-file", str(out / f"{sample}.hisat2.txt"),
        *([] if options.spliced else ["--no-spliced-alignment"]), "-S", str(sam),
    ], cwd=out, log=out / LOGS / f"{sample}.hisat2.log")  # fmt: skip
    bam = out / bam_name(sample)
    tools.run("samtools", ["sort", "-@", str(options.threads), "-o", str(bam), str(sam)],
              cwd=out, log=out / LOGS / f"{sample}.sort.log")  # fmt: skip
    # The unsorted SAM is this build's scratch; the sorted BAM holds every record of it.
    sam.unlink()
    return bam


def bam_stats(bam: Path) -> None:
    """Index a BAM and write its flagstat and idxstats beside it, named by its sample."""
    # Absolute: samtools runs in the BAM's folder, where a relative path misses (exit 1).
    bam = bam.resolve()
    out, stem = bam.parent, bam.name.removesuffix(".bam")
    (out / LOGS).mkdir(parents=True, exist_ok=True)
    tools.run("samtools", ["index", str(bam)], cwd=out, log=out / LOGS / f"{stem}.index.log")
    for command in ("flagstat", "idxstats"):
        tools.run("samtools", [command, str(bam)], cwd=out, stdout=out / f"{stem}.{command}",
                  log=out / LOGS / f"{stem}.{command}.log")  # fmt: skip


def align(reads: Mapping[str, Sequence[Sequence[Path]]], parent: Path,
          options: AlignOptions) -> Alignment:  # fmt: skip
    """Every sample aligned, sorted, indexed and summarised, in one keyed folder.

    Args:
        reads: Each sample's runs, each run's files read 1 first (``by_sample``).
        parent: Where keyed folders sit.
        options: The aligner and its setting.

    Returns:
        The folder and its samples, in the order given.

    Raises:
        ProcessError: If single-end and paired-end runs are mixed, or a setting has no switch in
            the aligner: STAR cannot turn off spliced alignment (M3 plan, section 2), and
            HISAT2 here takes no multimapper limit.
    """
    star = options.aligner is Aligner.STAR
    if star and not options.spliced:
        raise ProcessError("spliced=False is HISAT2's --no-spliced-alignment; STAR has no "
                           "switch that turns off spliced alignment")  # fmt: skip
    if not star and options.multimap is not None:
        raise ProcessError(f"multimap={options.multimap} is STAR's --outFilterMultimapNmax; "
                           "HISAT2 takes no multimapper limit here")  # fmt: skip
    mates = _mates(reads)
    tool = "STAR" if star else "hisat2"
    inputs = {
        "tool": tool, "version": tools.TOOLS[tool].version,
        "samtools": tools.TOOLS["samtools"].version, "index": options.index.resolve().name,
        "options": list(STAR_OPTIONS if star else HISAT2_OPTIONS), "multimap": options.multimap,
        "spliced": options.spliced,
        "reads": {s: [[sha256(Path(p)) for p in files] for files in rs] for s, rs in reads.items()},
    }  # fmt: skip

    def build(dest: Path) -> None:
        """Align each sample, then summarise its BAM."""
        for sample, sample_runs in reads.items():
            if star:
                bam = star_align(sample, sample_runs, options.index, dest, options.multimap,
                                 options.threads)  # fmt: skip
            else:
                bam = hisat2_align(sample, sample_runs, options, dest)
            bam_stats(bam)

    folder, built = cache.keyed(parent, f"align_{options.aligner}", inputs, build)
    return Alignment(folder, built, options.aligner, tuple(reads), mates == 2)  # noqa: PLR2004 - a pair's two files


def summarise(alignment: Alignment, dest: Path, config: Path,
              trimmed: Trimmed | None = None) -> mqc_io.MultiqcRun:  # fmt: skip
    """MultiQC on an alignment's logs, and on fastp's reports when ``trimmed`` is given.

    Args:
        alignment: Its aligner's summaries and samtools' flagstat and idxstats are read.
        dest: MultiQC's output folder.
        config: The MultiQC config; RAW_COUNTS keeps counts unscaled (D-029).
        trimmed: The trim folder whose fastp reports go into the same report.

    Returns:
        The run's data folder and MultiQC's version.
    """
    inputs = [*alignment.logs(), *([trimmed.reports] if trimmed else [])]
    modules = (str(alignment.aligner), "samtools", *(["fastp"] if trimmed else []))
    return mqc_io.run(inputs, dest, config=config, modules=modules)


def reads_per_gene(tab: Path) -> pl.DataFrame:
    """STAR's ReadsPerGene.out.tab, every row: the four N_ rows, then one row per gene.

    Columns: gene_id, then STAR columns 2 to 4 named by the strand they count (``STRAND``).
    """
    return pl.read_csv(tab, separator="\t", has_header=False, new_columns=RPG_COLUMNS,
                       schema_overrides={c: pl.Int64 for c in RPG_COLUMNS[1:]})  # fmt: skip


def derived_limit_one(star: Mapping[str, float]) -> dict[str, float]:
    """A STAR sample's Log.final.out numbers had --outFilterMultimapNmax been 1.

    A read with more than one locus becomes "mapped to too many loci"; unique reads stay
    unique, and no gene count changes, because GeneCounts never counts a multimapper.

    Args:
        star: One row of MultiQC's STAR table, as numbers.

    Returns:
        The row with multimapped moved into multimapped_toomany, percentages recomputed as
        STAR prints them: reads / total * 100 to two decimals. 100 * 49 / 160 rounds to 30.62,
        STAR's 49 / 160 * 100 to 30.63 (measured).
    """
    total = star["total_reads"]
    moved = {"multimapped": 0.0, "multimapped_toomany": star["multimapped"]
             + star["multimapped_toomany"], "mapped": star["uniquely_mapped"]}  # fmt: skip
    percents = {f"{k}_percent": float(f"{v / total * 100:.2f}") for k, v in moved.items()}
    return {**star, **moved, **percents}


def fastp_reads(run: mqc_io.MultiqcRun, names: Sequence[str]) -> pl.DataFrame:
    """Each run's reads before and after fastp, from MultiQC's fastp table, mates counted apart.

    MultiQC 1.35 keeps fastp's summary and filtering result as Python literals in one cell
    each, so their integers are read exactly.

    Raises:
        QcError: If MultiQC wrote no fastp table, or a row names no run.
        ProcessError: If a run has no fastp row.
    """
    table = mqc_io.table(run, "fastp")
    mqc_io.map_names(table["Sample"], {n: n for n in names})
    if missing := sorted(set(names) - set(table["Sample"])):
        raise ProcessError(f"no fastp report for {missing}")
    cells = table.select("Sample", "summary", "filtering_result").rows()
    rows = [(name, ast.literal_eval(summary)["before_filtering"]["total_reads"],
             ast.literal_eval(result)["passed_filter_reads"])
            for name, summary, result in cells]  # fmt: skip
    return pl.DataFrame(rows, schema={"run_accession": pl.String, "raw_reads": pl.Int64,
                                      "passed_reads": pl.Int64}, orient="row")  # fmt: skip


def aligner_reads(run: mqc_io.MultiqcRun, alignment: Alignment) -> pl.DataFrame:
    """Each sample's input reads, mapped reads and mapped mates as its aligner reports them.

    A pair is one read to both aligners. Mapped is unique plus multimapped reads: pairs under
    STAR, which counts a pair mapped when one mate maps; null for HISAT2 pairs, whose summary
    counts mates. Mapped mates are what flagstat should count as mapped records (I6): STAR's
    mapped reads, twice for pairs; under HISAT2, twice the pairs aligned concordantly or
    discordantly, plus the mates of the other pairs aligned alone.

    Raises:
        ProcessError: If a sample has no row, or a column is absent: a zero from a column the
            table does not carry is not a result.
    """
    mapped: tuple[str, ...]
    pairs: tuple[str, ...]
    if alignment.aligner is Aligner.STAR:
        given, mapped, pairs = "total_reads", ("uniquely_mapped", "multimapped"), ()
    elif alignment.paired:
        given, mapped, pairs = "paired_total", (), HISAT2_PAIRS
    else:
        given, mapped, pairs = "unpaired_total", HISAT2_ALONE, ()
    columns = {given, *mapped, *pairs, *(HISAT2_ALONE if pairs else ())}
    table = mqc_io.table(run, str(alignment.aligner))
    if missing := sorted(columns - set(table.columns)):
        raise ProcessError(f"MultiQC's {alignment.aligner} table has no {missing}")
    if absent := sorted(set(alignment.samples) - set(table["Sample"])):
        raise ProcessError(f"no {alignment.aligner} report for {absent}")
    table = mqc_io.numeric(table, columns).filter(pl.col("Sample").is_in(alignment.samples))
    mates = (2 * pl.sum_horizontal(pairs) + pl.sum_horizontal(HISAT2_ALONE) if pairs
             else (2 if alignment.paired else 1) * pl.sum_horizontal(mapped))  # fmt: skip
    return table.select(
        sample=pl.col("Sample"), input_reads=pl.col(given).cast(pl.Int64),
        mapped_reads=(pl.sum_horizontal(mapped) if mapped else pl.lit(None)).cast(pl.Int64),
        mapped_mates=mates.cast(pl.Int64),
    )  # fmt: skip


def _probe_inputs(work: Path) -> tuple[Path, Path, Path]:
    """The probe's genome, annotation and gzipped reads, from seed 42 and the standard library."""
    rng = random.Random(SEED)
    genome = "".join(rng.choices("ACGT", k=PROBE_LENGTH))
    (work / "probe.fa").write_text(f">probe\n{genome}\n")
    lines = []
    for gene, strand, exons in PROBE_GENES:
        ids = f'gene_id "{gene}"; transcript_id "{gene}.1";'
        spans = [("gene", exons[0][0], exons[-1][1]), ("transcript", exons[0][0], exons[-1][1]),
                 *(("exon", a, b) for a, b in exons)]  # fmt: skip
        lines += [f"probe\tctx\t{kind}\t{a}\t{b}\t.\t{strand}\t.\t{ids}\n" for kind, a, b in spans]
    (work / "probe.gtf").write_text("".join(lines))
    records = []
    for n in range(PROBE_READS):
        exons = PROBE_GENES[n % len(PROBE_GENES)][2]
        start, end = exons[n // len(PROBE_GENES) % len(exons)]
        at = rng.randrange(start - 1, end - PROBE_READ_LENGTH)
        read = genome[at : at + PROBE_READ_LENGTH]
        records.append(f"@probe{n}\n{read}\n+\n{'I' * len(read)}\n")
    with gzip.open(work / "probe.fastq.gz", "wt") as handle:
        handle.write("".join(records))
    return work / "probe.fa", work / "probe.gtf", work / "probe.fastq.gz"


def refusal() -> str:
    """D-036's routes, and the osx-arm64 build's known faults where they apply: darwin-arm64."""
    mac = platform.system() == "Darwin" and platform.machine() == "arm64"
    return f"{D036} {D036_MAC}" if mac else D036


def _map_probe(work: Path) -> pl.DataFrame:
    """Index the probe's genome and map its reads in ``work``; its ReadsPerGene rows."""
    genome, gtf, reads = _probe_inputs(work)
    (index := work / "index").mkdir()
    tools.run("STAR", [
        "--runMode", "genomeGenerate", "--runThreadN", "1", "--genomeDir", str(index),
        "--genomeFastaFiles", str(genome), "--sjdbGTFfile", str(gtf),
        "--sjdbOverhang", str(OVERHANG), "--genomeSAindexNbases", str(sa_index_bases(PROBE_LENGTH)),
        "--outFileNamePrefix", f"{index}/",
    ], cwd=index, log=work / "genomeGenerate.log")  # fmt: skip
    star_align("probe", [(reads,)], index, work, None, 1)
    return reads_per_gene(work / "probe.ReadsPerGene.out.tab")


@functools.cache
def _probe(star: Path) -> str | None:
    """Why the STAR at ``star`` failed the probe, None when it mapped; once per executable.

    A failure is cached too, so a process probes once. A failing probe's files are kept, and
    the message names their folder; a passing probe leaves nothing behind.
    """
    work = Path(tempfile.mkdtemp(prefix="ctx_star_probe_")).resolve()
    kept = f"the probe's files are kept in {work}. {refusal()}"
    try:
        counts = _map_probe(work)
    except (tools.ToolError, pl.exceptions.PolarsError) as error:
        return f"{star} cannot map the probe's reads; {kept}\n{error}"
    totals = {c: int(counts[c].sum()) for c in RPG_COLUMNS[1:]}
    on_genes = int(counts.filter(~pl.col("gene_id").str.starts_with("N_"))[RPG_COLUMNS[1]].sum())
    if set(totals.values()) != {PROBE_READS} or on_genes != PROBE_READS:
        return (f"{star} exited 0 but did not map the probe's {PROBE_READS} reads: each "
                f"ReadsPerGene column sums to {totals}, {on_genes} on genes; {kept}")  # fmt: skip
    # The probe's own scratch, never anyone's data.
    shutil.rmtree(work)
    return None


def star_probe() -> None:
    """Refuse a STAR that cannot map: 10 exact reads through the production options.

    Each ReadsPerGene column counts every read once (genes and the N_ rows), so every column
    must sum to 10, and all 10 must fall on the two genes. Exit 0 and ``--version`` both pass
    on the failing osx-arm64 build, so neither is the check (D-036).

    Raises:
        ToolError: If STAR is missing, reports another version or does not map the reads;
            the message names D-036 and both routes.
    """
    try:
        star = tools.require("STAR")
    except tools.ToolError as error:
        raise tools.ToolError(f"{error}. {refusal()}") from error
    if (failure := _probe(star)) is not None:
        raise tools.ToolError(failure)
