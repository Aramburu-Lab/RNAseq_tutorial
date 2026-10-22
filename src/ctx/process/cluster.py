"""M3 from a cluster's alignment: ``ctx process --align-dir``, never the default (ED.6e-4, ED.6e-5).

Block D of ``docs/plans/2026-10-22_ed6_m2_m3_core_cluster_blocks.md`` copies each sample's
``runs/<sample>/{fastp,star,featurecounts,time}/`` files, ``tasks.tsv``, ``versions.txt`` and
``refs.txt`` into ``data/raw/<slug>_align/``. This module lays them out as the local route's keyed
folders, a ``Trimmed``, an ``Alignment`` and ``Counts``, so MultiQC, the checks I1 to I6, the
tables, the figures and the page read both routes through one code path.

featureCounts ran once per sample and setting: ``s<k>`` and ``O_s<k>``, ``O_fraction_s<k>``,
``M_s<k>`` at each strand k, each kept compact (``fc_compact``: Geneid and counts, gzipped, the
annotation columns once). Each file's options are read back from its first line and must be
``count.FcOptions.args()``'s for its name (I5); the matrices are then joined in sample order. The
variants are taken at each sample's strand as ``count.called_strands`` calls it, the one
definition, so the cluster holds no copy of the strand rule.
"""

import csv
import re
import shutil
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

import polars as pl

from ctx import tools
from ctx.process import Aligner, ProcessError, cache, count, fc_compact
from ctx.process.align import Alignment, Trimmed, bam_name
from ctx.process.reference import sha256
from ctx.qc import Fastp
from ctx.qc.reads import FILE_COLUMNS

TASKS, VERSIONS, REFS, RUNS = "tasks.tsv", "versions.txt", "refs.txt", "runs"
# Block E's outputs: the featured sample's settings, and the references ctx prepared there.
FEATURED, REFERENCE, INPUTS = "featured", "reference", "inputs.txt"
TASK_COLUMNS = ("sample", "mates", "runs", "r1", "r2")
STAR_KEPT = ("Log.final.out", "ReadsPerGene.out.tab", "flagstat", "idxstats")
FASTP, STAR, COUNTS, TIME = "fastp", "star", "featurecounts", "time"
# versions.txt as block B1 writes it: STAR's line is prefixed with its name; the others print
# their own, which ctx.tools' patterns read.
VERSION_PATTERNS = {"STAR": r"^STAR (\S+)$",
                    **{t: tools.TOOLS[t].pattern for t in ("fastp", "fastqc", "seqkit", "samtools",
                                                           "featureCounts")}}  # fmt: skip
QUOTED = re.compile(r'"([^"]*)"')
# featureCounts options that name a path or a thread count, not a counting setting.
PATH_OPTIONS = ("-T", "-a", "-o")
# GNU time -v's two lines the page reads.
WALL = re.compile(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)")
RSS = re.compile(r"Maximum resident set size \(kbytes\): (\d+)")
SHOWN = 5


@dataclass(frozen=True)
class Source:
    """A copied cluster alignment: its folder, its task list and what the cluster recorded.

    Attributes:
        folder: ``data/raw/<slug>_align/`` as block D writes it.
        tasks: tasks.tsv, one row per sample.
        versions: Each tool's version as versions.txt states it.
        index: The STAR index folder refs.txt names.
        gtf: The annotation file featureCounts read, as refs.txt names it.
        genome: The genome FASTA ctx prepared on the cluster, as block E records it.
        cut_index: The STAR index of the half reference block E built.
    """

    folder: Path
    tasks: pl.DataFrame
    versions: dict[str, str]
    index: str
    gtf: str
    genome: str = ""
    cut_index: str = ""

    def sample_dir(self, sample: str) -> Path:
        """One sample's folder of tool outputs."""
        return self.folder / RUNS / sample

    def featured_dir(self, sample: str, setting: str) -> Path:
        """The featured sample's folder for one setting: ``featured/<sample>/<setting>/``."""
        return self.folder / FEATURED / sample / setting

    def reference(self, name: str) -> Path:
        """The folder of the reference ctx prepared on the cluster: ``main`` or ``cut``."""
        return self.folder / REFERENCE / name


def read(folder: Path) -> tuple[Source | None, list[str]]:
    """The copied folder's task list, versions and references, or every reason it cannot be read.

    A version other than ctx's pin is a problem, as ``ctx.tools`` refuses one here.
    """
    if not folder.is_dir():
        return None, [f"--align-dir {folder} is not a folder"]
    missing = [n for n in (TASKS, VERSIONS, REFS) if not (folder / n).is_file()]
    if missing:
        return None, [f"--align-dir {folder} holds no {missing}; block D copies them"]
    tasks = pl.read_csv(folder / TASKS, separator="\t", infer_schema=False)
    if absent := [c for c in TASK_COLUMNS if c not in tasks.columns]:
        return None, [f"{folder / TASKS} has no column {absent}"]
    text = (folder / VERSIONS).read_text()
    versions, problems = {}, []
    for tool, pattern in VERSION_PATTERNS.items():
        found = re.search(pattern, text, re.M)
        versions[tool] = found[1] if found else ""
        if versions[tool] != tools.TOOLS[tool].version:
            problems.append(f"{folder / VERSIONS}: {tool} {versions[tool] or 'absent'}, not "
                            f"the pin {tools.TOOLS[tool].version}")  # fmt: skip
    refs = _pairs(folder / REFS)
    if missing := [k for k in ("IDX", "GTF") if not refs.get(k)]:
        problems.append(f"{folder / REFS} names no {missing}")
    made = _pairs(folder / REFERENCE / INPUTS)
    source = Source(folder, tasks, versions, *(Path(v).name for v in (
        refs.get("IDX", ""), refs.get("GTF", ""), made.get("GENOME", ""),
        made.get("CUT_INDEX", ""))))  # fmt: skip
    return source, problems


def _pairs(path: Path) -> dict[str, str]:
    """A file of KEY=value lines, as the cluster blocks write them; empty when absent."""
    if not path.is_file():
        return {}
    return dict(line.partition("=")[::2] for line in path.read_text().splitlines() if "=" in line)


def files_of(sheet: pl.DataFrame, source: Source) -> tuple[pl.DataFrame, list[str]]:
    """Each sheet run's FASTQ names from tasks.tsv, as ``reads.fastq_files`` lists files.

    The sheet decides which runs make which sample; a run the cluster grouped otherwise, a run
    it did not align, a run it aligned that the sheet lacks, and a sample with no copied
    folder are each a problem. Paths are the cluster's file names: nothing here opens them.
    """
    given: dict[str, tuple[str, list[str]]] = {}
    for task in source.tasks.iter_rows(named=True):
        runs = task["runs"].split(",")
        mates = [task["r1"].split(",")] + ([task["r2"].split(",")] if task["mates"] == "2" else [])
        for i, run in enumerate(runs):
            given[run] = (task["sample"], [m[i] for m in mates])
    rows, problems = [], []
    for run, sample in sheet.select("run_accession", "sample").rows():
        if run not in given:
            problems.append(f"run {run} of the sheet is in no task of {source.folder / TASKS}")
            continue
        if given[run][0] != sample:
            problems.append(f"run {run}: sample {sample} in the sheet, {given[run][0]} on the "
                            "cluster")  # fmt: skip
        rows += [(run, sample, m, name) for m, name in enumerate(given[run][1], 1)]
    if extra := sorted(set(given) - set(sheet["run_accession"])):
        problems.append(f"{len(extra)} aligned runs are not in the sheet: {extra[:SHOWN]}")
    files = pl.DataFrame(rows, schema=FILE_COLUMNS, orient="row")
    absent = [s for s in files["sample"].unique(maintain_order=True)
              if not source.sample_dir(s).is_dir()]  # fmt: skip
    if absent:
        problems.append(f"{len(absent)} samples have no {source.folder / RUNS}/<sample>: "
                        f"{absent[:SHOWN]}; their task is unfinished or uncopied")  # fmt: skip
    return files, problems


def contigs(source: Source, samples: Sequence[str]) -> tuple[dict[str, int], list[str]]:
    """The genome's contigs and lengths from samtools idxstats, the same for every sample."""
    seen: dict[tuple[tuple[str, int], ...], list[str]] = {}
    for sample in samples:
        path = source.sample_dir(sample) / STAR / f"{sample}.idxstats"
        if not path.is_file():
            return {}, [f"no {path}"]
        rows = [line.split("\t")[:2] for line in path.read_text().splitlines()]
        found = tuple((name, int(length)) for name, length in rows if name != "*")
        seen.setdefault(found, []).append(sample)
    if len(seen) != 1:
        return {}, [f"samtools idxstats names {len(seen)} different genomes over the samples: "
                    f"{[v[:SHOWN] for v in seen.values()]}"]  # fmt: skip
    return dict(next(iter(seen))), []


def _copied(paths: Mapping[str, Path], parent: Path, prefix: str,
            inputs: Mapping[str, object]) -> tuple[Path, bool]:  # fmt: skip
    """The files copied under their names into one keyed folder, keyed by every file's SHA-256."""
    keyed = {**inputs, "source": "cluster", "files": {n: sha256(p) for n, p in paths.items()}}

    def build(dest: Path) -> None:
        """Copy each file in."""
        for name, path in paths.items():
            (dest / name).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, dest / name)

    return cache.keyed(parent, prefix, keyed, build)


def fastp_problems(json_text: str, run: str, fastp: Fastp) -> list[str]:
    """Whether a run's fastp report was made with ``fastp``'s filter and titled by the run.

    fastp records its command line in the JSON; ctx's M2 and M3 run one filter (M3.9).
    """
    command = re.search(r'"command":\s*"([^"]*)"', json_text)
    words = command[1].split() if command else []
    if not (_inside(words, ["-R", run]) and _inside(words, fastp.args())):
        return [f"{run}: fastp ran as {command[1] if command else 'no command'!r}, not with -R "
                f"{run} {' '.join(fastp.args())}"]  # fmt: skip
    return []


def _inside(words: Sequence[str], pair: Sequence[str]) -> bool:
    """Whether an option and its value stand side by side among the words."""
    return any(list(words[i : i + len(pair)]) == list(pair) for i in range(len(words)))


def trimmed(source: Source, files: pl.DataFrame, parent: Path,
            fastp: Fastp) -> tuple[Trimmed | None, list[str]]:  # fmt: skip
    """Each run's fastp report in one keyed folder, laid out as ``align.trim`` writes it."""
    pairs = files.unique("run_accession", keep="first", maintain_order=True)
    reports, problems = {}, []
    for run, sample in pairs.select("run_accession", "sample").rows():
        path = source.sample_dir(sample) / FASTP / f"{run}.fastp.json"
        if not path.is_file():
            problems.append(f"no fastp report {path}")
            continue
        problems += fastp_problems(path.read_text(), run, fastp)
        reports[f"{FASTP}/{run}.fastp.json"] = path
    if problems:
        return None, problems
    folder, built = _copied(reports, parent, "trim", {"fastp_args": fastp.args()})
    given = {r: tuple(Path(p) for p in g["path"]) for (r,), g in
             files.sort("mate").group_by("run_accession", maintain_order=True)}  # fmt: skip
    return Trimmed(folder, built, given), []


def alignment(dirs: Mapping[str, Path], parent: Path, paired: bool,
              index: str) -> tuple[Alignment | None, list[str]]:  # fmt: skip
    """STAR's two tables and samtools' flagstat and idxstats per sample, in one keyed folder.

    Args:
        dirs: Each sample's folder holding ``star/``, in sample order.
        parent: Where keyed folders sit.
        paired: Whether the runs are pairs.
        index: The STAR index folder's name, keyed with the files.

    The BAMs stayed on the cluster: M3 reads a BAM only through featureCounts and samtools,
    whose outputs are here.
    """
    paths = {f"{s}.{n}": d / STAR / f"{s}.{n}" for s, d in dirs.items() for n in STAR_KEPT}
    if missing := [str(p) for p in paths.values() if not p.is_file()]:
        return None, [f"{len(missing)} STAR or samtools files absent: {missing[:SHOWN]}"]
    folder, built = _copied(paths, parent, "align_star", {"index": index})
    return Alignment(folder, built, Aligner.STAR, tuple(dirs), paired), []


def file_name(setting: str, strand: int) -> str:
    """The per-sample file of a setting at a strand: ``s<k>``, else ``<setting>_s<k>``."""
    return setting if setting in count.PLAIN else f"{setting}_s{strand}"


def summary_of(compact: Path) -> Path:
    """The summary featureCounts wrote beside a matrix, kept as written."""
    return compact.with_name(compact.name.removesuffix(fc_compact.SUFFIX) + fc_compact.MATRIX
                             + ".summary")  # fmt: skip


def options_problems(path: Path, sample: str, wanted: count.FcOptions, gtf: str) -> list[str]:
    """I5 on one per-sample file: its options, its BAM and its annotation, from its first line."""
    command, _ = fc_compact.read(str(path))
    words = QUOTED.findall(command)[1:]
    named = {words[i]: words[i + 1] for i in range(len(words) - 1) if words[i] in PATH_OPTIONS}
    rest = [w for i, w in enumerate(words) if w not in PATH_OPTIONS
            and (i == 0 or words[i - 1] not in PATH_OPTIONS)]  # fmt: skip
    problems = []
    if rest != [*wanted.args(), bam_name(sample)]:
        problems.append(f"I5: {path} ran with {rest}, not {[*wanted.args(), bam_name(sample)]}")
    if Path(named.get("-a", "")).name != gtf:
        problems.append(f"I5: {path} counted with {named.get('-a')}, not {gtf}")
    return problems


def join(counts: "Counter", settings: Mapping[str, count.FcOptions], dest: Path) -> list[str]:
    """Each setting's per-sample files joined in sample order into ``dest``, as one matrix.

    Returns:
        I4 and I5 problems; the joined files are written only for settings without one.
    """
    problems = []
    for name, options in settings.items():
        parts = []
        for sample, strand in zip(counts.alignment.samples, options.strands, strict=True):
            wanted = count.FcOptions((strand,), options.paired, options.overlap, options.fraction,
                                     options.multimap)  # fmt: skip
            path = counts.dirs[sample] / COUNTS / f"{file_name(name, strand)}{fc_compact.SUFFIX}"
            problems += options_problems(path, sample, wanted, counts.gtf)
            parts.append(path)
        if not problems:
            out = dest / f"{name}.featureCounts.txt"
            problems += _write(parts, options, counts.annotation, out)
    return problems


def _write(parts: Sequence[Path], options: count.FcOptions, annotation: Path,
           out: Path) -> list[str]:  # fmt: skip
    """One matrix and summary from per-sample files, the annotation columns from their one copy.

    The matrix is featureCounts' own layout, so ``count``'s readers and I5 read it unchanged.
    """
    meta = pl.DataFrame(fc_compact.annotation_rows(str(annotation)), orient="row")
    meta = meta.rename(dict(zip(meta.columns, meta.row(0), strict=True))).slice(1)
    columns, problems = [], []
    for path in parts:
        _, rows = fc_compact.read(str(path))
        genes = [g for g, _ in rows[1:]]
        if genes != meta["Geneid"].to_list():
            problems.append(f"I4: {path} holds other genes, or in another order, than {annotation}")
        columns.append(pl.Series(rows[0][1], [c for _, c in rows[1:]]))
    if problems:
        return problems
    strands = ",".join(map(str, options.strands))
    with out.open("w") as handle:
        handle.write(f"# ctx process: {len(parts)} per-sample featureCounts files joined in sample "
                     f'order, each file\'s options checked (I5); "-s" "{strands}"\n')  # fmt: skip
        meta.hstack(columns).write_csv(handle, separator="\t", quote_style="never")
    summaries = [list(csv.reader(summary_of(p).open(), delimiter="\t")) for p in parts]
    statuses = [[row[0] for row in s] for s in summaries]
    if any(s != statuses[0] for s in statuses):
        return [f"the summaries of {out.name} list different statuses"]
    with out.with_name(out.name + ".summary").open("w") as handle:
        for i, status in enumerate(statuses[0]):
            handle.write("\t".join([status, *(s[i][1] for s in summaries)]) + "\n")
    return []


@dataclass(frozen=True)
class Counter:
    """Counts from the cluster's per-sample featureCounts files, as ``count.count`` returns them.

    Attributes:
        dirs: Each sample's folder holding ``featurecounts/``, in sample order.
        alignment: The alignment those files counted.
        parent: Where keyed folders sit.
        gtf: The annotation file's name featureCounts must have read.
        annotation: The one copy of that annotation's columns (``fc_compact.ANNOTATION``).
    """

    dirs: Mapping[str, Path]
    alignment: Alignment
    parent: Path
    gtf: str
    annotation: Path

    def __call__(self, settings: Mapping[str, count.FcOptions]) -> count.Counts:
        """The settings joined into one keyed folder, keyed by every per-sample file read.

        Raises:
            ProcessError: Listing every per-sample file that is absent or fails I4 or I5.
        """
        files = {f"{s}/{file_name(n, k)}": self.dirs[s] / COUNTS
                 / f"{file_name(n, k)}{fc_compact.SUFFIX}" for n, o in settings.items()
                 for s, k in zip(self.alignment.samples, o.strands, strict=True)}  # fmt: skip
        absent = [str(p) for p in (*files.values(), *map(summary_of, files.values()),
                                   self.annotation) if not p.is_file()]  # fmt: skip
        if absent:
            raise ProcessError(f"{len(absent)} featureCounts files absent: {absent[:SHOWN]}")
        inputs = {"source": "cluster", "settings": {n: asdict(o) for n, o in settings.items()},
                  "samples": list(self.alignment.samples), "annotation": sha256(self.annotation),
                  "files": {n: [sha256(p), sha256(summary_of(p))]
                            for n, p in files.items()}}  # fmt: skip
        found: list[str] = []

        def build(dest: Path) -> None:
            """Join each setting; a failed check leaves the .part folder unfinished."""
            found.extend(join(self, settings, dest))
            if found:
                raise ProcessError("checks failed (D-038):\n  " + "\n  ".join(found))

        folder, built = cache.keyed(self.parent, "count", inputs, build)
        return count.Counts(folder, built, dict(settings))


def wall_seconds(value: str) -> float:
    """GNU time's elapsed wall clock, ``h:mm:ss`` or ``m:ss.ss``, in seconds."""
    seconds = 0.0
    for part in value.split(":"):
        seconds = seconds * 60 + float(part)
    return seconds


def steps(dirs: Sequence[Path]) -> pl.DataFrame:
    """Each tool step's wall time and peak memory over the folders' time/ files (GNU time -v).

    A step is the file's name before the first underscore or dot: fastp, fastqc, seqkit, star,
    samtools or fc (featureCounts). One row per step: files read, summed and largest wall
    time, and the largest maximum resident set size.
    """
    rows = []
    for folder in dirs:
        for path in sorted((folder / TIME).glob("*.time")):
            text = path.read_text()
            wall, rss = WALL.search(text), RSS.search(text)
            if wall and rss:
                step = re.split(r"[_.]", path.name, maxsplit=1)[0]
                rows.append((step, wall_seconds(wall[1]), int(rss[1]) * 1024))
    frame = pl.DataFrame(rows, schema={"step": pl.String, "wall_s": pl.Float64,
                                       "max_rss_bytes": pl.Int64}, orient="row")  # fmt: skip
    return frame.group_by("step", maintain_order=True).agg(
        files=pl.len(), wall_s_sum=pl.col("wall_s").sum(), wall_s_max=pl.col("wall_s").max(),
        max_rss_bytes=pl.col("max_rss_bytes").max())  # fmt: skip
