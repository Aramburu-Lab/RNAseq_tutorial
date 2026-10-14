"""`ctx reads`: a student's own FASTQ files laid out as fetchngs.sbatch lays out a download.

Public reads enter through nf-core/fetchngs; a student's own files enter here (D-006 amended).
The run folder gets the two outdirs fetchngs.sbatch writes, ``nfcore_fetchngs`` (with
``fastq/``) and ``nfcore_fetchngs_meta`` (without), each with its ``results.complete``, so
rnaseq.sbatch takes the samplesheet unchanged and ``ctx gather`` reads either outdir through
``ctx.gather.fetchngs``, the same code and checks, the samplesheet and metadata headed by
fetchngs 1.12.0's columns. The FASTQ files are never moved or copied: ``fastq/`` holds a
symbolic link per file, and the samplesheet names each file where it is.

Standard library only: on Kebnekaise this file runs on the fetch environment's Python, where
ctx is not installed. Usage, on Kebnekaise after kebnekaise_env.sh:
    python3 "$CTX_REPO/src/ctx/gather/own_fastq.py" --run-folder "$RUN" --fastq-dir <folder>
        --table <table> --sample-column sample
"""

import argparse
import csv
import gzip
import hashlib
import io
import itertools
import json
import re
import sys
import tempfile
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

# The layout's one version; ctx.gather.fetchngs reads this version only.
LAYOUT = "1"
VERSIONS_FILE = "ctx_reads_versions.yml"
PARAMS_FILE = "params_ctx_reads.json"
IDS = "ids.csv"
# The outdirs fetchngs.sbatch writes in a run folder, the download's and the metadata run's.
DOWNLOAD, METADATA = "nfcore_fetchngs", "nfcore_fetchngs_meta"
MARKER = "results.complete"
# fetchngs 1.12.0's headers, as tests/fixtures/fetchngs holds them (tests/test_own_fastq.py).
SHEET_COLUMNS = (
    "sample", "fastq_1", "fastq_2", "strandedness", "run_accession", "experiment_accession",
    "sample_accession", "secondary_sample_accession", "study_accession",
    "secondary_study_accession", "submission_accession", "run_alias", "experiment_alias",
    "sample_alias", "study_alias", "library_layout", "library_selection", "library_source",
    "library_strategy", "library_name", "instrument_model", "instrument_platform", "base_count",
    "read_count", "tax_id", "scientific_name", "sample_title", "experiment_title", "study_title",
    "sample_description", "fastq_md5", "fastq_bytes", "fastq_ftp", "fastq_galaxy", "fastq_aspera",
)  # fmt: skip
METADATA_COLUMNS = (
    "id", *SHEET_COLUMNS[4:], "fastq_1", "fastq_2", "md5_1", "md5_2", "single_end",
)  # fmt: skip
MAPPING_COLUMNS = (
    "sample", "experiment_accession", "run_accession", "sample_accession", "experiment_alias",
    "run_alias", "sample_alias", "experiment_title", "sample_title", "sample_description",
)  # fmt: skip
# Measured from the files or named by them; a table column of one of these names is refused.
COMPUTED = frozenset({
    "id", "sample", "fastq_1", "fastq_2", "strandedness", "run_accession",
    "experiment_accession", "library_layout", "base_count", "read_count", "fastq_md5",
    "fastq_bytes", "fastq_ftp", "fastq_galaxy", "fastq_aspera", "md5_1", "md5_2", "single_end",
})  # fmt: skip
# nf-core/rnaseq takes gzipped FASTQ under these endings alone.
SUFFIXES = (".fastq.gz", ".fq.gz")
UNCOMPRESSED = (".fastq", ".fq")
TABLE_SUFFIXES = (".csv", ".tsv", ".txt")
DEFAULT_COLUMN = "sample"
# The mate is the last _1, _2, _R1 or _R2 of the name, before Illumina's _001 chunk if any.
MATE = re.compile(r"^(?P<head>.+?)_R?(?P<mate>[12])(?P<tail>_\d{3})?$")
# Lines per FASTQ record, and files per paired run.
RECORD, PAIR = 4, 2


class ReadsError(ValueError):
    """The files and table cannot be laid out; ``problems`` lists every one found."""

    def __init__(self, problems: list[str]) -> None:
        """Keep the problems and list them in the message."""
        self.problems = problems
        super().__init__(f"{len(problems)} problems:\n  " + "\n  ".join(problems))


@dataclass(frozen=True)
class Run:
    """One run: a single-end file or a pair of mates, of one sample.

    Attributes:
        sample: The table's sample the file names start with.
        key: The file name without its ending and mate token; the run's accession here.
        files: The FASTQ files, mate 1 first, resolved.
    """

    sample: str
    key: str
    files: tuple[Path, ...]


def table_rows(path: Path) -> tuple[list[list[str]], str]:
    """A table as a spreadsheet saves it: rows of cells, and the SHA-256 of its bytes.

    The byte-order mark and Windows line ends are removed; a ``.csv`` is read with commas,
    quotes honoured, every other ending with tabs. ``ctx.gather.curation.read_table`` reads
    an own table through this function, so the two readers cannot differ.

    Args:
        path: The table.

    Returns:
        The rows, the header first, and the SHA-256 of the file's bytes as saved.
    """
    data = path.read_bytes()
    body = data.decode("utf-8").removeprefix(chr(0xFEFF))
    # newline="" keeps a CR inside a quoted cell for csv, which strips the rest itself.
    delimiter = "," if path.suffix == ".csv" else "\t"
    rows = list(csv.reader(io.StringIO(body, newline=""), delimiter=delimiter))
    return rows, hashlib.sha256(data).hexdigest()


def _samples(path: Path, column: str) -> tuple[dict[str, dict[str, str]], list[str]]:
    """Each sample's row of the table, by the sample column, and every problem with them."""
    if path.suffix not in TABLE_SUFFIXES:
        return {}, [f"{path}: a table ends in {', '.join(TABLE_SUFFIXES)}; save it as CSV"]
    if not path.is_file():
        return {}, [f"{path}: no such table"]
    rows, _ = table_rows(path)
    header, body = rows[0] if rows else [], [r for r in rows[1:] if any(c.strip() for c in r)]
    if header.count(column) != 1:
        return {}, [f"{path}: column {column!r} appears {header.count(column)} times, not once"]
    problems = [f"{path}: column {c!r} is measured from the files; remove it from the table"
                for c in header if c in COMPUTED and c != column]  # fmt: skip
    samples: dict[str, dict[str, str]] = {}
    for n, row in enumerate(body, start=2):
        cells = dict(zip(header, (c.strip() for c in row), strict=False))
        name = cells.get(column, "")
        if not name or re.search(r"\s", name):
            problems.append(f"{path} line {n}: sample {name!r} is empty or holds a space")
        elif name in samples:
            problems.append(f"{path} line {n}: sample {name} is in an earlier row too")
        else:
            samples[name] = cells
    return samples, problems


def _fastq_files(folder: Path) -> tuple[list[Path], list[str]]:
    """Every gzipped FASTQ under the folder, and a problem per unusable or repeated name.

    Names starting with a dot (macOS's ``._`` files on external drives) are not FASTQ files.
    """
    if not folder.is_dir():
        return [], [f"{folder}: no such folder"]
    found: dict[str, Path] = {}
    problems = []
    for path in sorted(folder.rglob("*")):
        name = path.name
        if name.startswith(".") or not path.is_file():
            continue
        if name.endswith(UNCOMPRESSED):
            problems.append(f"{path}: not gzipped; nf-core/rnaseq takes .fastq.gz or .fq.gz")
        elif name.endswith(SUFFIXES):
            if re.search(r"\s", str(path.resolve())):
                problems.append(f"{path}: a space in its path, which rnaseq.sbatch refuses")
            if name in found:
                problems.append(f"{path}: the same name as {found[name]}; one name, one file")
            found[name] = path
    if not found and not problems:
        problems.append(f"{folder}: no .fastq.gz or .fq.gz file")
    return list(found.values()), problems


def _stem(name: str) -> str:
    """The name without its FASTQ ending."""
    return next(name[: -len(s)] for s in SUFFIXES if name.endswith(s))


def _owner(key: str, samples: dict[str, dict[str, str]]) -> list[str]:
    """The samples a run key belongs to: the sample's name, alone or followed by ``_``."""
    return [s for s in samples if key == s or key.startswith(f"{s}_")]


def _mates(files: list[Path]) -> tuple[dict[str, dict[str, Path]], list[str]]:
    """Files by run key and mate (``0`` single-end, ``1``, ``2``), and every pairing problem."""
    runs: dict[str, dict[str, Path]] = {}
    problems = []
    for path in files:
        stem = _stem(path.name)
        found = MATE.match(stem)
        key, mate = (found["head"] + (found["tail"] or ""), found["mate"]) if found else (stem, "0")
        if mate in runs.setdefault(key, {}):
            problems.append(f"{path.name} and {runs[key][mate].name}: both are run {key} "
                            f"{'single-end' if mate == '0' else f'mate {mate}'}")  # fmt: skip
        runs[key][mate] = path
    for key, mates in runs.items():
        names = ", ".join(sorted(p.name for p in mates.values()))
        if "0" in mates and len(mates) > 1:
            problems.append(f"{names}: run {key} is both single-end and paired")
        elif "0" not in mates and set(mates) != {"1", "2"}:
            have = next(iter(mates))
            other = "2" if have == "1" else "1"
            problems.append(f"{names}: mate {have} of run {key} has no mate {other}; "
                            "a single-end file has no _1 or _R1 in its name")  # fmt: skip
    return runs, problems


def plan(fastq_dir: Path, table: Path, column: str) -> tuple[list[Run], dict[str, dict[str, str]]]:
    """Every run of every sample, from file names alone, before any file is read.

    Args:
        fastq_dir: The folder holding the FASTQ files, searched with its subfolders.
        table: The student's sample table, CSV or TSV.
        column: The table's column naming each sample; each file name starts with one.

    Returns:
        The runs, sorted by sample and run, and each sample's table row.

    Raises:
        ReadsError: With every problem of the table, the names and the pairing.
    """
    samples, problems = _samples(table, column)
    files, found = _fastq_files(fastq_dir)
    problems += found
    runs, paired = _mates(files)
    problems += paired
    planned = []
    for key, mates in sorted(runs.items()):
        names = ", ".join(sorted(p.name for p in mates.values()))
        owners = _owner(key, samples)
        if len(owners) != 1:
            which = " and ".join(owners) if owners else "no sample"
            problems.append(f"{names}: names {which} of {table.name}; a file starts with "
                            "exactly one sample's name, then _ or its ending")  # fmt: skip
            continue
        ordered = tuple(mates[m].resolve() for m in sorted(mates))
        planned.append(Run(owners[0], key, ordered))
    layouts: dict[str, set[int]] = {}
    for run in planned:
        layouts.setdefault(run.sample, set()).add(len(run.files))
    problems += [f"sample {s}: no FASTQ file starts with its name" for s in samples
                 if s not in layouts]  # fmt: skip
    problems += [f"sample {s}: single-end and paired runs together, which nf-core/rnaseq "
                 "refuses" for s, n in layouts.items() if len(n) > 1]  # fmt: skip
    if problems:
        raise ReadsError(problems)
    return sorted(planned, key=lambda r: (r.sample, r.key)), samples


def measure(path: Path) -> tuple[int, int, str]:
    """Reads and bases of one gzipped FASTQ, and the MD5 of its bytes as stored.

    Raises:
        ReadsError: If the file is not gzip, or holds a record that is not whole.
    """
    with path.open("rb") as handle:
        md5 = hashlib.file_digest(handle, "md5").hexdigest()
    reads = bases = 0
    try:
        with gzip.open(path, "rb") as handle:
            lines = (line.rstrip(b"\r\n") for line in handle)
            for reads, record in enumerate(itertools.batched(lines, RECORD), start=1):  # noqa: B007
                if not _whole(record):
                    raise ReadsError([f"{path.name} record {reads}: not a whole FASTQ record"])
                bases += len(record[1])
    except (OSError, EOFError) as error:
        raise ReadsError([f"{path.name}: {error}"]) from error
    if reads == 0:
        raise ReadsError([f"{path.name}: no FASTQ record"])
    return reads, bases, md5


def _whole(record: tuple[bytes, ...]) -> bool:
    """Four lines: a header, a sequence, a separator, and one quality per base."""
    return (len(record) == RECORD and record[0].startswith(b"@")
            and record[2].startswith(b"+") and len(record[3]) == len(record[1]))  # fmt: skip


def _row(run: Run, cells: dict[str, str], measured: list[tuple[int, int, str]]) -> dict[str, str]:
    """The run's metadata row, under ``METADATA_COLUMNS``."""
    paired = len(run.files) == PAIR
    files = [str(p) for p in run.files] + [""] * (PAIR - len(run.files))
    md5s = [m[2] for m in measured] + [""] * (PAIR - len(measured))
    row = {c: cells.get(c, "") for c in METADATA_COLUMNS if c not in COMPUTED}
    return row | {
        "id": run.key, "run_accession": run.key, "experiment_accession": run.sample,
        "library_layout": "PAIRED" if paired else "SINGLE",
        "base_count": str(sum(m[1] for m in measured)), "read_count": str(measured[0][0]),
        "fastq_md5": ";".join(m[2] for m in measured),
        "fastq_bytes": ";".join(str(p.stat().st_size) for p in run.files),
        "fastq_ftp": ";".join(p.name for p in run.files), "fastq_galaxy": "", "fastq_aspera": "",
        "fastq_1": files[0], "fastq_2": files[1], "md5_1": md5s[0], "md5_2": md5s[1],
        "single_end": "false" if paired else "true",
    }  # fmt: skip


def rows(runs: list[Run], samples: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    """One metadata row per run, every file read once; mates must hold equal reads.

    Raises:
        ReadsError: With every unreadable file and every pair of unequal mates.
    """
    out, problems = [], []
    for run in runs:
        measured = []
        for path in run.files:
            try:
                measured.append(measure(path))
            except ReadsError as error:
                problems += error.problems
        if len(measured) < len(run.files):
            continue
        if len({m[0] for m in measured}) > 1:
            counts = " and ".join(
                f"{p.name} {m[0]:,}" for p, m in zip(run.files, measured, strict=True)
            )
            problems.append(f"run {run.key}: mates hold {counts} reads; a pair holds equal reads")
            continue
        out.append(_row(run, samples[run.sample], measured))
    if problems:
        raise ReadsError(problems)
    return out


def _csv(path: Path, columns: tuple[str, ...], data: list[dict[str, str]], quoted: bool) -> None:
    """A table as fetchngs writes it: samplesheet and mappings quoted, metadata tab-separated."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, lineterminator="\n", extrasaction="ignore",
            quoting=csv.QUOTE_ALL if quoted else csv.QUOTE_MINIMAL,
            delimiter="," if quoted else "\t",
        )  # fmt: skip
        writer.writeheader()
        writer.writerows(data)


def _outdir(
    launch: Path, data: list[dict[str, str]], download: bool, given: dict[str, str]
) -> None:
    """One outdir, ``<launch>/results``, and the launch folder's input list and marker."""
    out = launch / "results"
    ids = "".join(f"{s}\n" for s in dict.fromkeys(r["experiment_accession"] for r in data))
    for folder in (launch, out):
        folder.mkdir(parents=True, exist_ok=True)
        (folder / IDS).write_text(ids)
    sheet = [r | {"sample": r["experiment_accession"], "strandedness": "auto"} for r in data]
    _csv(out / "samplesheet" / "samplesheet.csv", SHEET_COLUMNS, sheet, quoted=True)
    _csv(out / "samplesheet" / "id_mappings.csv", MAPPING_COLUMNS, sheet, quoted=True)
    for sample in dict.fromkeys(r["experiment_accession"] for r in data):
        mine = [r for r in data if r["experiment_accession"] == sample]
        _csv(out / "metadata" / f"{sample}.runinfo_ftp.tsv", METADATA_COLUMNS, mine, quoted=False)
    info = out / "pipeline_info"
    info.mkdir()
    (info / VERSIONS_FILE).write_text(
        f"ctx reads: layout {LAYOUT}\nPython: {sys.version.split()[0]}\n"
    )
    params = {"input": IDS, "nf_core_pipeline": "rnaseq", "skip_fastq_download": not download}
    (info / PARAMS_FILE).write_text(json.dumps(params | given, indent=4) + "\n")
    if download:
        (out / "fastq" / "md5").mkdir(parents=True)
        for r in data:
            for key in ("1", "2"):
                if source := r[f"fastq_{key}"]:
                    name = Path(source).name
                    (out / "fastq" / name).symlink_to(source)
                    (out / "fastq" / "md5" / f"{name}.md5").write_text(
                        f"{r[f'md5_{key}']}  {name}\n"
                    )
    (launch / MARKER).write_text(datetime.now(UTC).isoformat(timespec="seconds") + "\n")


def write(run_folder: Path, fastq_dir: Path, table: Path, column: str) -> list[Path]:
    """Lay out the student's FASTQ files in the run folder, both outdirs or neither.

    The names and the table are checked first, every problem listed together; then every
    file is read once, for its MD5, reads and bases. Both outdirs are built in a fresh hidden
    folder of the run folder and moved into place only when both are whole.

    Args:
        run_folder: The run folder, as fetchngs.sbatch and rnaseq.sbatch take it.
        fastq_dir: The folder holding the FASTQ files, searched with its subfolders.
        table: The student's sample table, CSV or TSV, one row per sample.
        column: The table's column naming each sample.

    Returns:
        The two outdirs: the download's, then the metadata run's.

    Raises:
        ReadsError: With every problem found; the run folder is then left as it was.
    """
    run_folder, fastq_dir, table = (
        p.expanduser().absolute() for p in (run_folder, fastq_dir, table)
    )
    launches = [run_folder / name for name in (DOWNLOAD, METADATA)]
    refused = [f"{p} exists; move it aside, then run this again" for p in launches if p.exists()]
    if re.search(r"\s", str(run_folder)):
        refused.append(f"{run_folder}: a space in the run folder, which rnaseq.sbatch refuses")
    try:
        planned, samples = plan(fastq_dir, table, column)
    except ReadsError as error:
        raise ReadsError(refused + error.problems) from None
    if refused:
        raise ReadsError(refused)
    data = rows(planned, samples)
    given = {"fastq_dir": str(fastq_dir), "table": str(table), "sample_column": column,
             f"sha256:{table.name}": table_rows(table)[1]}  # fmt: skip
    run_folder.mkdir(parents=True, exist_ok=True)
    # A fresh folder each time, so an attempt cut short never blocks the next one.
    building = Path(tempfile.mkdtemp(prefix=".ctx_reads_", dir=run_folder))
    for launch, download in zip(launches, (True, False), strict=True):
        _outdir(building / launch.name, data, download, given)
    for launch in launches:
        (building / launch.name).rename(launch)
    building.rmdir()
    return [launch / "results" for launch in launches]


def summary(outdir: Path) -> str:
    """Samples, runs and reads of a written outdir, as the command prints them last."""
    with (outdir / "samplesheet" / "samplesheet.csv").open(newline="") as handle:
        sheet = list(csv.DictReader(handle))
    reads = sum(int(r["read_count"]) for r in sheet)
    paired = sum(r["library_layout"] == "PAIRED" for r in sheet)
    return (f"{len({r['sample'] for r in sheet})} samples; {len(sheet)} runs, {paired} paired; "
            f"{reads:,} reads or pairs")  # fmt: skip


def parser() -> argparse.ArgumentParser:
    """The options ``ctx reads`` takes, under the same names, for the file run as a script."""
    found = argparse.ArgumentParser(prog="own_fastq.py", description=__doc__.splitlines()[0])
    found.add_argument("--run-folder", type=Path, required=True)
    found.add_argument("--fastq-dir", type=Path, required=True)
    found.add_argument("--table", type=Path, required=True)
    found.add_argument("--sample-column", default=DEFAULT_COLUMN)
    return found


def main(argv: list[str] | None = None) -> int:
    """Run as a script: lay out the files, print both outdirs and the summary, or every problem.

    Returns:
        0 on success, 1 when a problem was found.
    """
    args = parser().parse_args(argv)
    try:
        written = write(args.run_folder, args.fastq_dir, args.table, args.sample_column)
    except ReadsError as error:
        print(error, file=sys.stderr)
        return 1
    print(*written, summary(written[0]), sep="\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
