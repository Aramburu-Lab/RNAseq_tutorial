"""M2a from a cluster's reports (ED.6): fastp, FastQC and seqkit as Kebnekaise ran them.

The folder is ED.6's block D copy of one project's run folder: ``tasks.tsv``, one row per
sample with its runs and FASTQ names, and per sample ``runs/<sample>/fastp/<run>.fastp.json``,
``fastqc/<stem>_fastqc.zip`` and ``seqkit/<sample>.seqkit_stats.tsv``, with the fastp grid in
``grid/``. The FASTQ files stay on the cluster, so the sheet names none and the task list names
them all. MultiQC reads every report, as it reads ``reads.run_tools``' folders, and nothing in
the folder is written to.
"""

from collections.abc import Mapping, Sequence
from dataclasses import replace
from pathlib import Path

import polars as pl

from ctx.qc import Fastp, QcError, multiqc_io, reads

TASKS = "tasks.tsv"
TASK_COLUMNS = ("sample", "mates", "runs", "r1", "r2")
GRID = "grid"
# fastp's options naming a file, the title or threads; every other word is the read filter.
NAMED = ("-i", "-I", "-R", "-j", "-h", "-w", "-o", "-O")
LAYOUT_MATES = {"SINGLE": 1, "PAIRED": 2}


def tasks(align_dir: Path) -> pl.DataFrame:
    """One row per run of the cluster's task list, in its order: the run's sample and files.

    Raises:
        QcError: If the list is absent, lacks a column, or a row's runs and files disagree.
    """
    path = align_dir / TASKS
    if not path.is_file():
        raise QcError(f"no {TASKS} in --align-dir {align_dir}")
    table = pl.read_csv(path, separator="\t", infer_schema=False)
    if missing := [c for c in TASK_COLUMNS if c not in table.columns]:
        raise QcError(f"{path} has no {missing}")
    rows, bad = [], []
    for task in table.iter_rows(named=True):
        runs, first = task["runs"].split(","), task["r1"].split(",")
        second = task["r2"].split(",") if task["mates"] == "2" else [""] * len(runs)
        if task["mates"] not in {"1", "2"} or not len(runs) == len(first) == len(second):
            bad.append(task["sample"])
            continue
        rows += [(r, task["sample"], int(task["mates"]), f1, f2)
                 for r, f1, f2 in zip(runs, first, second, strict=True)]  # fmt: skip
    if bad:
        raise QcError(f"{path}: samples whose mates, runs and files disagree: {bad}")
    schema = {"run_accession": pl.String, "task": pl.String, "mates": pl.Int32,
              "fastq_1": pl.String, "fastq_2": pl.String}  # fmt: skip
    return pl.DataFrame(rows, schema=schema, orient="row")


def agree(sheet: pl.DataFrame, runs: pl.DataFrame, path: Path) -> None:
    """Stop unless the task list holds the sheet's runs, each once, in the sheet's layout.

    Raises:
        QcError: Naming the runs only in the sheet, only in the list, listed twice, or laid out
            otherwise than the sheet's ``library_layout`` says.
    """
    listed = runs["run_accession"].to_list()
    in_sheet, in_list = set(sheet["run_accession"]), set(listed)
    twice = sorted({r for r in listed if listed.count(r) > 1})
    only_sheet, only_list = sorted(in_sheet - in_list), sorted(in_list - in_sheet)
    if only_sheet or only_list or twice:
        raise QcError(f"{path} and the gather sheet disagree: runs only in the sheet "
                      f"{only_sheet}; only in the task list {only_list}; "
                      f"listed twice {twice}")  # fmt: skip
    if "library_layout" not in sheet.columns:
        return
    layout = sheet.select("run_accession", sheet_mates=pl.col("library_layout")
                          .replace_strict(LAYOUT_MATES, default=None))  # fmt: skip
    joined = runs.join(layout, on="run_accession").filter(pl.col("sheet_mates").is_not_null())
    other = joined.filter(pl.col("mates") != pl.col("sheet_mates"))["run_accession"].to_list()
    if other:
        raise QcError(f"{path}: runs whose FASTQ files are not the sheet's library_layout: "
                      f"{other}")  # fmt: skip


def files(sheet: pl.DataFrame, align_dir: Path) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """``reads.fastq_files``' two tables from the task list, and the list by run.

    Each file's path is its name on the cluster: only its stem is read, to name its reports.
    Every run of the sheet is read, in the sheet's order, so none is left unread.

    Raises:
        QcError: If the task list and the sheet disagree (``agree``).
    """
    runs = tasks(align_dir)
    agree(sheet, runs, align_dir / TASKS)
    ordered = sheet.select("run_accession", "sample").join(runs, on="run_accession",
                                                            maintain_order="left")  # fmt: skip
    rows = [(r["run_accession"], r["sample"], m, r[f"fastq_{m}"])
            for r in ordered.iter_rows(named=True) for m in range(1, r["mates"] + 1)]  # fmt: skip
    found = pl.DataFrame(rows, schema=reads.FILE_COLUMNS, orient="row")
    return found, pl.DataFrame(schema={"run_accession": pl.String, "status": pl.String}), ordered


def reports(align_dir: Path, ordered: pl.DataFrame) -> list[Path]:
    """Every report the reads stage reads, as block C names them.

    Raises:
        QcError: If any is absent; the message lists them all.
    """
    paths = []
    for task, group in ordered.group_by("task", maintain_order=True):
        folder = align_dir / "runs" / str(task[0])
        paths.append(folder / "seqkit" / f"{task[0]}.seqkit_stats.tsv")
        paths += [folder / "fastp" / f"{r}.fastp.json" for r in group["run_accession"]]
        names = [n for c in ("fastq_1", "fastq_2") for n in group[c] if n]
        paths += [folder / "fastqc" / f"{reads.stem(n)}_fastqc.zip" for n in names]
    if absent := [str(p.relative_to(align_dir)) for p in paths if not p.is_file()]:
        raise QcError(f"{len(absent)} reports absent from --align-dir {align_dir}: {absent}")
    return paths


def split_command(command: str) -> tuple[dict[str, str], list[str]]:
    """The named options of a fastp command, and the words left: the read filter."""
    words, named, rest = command.split()[1:], {}, []
    i = 0
    while i < len(words):
        if words[i] in NAMED and i + 1 < len(words):
            named[words[i]] = words[i + 1]
            i += 2
        else:
            rest.append(words[i])
            i += 1
    return named, rest


def check_commands(
    mqc: multiqc_io.MultiqcRun, expected: Mapping[str, tuple[Sequence[str], Sequence[str]]]
) -> None:
    """Stop unless each fastp report ran on its FASTQ files with its filter.

    Args:
        mqc: The MultiQC run that read the reports.
        expected: Title to (FASTQ file names, read 1 first; the filter's options).

    Raises:
        QcError: Naming each title whose report is absent, read other files, or ran another
            filter.
    """
    commands = dict(multiqc_io.table(mqc, "fastp").select("Sample", "command").rows())
    wrong = []
    for title, (names, filter_args) in expected.items():
        if title not in commands:
            wrong.append(f"{title}: no report")
            continue
        named, rest = split_command(commands[title])
        given = [Path(named[f]).name for f in ("-i", "-I") if f in named]
        if named.get("-R") != title or given != list(names) or rest != list(filter_args):
            wrong.append(f"{title}: ran {given} with {rest}, "
                         f"not {list(names)} with {list(filter_args)}")  # fmt: skip
    if wrong:
        raise QcError(f"fastp reports not run as this ctx qc asks: {'; '.join(wrong)}")


def read(ordered: pl.DataFrame, found: pl.DataFrame, paths: Sequence[Path], dest: Path,
         config: Path, fastp: Fastp) -> multiqc_io.MultiqcRun:  # fmt: skip
    """MultiQC on the reports, each run's and each file's found, and fastp's filter checked.

    Raises:
        QcError: If a report MultiQC did not read, or one ran another filter or other files.
    """
    mqc = multiqc_io.run(paths, dest, config=config, modules=reads.MODULES)
    stems = [reads.stem(p) for p in found["path"]]
    for module, units in (("fastqc", stems), ("seqkit", stems)):
        multiqc_io.require_sources(mqc, module, units, {u: u for u in units})
    names = {r: [n for n in (f1, f2) if n] for r, f1, f2 in
             ordered.select("run_accession", "fastq_1", "fastq_2").rows()}  # fmt: skip
    check_commands(mqc, {r: (n, fastp.args()) for r, n in names.items()})
    return mqc


def grid(align_dir: Path, dest: Path, names: Sequence[str], fastp: Fastp) -> pl.DataFrame:
    """``reads.fastp_grid``'s table from the cluster's grid reports on the run named ``names``.

    Raises:
        QcError: If a cell's report is absent, or one ran on other files or another filter.
    """
    cells = reads.grid_cells((reads.GRID_Q, reads.GRID_L))
    paths = [align_dir / GRID / f"{t}.fastp.json" for t in cells]
    if absent := [p.name for p in paths if not p.is_file()]:
        raise QcError(f"{len(absent)} grid reports absent from {align_dir / GRID}: {absent}")
    dest.mkdir(parents=True, exist_ok=True)
    mqc = multiqc_io.run(paths, dest, config=reads.raw_config(dest / "multiqc_config.json"),
                         modules=("fastp",))  # fmt: skip
    filters = {t: (names, replace(fastp, q=q, l=n).args()) for t, (q, n) in cells.items()}
    check_commands(mqc, filters)
    return reads.grid_table(mqc, cells)
