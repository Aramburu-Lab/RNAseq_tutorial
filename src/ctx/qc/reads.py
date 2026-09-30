"""M2a on reads: FastQC, fastp and seqkit on the gathered FASTQs, read back through MultiQC.

The harmonised sheet names each run's FASTQ files relative to the fetchngs outdir, as gather
wrote them after checking their MD5s. This module opens those files and nothing else of the
outdir's layout (D-031). Every number comes back through MultiQC (`ctx.qc.multiqc_io`).
"""

import json
import os
import re
from collections.abc import Sequence
from dataclasses import replace
from pathlib import Path

import polars as pl

from ctx import tools
from ctx.qc import Fastp, QcError, multiqc_io
from ctx.qc.registry import RAW_COUNTS, Registry

# FASTQ statuses (ctx.gather.fetchngs.STATUSES) of a file worth reading: on disk, and not
# known to differ from what ENA published.
READABLE = ("verified", "no_md5")
SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
MODULES = ("fastqc", "fastp", "seqkit")
FILE_COLUMNS = {"run_accession": pl.String, "sample": pl.String, "mate": pl.Int32,
                "path": pl.String}  # fmt: skip
# A run is judged on its own row (rank 0) and its read-1 row (1); read 2 (2) is not judged.
READ_2 = 2
# The -q × -l grid on one run: fastp's defaults (15, 15) inside it, and -l past the longest
# Step-1 read (101 nt), so the map always shows where the length floor bites.
GRID_Q = (10, 15, 20, 25, 30, 35)
GRID_L = (15, 25, 36, 50, 75, 102)
# A cell of MultiQC's FastQC plot data: "(position, value)".
POINT = re.compile(r"\(\s*([-0-9.eE]+)\s*,\s*([-0-9.eE]+)\s*\)")
GRID_COLUMNS = ("fastp-filtering_result_passed_filter_reads", "fastp-pct_surviving_reads")


def stem(path: str) -> str:
    """A FASTQ's name without its FASTQ suffix, as FastQC and seqkit name its sample."""
    name = Path(path).name
    return next((name[: -len(s)] for s in SUFFIXES if name.endswith(s)), name)


def fastq_files(sheet: pl.DataFrame, fetchngs_dir: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    """The FASTQ files of every readable run, and the runs left unread with their status.

    A run is read when each of its files is verified or has no MD5 to check; a mismatched,
    missing or undownloaded file leaves its whole run unread.

    Raises:
        QcError: If a file the sheet calls readable is not on disk (all such files listed), or
            no run is readable.
    """
    rows, skipped, absent = [], [], []
    for row in sheet.iter_rows(named=True):
        mates = [(m, row[f"fastq_{m}"], row[f"status_{m}"]) for m in (1, 2) if row[f"fastq_{m}"]]
        bad = [s for _, _, s in mates if s not in READABLE]
        if not mates or bad:
            skipped.append((row["run_accession"], ", ".join(bad) or "not_downloaded"))
            continue
        for mate, relative, _ in mates:
            # Absolute: each tool runs in its own log folder, where a relative path misses.
            path = (fetchngs_dir / relative).resolve()
            absent += [str(path)] * (not path.is_file())
            rows.append((row["run_accession"], row["sample"], mate, str(path)))
    if absent:
        raise QcError(f"{len(absent)} FASTQ files the sheet calls readable are not on disk: "
                      f"{absent}; is --fetchngs-dir the outdir gather read?")  # fmt: skip
    if not rows:
        raise QcError(f"no run has a readable FASTQ ({READABLE}); unread: {skipped}")
    return (pl.DataFrame(rows, schema=FILE_COLUMNS, orient="row"),
            pl.DataFrame(skipped, schema=["run_accession", "status"], orient="row"))  # fmt: skip


def _fastp(paths: list[str], out: Path, title: str, params: Fastp, threads: int) -> None:
    """One fastp run on a run's files, titled so MultiQC names its sample by the title.

    No -o: fastp then filters and reports without writing the reads.
    """
    pair = ["-I", paths[1]] if len(paths) > 1 else []
    args = ["-i", paths[0], *pair, "-R", title, "-j", str(out / f"{title}.fastp.json"),
            "-h", os.devnull, "-w", str(threads), *params.args()]  # fmt: skip
    tools.run("fastp", args, cwd=out, log=out / f"{title}.fastp.log")


def run_tools(files: pl.DataFrame, logs: Path, *, threads: int, fastp: Fastp) -> list[Path]:
    """FastQC and seqkit on every file, fastp on every run; returns the folders MultiQC reads.

    Args:
        files: ``fastq_files``'s first table.
        logs: Where each tool's folder goes.
        threads: CPUs per tool.
        fastp: The fastp filter.

    Returns:
        The FastQC, fastp and seqkit folders.
    """
    folders = [logs / m for m in MODULES]
    for folder in folders:
        folder.mkdir(parents=True, exist_ok=True)
    paths = files["path"].to_list()
    tools.run("fastqc", ["--quiet", "--threads", str(threads), "--outdir", str(folders[0]),
                         *paths], cwd=folders[0], log=folders[0] / "fastqc.log")  # fmt: skip
    for (run,), mates in files.sort("mate").group_by("run_accession", maintain_order=True):
        _fastp(mates["path"].to_list(), folders[1], str(run), fastp, threads)
    tools.run("seqkit", ["stats", "-a", "-T", "-j", str(threads), "-o",
                         str(folders[2] / "seqkit_stats.tsv"), *paths],
              cwd=folders[2], log=folders[2] / "seqkit.log")  # fmt: skip
    return folders


def raw_config(path: Path) -> Path:
    """A MultiQC config that keeps every read and base count raw (D-029)."""
    path.write_text(json.dumps(RAW_COUNTS, indent=1) + "\n")
    return path


def run_values(mqc: multiqc_io.MultiqcRun, files: pl.DataFrame, registry: Registry) -> pl.DataFrame:
    """One row per run: every read-stage registry metric MultiQC wrote, by the metric's column.

    fastp's row is named by the run; FastQC's and seqkit's by each file. A run takes a column
    from its own row, else from its read-1 file's row; read 2 is not judged.

    Raises:
        QcError: If a MultiQC sample maps to no file, or a metric's column is missing.
    """
    ranks = {str(r): (str(r), 0) for r in files["run_accession"]}
    ranks |= {
        stem(p): (str(r), m) for r, m, p in files.select("run_accession", "mate", "path").rows()
    }
    metrics = registry.stage("reads")
    frames = []
    for name in dict.fromkeys(m.table for m in metrics):
        table = multiqc_io.table(mqc, name)
        wanted = [m for m in metrics if m.table == name]
        if missing := [m.column for m in wanted if m.column not in table.columns]:
            raise QcError(f"MultiQC's {name} table has no {missing}")
        numbers = [m.column for m in wanted if m.kind == "numeric"]
        table = multiqc_io.numeric(table, numbers).select("Sample", *(m.column for m in wanted))
        frames.append(multiqc_io.by_unit(table, ranks, "run_accession", READ_2))
    values = files.unique("run_accession", keep="first", maintain_order=True).select(
        "run_accession", "sample"
    )
    for frame in frames:
        values = values.join(frame, on="run_accession", how="left")
    return values


def quality_curves(mqc: multiqc_io.MultiqcRun, files: pl.DataFrame) -> pl.DataFrame:
    """Mean base quality along each run's read-1 file, from MultiQC's FastQC plot data.

    MultiQC writes one row per FASTQ, a cell per position as "(position, mean)"; its
    positions are FastQC's bins, one base each at first, then wider.
    """
    table = pl.read_csv(mqc.data / "fastqc_per_base_sequence_quality_plot.txt", separator="\t",
                        infer_schema=False)  # fmt: skip
    multiqc_io.map_names(table["Sample"], {stem(p): stem(p) for p in files["path"]})
    first = {stem(p): r for r, m, p in files.select("run_accession", "mate", "path").rows()
             if m == 1}  # fmt: skip
    rows = [(first[row[0]], float(x), float(y))
            for row in table.rows() if row[0] in first
            for cell in row[1:] if cell for x, y in POINT.findall(cell)]  # fmt: skip
    return pl.DataFrame(rows, schema=["run_accession", "position", "quality"], orient="row")


def fastp_grid(
    paths: Sequence[str], logs: Path, dest: Path, *, grid: tuple[Sequence[int], Sequence[int]],
    fastp: Fastp, threads: int,
) -> pl.DataFrame:  # fmt: skip
    """Reads fastp keeps from one run at every -q and -l of a grid, the rest of ``fastp`` fixed.

    Args:
        paths: The run's FASTQ files, read 1 first.
        logs: Where the grid's fastp reports go.
        dest: MultiQC's output folder for the grid.
        grid: The -q values and the -l values.
        fastp: The other fastp settings.
        threads: CPUs for fastp.

    Returns:
        One row per cell: q, l, the reads that passed, and their percentage of the reads in.
    """
    logs.mkdir(parents=True, exist_ok=True)
    cells = {f"grid_q{q}_l{length}": (q, length) for q in grid[0] for length in grid[1]}
    for title, (q, length) in cells.items():
        _fastp(list(paths), logs, title, replace(fastp, q=q, l=length), threads)
    mqc = multiqc_io.run([logs], dest, config=raw_config(logs / "multiqc_config.json"),
                         modules=("fastp",))  # fmt: skip
    table = multiqc_io.numeric(multiqc_io.table(mqc, "general_stats"), GRID_COLUMNS)
    multiqc_io.map_names(table["Sample"], {t: t for t in cells})
    return table.select(
        q=pl.col("Sample").replace_strict({t: q for t, (q, _) in cells.items()}),
        l=pl.col("Sample").replace_strict({t: n for t, (_, n) in cells.items()}),
        passed_reads=pl.col(GRID_COLUMNS[0]).cast(pl.Int64),
        passed_pct=pl.col(GRID_COLUMNS[1]),
    ).sort("q", "l")
