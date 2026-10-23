"""The gathered FASTQ files and fastp, for M3's laptop route (`ctx process --genome`).

The harmonised sheet names each run's FASTQ files relative to the fetchngs outdir, as gather
wrote them after checking their MD5s. This module opens those files and nothing else of the
outdir's layout (D-031). M2's read stage reads nf-core's FastQC and Trim Galore tables instead
(D-116); FastQC, seqkit and fastp's grid no longer run here.
"""

import json
import os
from pathlib import Path

import polars as pl

from ctx import tools
from ctx.qc import Fastp, QcError
from ctx.qc.registry import RAW_COUNTS

# FASTQ statuses (ctx.gather.fetchngs.STATUSES) of a file worth reading: on disk, and not
# known to differ from what ENA published.
READABLE = ("verified", "no_md5")
SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
FILE_COLUMNS = {"run_accession": pl.String, "sample": pl.String, "mate": pl.Int32,
                "path": pl.String}  # fmt: skip


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


def trimmed_names(title: str, mates: int) -> list[str]:
    """File names of one run's reads as ``run_fastp`` writes them, read 1 first."""
    return [f"{title}_{m}.fastq.gz" for m in (1, 2)] if mates > 1 else [f"{title}.fastq.gz"]


def run_fastp(
    paths: list[str], out: Path, title: str, params: Fastp, threads: int,
    trimmed: Path | None = None,
) -> list[Path]:  # fmt: skip
    """One fastp run on a run's files, titled so MultiQC names its sample by the title.

    Args:
        paths: The run's FASTQ files, read 1 first.
        out: Where the JSON report and the log go.
        title: The report's title, the run's name in MultiQC.
        params: The filter.
        threads: fastp's worker threads.
        trimmed: A folder for the reads fastp keeps, ``<title>.fastq.gz`` or ``<title>_1`` and
            ``_2`` for a pair; None (M2) filters and reports without writing reads.

    Returns:
        The trimmed FASTQ files written, read 1 first, absolute; none without ``trimmed``.
    """
    # Absolute: fastp runs in `out`, where a relative path misses, and a missed -j is silent.
    paths, out = [str(Path(p).resolve()) for p in paths], out.resolve()
    trimmed = trimmed.resolve() if trimmed is not None else None
    pair = ["-I", paths[1]] if len(paths) > 1 else []
    args = ["-i", paths[0], *pair, "-R", title, "-j", str(out / f"{title}.fastp.json"),
            "-h", os.devnull, "-w", str(threads), *params.args()]  # fmt: skip
    written = []
    if trimmed is not None:
        trimmed.mkdir(parents=True, exist_ok=True)
        written = [trimmed / name for name in trimmed_names(title, len(paths))]
        for flag, path in zip(("-o", "-O"), written, strict=False):
            args += [flag, str(path)]
    tools.run("fastp", args, cwd=out, log=out / f"{title}.fastp.log")
    return written


def raw_config(path: Path) -> Path:
    """A MultiQC config that keeps every read and base count raw (D-029)."""
    path.write_text(json.dumps(RAW_COUNTS, indent=1) + "\n")
    return path
