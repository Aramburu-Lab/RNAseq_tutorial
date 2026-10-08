"""What a pipeline run cost: peak memory and real time per process, from Nextflow's trace.

Every launch writes one ``pipeline_info/execution_trace_<time>.txt``, read by its own column
names. A ``CACHED`` row repeats a task an earlier launch ran, with that launch's numbers, so
it is left out: summed over every launch's trace, the rest is each task once. The sample is
the text in parentheses at the end of ``name``. Nothing here writes into the outdir.
"""

import csv
import functools
import re
from collections.abc import Sequence
from pathlib import Path

import polars as pl

from ctx.report import text
from ctx.report.builder import Group, ResultTable

PAGE = "gather"
TRACES = "pipeline_info/execution_trace_*.txt"
# Nextflow's MemoryUnit prints powers of 1024 under these names (MemoryUnit.groovy).
MEMORY = {"B": 1, "KB": 1024, "MB": 1024**2, "GB": 1024**3, "TB": 1024**4, "PB": 1024**5}
# Nextflow's Duration prints days, hours, minutes, seconds and milliseconds, e.g. "1m 3s".
DURATION = {"d": 86400.0, "h": 3600.0, "m": 60.0, "s": 1.0, "ms": 0.001}
CACHED = "CACHED"
_PART = re.compile(r"(\d+(?:\.\d+)?)(ms|d|h|m|s)")
_NAME = re.compile(r"^(?P<process>[^(]+?)(?: \((?P<sample>.*)\))?$")


class TraceError(ValueError):
    """A trace without the columns this reader needs."""


def memory(value: str) -> int | None:
    """Bytes from a trace's memory, such as ``3 MB`` or ``1.1 GB``; None for ``-`` or empty."""
    value = value.strip()
    if value in {"", "-"}:
        return None
    number, _, unit = value.partition(" ")
    return round(float(number) * MEMORY[unit or "B"])


def seconds(value: str) -> float | None:
    """Seconds from a trace's duration, such as ``19ms`` or ``1h 2m 3s``; None for ``-``."""
    value = value.strip()
    if value in {"", "-"}:
        return None
    parts = _PART.findall(value)
    if "".join(n + u for n, u in parts) != value.replace(" ", ""):
        raise TraceError(f"unreadable duration {value!r}")
    return sum(float(n) * DURATION[u] for n, u in parts)


def read(outdir: Path) -> pl.DataFrame:
    """Every row of every trace in an outdir: process, sample, status, peak RSS and real time.

    Raises:
        TraceError: If a trace lacks ``name``, ``status``, ``peak_rss`` or ``realtime``.
    """
    rows = []
    for path in sorted(outdir.glob(TRACES)):
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            needed = {"name", "status", "peak_rss", "realtime"}
            if missing := sorted(needed - set(reader.fieldnames or ())):
                raise TraceError(f"{path}: no {missing} column")
            for row in reader:
                name = _NAME.fullmatch(row["name"].strip())
                assert name is not None  # noqa: S101 - the pattern accepts any name
                process = name["process"].rpartition(":")[2]
                rows.append((process, name["sample"] or "", row["status"],
                             memory(row["peak_rss"]), seconds(row["realtime"])))  # fmt: skip
    schema = {"process": pl.String, "sample": pl.String, "status": pl.String,
              "peak_rss": pl.Int64, "realtime": pl.Float64}  # fmt: skip
    return pl.DataFrame(rows, schema=schema, orient="row")


def costs(frame: pl.DataFrame) -> pl.DataFrame:
    """One row per process that ran in some launch: tasks, samples, peak RSS, summed time."""
    return (
        frame.filter(pl.col("status") != CACHED)
        .group_by("process")
        .agg(tasks=pl.len(), samples=pl.col("sample").filter(pl.col("sample") != "").n_unique(),
             peak_rss=pl.col("peak_rss").max(), realtime=pl.col("realtime").sum())
        .sort(["peak_rss", "process"], descending=[True, False], nulls_last=True)
    )  # fmt: skip


def table(outdirs: Sequence[tuple[str, Path | None]]) -> Group | None:
    """The gather page's cost table over each pipeline's outdir that holds a trace, or None."""
    word = functools.partial(text.word, PAGE)
    rows, cached, total = [], 0, 0
    for pipeline, outdir in outdirs:
        frame = read(outdir) if outdir is not None else pl.DataFrame()
        if frame.is_empty():
            continue
        total += frame.height
        cached += frame.filter(pl.col("status") == CACHED).height
        for process, tasks, samples, peak, real in costs(frame).rows():
            mib = text.fill(PAGE, "trace.mib", mib=peak / 1024**2) if peak is not None else "-"
            rows.append((pipeline, process, f"{tasks:,}", f"{samples:,}", mib,
                         text.fill(PAGE, "trace.minutes", minutes=real / 60)))  # fmt: skip
    if not total:
        return None
    note = text.fill(PAGE, "trace.note", cached=cached, total=total)
    result = ResultTable(word("trace.title"), text.words(PAGE, "trace.columns"), tuple(rows),
                         numeric=(2, 3, 4, 5))  # fmt: skip
    return Group(word("trace.group"), (note, result))
