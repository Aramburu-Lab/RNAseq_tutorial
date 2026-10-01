"""MultiQC, the one parser of tool QC logs (CONVENTIONS.md, Statistics rules): run and read.

MultiQC runs by subprocess from the `reads` environment (D-027), never imported, into a
module's `multiqc/` folder, which keeps MultiQC's own file names (naming lint exemption, D-028).
Its per-module tables are read with Polars as text, so no value is reinterpreted on the way
in; a table's sample names are mapped to the course's units (runs or samples) by an explicit
alias table, and one name that maps to nothing stops the run with every such name listed.
"""

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from ctx import tools
from ctx.qc import QcError

# The MultiQC modules M2 and M3 read; each runs only where its tool's output is among the inputs.
MODULES = ("fastqc", "fastp", "seqkit", "trim_galore", "salmon", "star", "hisat2",
           "featurecounts", "samtools")  # fmt: skip
DATA = "multiqc_data"


@dataclass(frozen=True)
class MultiqcRun:
    """One MultiQC run: where its tables are, and the version that wrote them."""

    data: Path
    version: str


def run(
    inputs: Sequence[Path], dest: Path, *, config: Path | None = None,
    modules: Sequence[str] = MODULES,
) -> MultiqcRun:  # fmt: skip
    """Run MultiQC on tool outputs and return where its tables are.

    Args:
        inputs: Folders or files holding the tools' outputs.
        dest: MultiQC's output folder, written with `--force`.
        config: A generated MultiQC config (`ctx.qc.registry`), passed with `-c`.
        modules: MultiQC modules to run.

    Returns:
        The run's data folder and MultiQC's version.

    Raises:
        QcError: If MultiQC fails, or writes no data folder.
    """
    # MultiQC runs in dest, so every path relative to the caller would point nowhere: with a
    # relative -o MultiQC exits 1, with a relative -c it exits 2 (measured, MultiQC 1.35).
    dest = dest.resolve()
    dest.mkdir(parents=True, exist_ok=True)
    args = ["--force", "--no-version-check", "--no-ansi", "-o", str(dest)]
    args += ["-c", str(config.resolve())] if config else []
    args += [f"--module={m}" for m in modules]
    try:
        tools.run("multiqc", [*args, *(str(Path(p).resolve()) for p in inputs)], cwd=dest,
                  log=dest / "multiqc_run.log")  # fmt: skip
    except tools.ToolError as error:
        raise QcError(str(error)) from error
    if not (dest / DATA).is_dir():
        raise QcError(f"MultiQC wrote no {DATA}/ in {dest}; see {dest / 'multiqc_run.log'}")
    return MultiqcRun(dest / DATA, tools.TOOLS["multiqc"].version)


def table(run: MultiqcRun, module: str) -> pl.DataFrame:
    """A module's table, `multiqc_<module>.txt`, all columns as text, MultiQC's "None" as null.

    Raises:
        QcError: If MultiQC wrote no table for the module: it found none of its inputs.
    """
    path = run.data / f"multiqc_{module}.txt"
    if not path.is_file():
        raise QcError(f"MultiQC found no {module} input: {path.name} was not written")
    frame = pl.read_csv(path, separator="\t", infer_schema=False)
    return frame.with_columns(pl.all().replace("None", None))


def plot_table(run: MultiqcRun, name: str) -> pl.DataFrame:
    """A plot's data, `<name>.txt` (e.g. `star_gene_counts_Unstranded`), every value a number.

    A module's `multiqc_<module>.txt` holds only what it parsed for the table; some totals
    exist only as plot data, such as STAR's reads per gene-count category.

    Raises:
        QcError: If MultiQC wrote no such plot, or a value is not a number.
    """
    path = run.data / f"{name}.txt"
    if not path.is_file():
        raise QcError(f"MultiQC wrote no {path.name}: its module found none of its inputs")
    frame = pl.read_csv(path, separator="\t", infer_schema=False)
    try:
        return numeric(frame, [c for c in frame.columns if c != "Sample"])
    except pl.exceptions.InvalidOperationError as error:
        raise QcError(f"{path.name} holds a value that is not a number: {error}") from error


def numeric(frame: pl.DataFrame, columns: Iterable[str]) -> pl.DataFrame:
    """The table with the named columns cast to Float64; a value that is not a number raises."""
    return frame.with_columns(pl.col(list(columns)).cast(pl.Float64, strict=True))


def map_names(names: Iterable[str], aliases: Mapping[str, str]) -> dict[str, str]:
    """Each MultiQC sample name to the unit it describes.

    Args:
        names: Sample names as a MultiQC table lists them.
        aliases: Every name a unit may appear under (its sample name, a FASTQ stem, the stem
            with `_1`/`_2` or `_trimmed_1`/`_2`), mapped to the unit.

    Returns:
        Name to unit.

    Raises:
        QcError: If any name maps to no unit; the message lists all of them.
    """
    names = list(names)
    if missing := sorted({n for n in names if n not in aliases}):
        raise QcError(f"MultiQC sample names that match no sample: {missing}")
    return {n: aliases[n] for n in names}


def by_unit(
    frame: pl.DataFrame, ranks: Mapping[str, tuple[str, int]], key: str, keep_below: int
) -> pl.DataFrame:
    """One row per unit from MultiQC's rows: each column from the best-ranked row that has it.

    A sample appears under several names (a run's own row, its read-1 and read-2 rows, a
    pair's summed row); ``ranks`` maps each name to (unit, rank), and rows ranked
    ``keep_below`` or worse are never read.

    Raises:
        QcError: If a MultiQC name is not in ``ranks``.
    """
    map_names(frame["Sample"], {name: name for name in ranks})
    keyed = (
        frame.with_columns(
            pl.col("Sample").replace_strict({n: u for n, (u, _) in ranks.items()}).alias(key),
            _rank=pl.col("Sample").replace_strict({n: r for n, (_, r) in ranks.items()}),
        )
        .filter(pl.col("_rank") < keep_below)
        .sort("_rank")
    )
    return keyed.group_by(key, maintain_order=True).agg(
        pl.exclude("Sample", "_rank", key).drop_nulls().first()
    )


def require_sources(run: MultiqcRun, module: str, expected: Iterable[str],
                    aliases: Mapping[str, str]) -> None:  # fmt: skip
    """Stop unless every expected unit has a source file in the module (D-022).

    A sample whose report the cluster copy missed must stop `ctx qc`, not draw a chart
    without it.

    Raises:
        QcError: If a unit has no source in the module; the message lists all of them.
    """
    sources = pl.read_csv(run.data / "multiqc_sources.txt", separator="\t", infer_schema=False)
    seen = {aliases.get(n) for n in sources.filter(pl.col("Module").str.to_lowercase()
            .str.replace_all(r"[^a-z]", "") == module.replace("_", ""))["Sample Name"]}  # fmt: skip
    if missing := sorted(set(expected) - seen):
        raise QcError(f"no {module} report for {missing}")
