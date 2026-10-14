"""M1c curation: join a cohort's SDRF to its fetchngs runs and fill the curation columns.

The order is fixed and tested: read the SDRF as text; refuse a configured header that
repeats; collapse whitespace, then map values; join both ways and cross-check; derive read
length and batch; rank replicates, one number per library; compare the declared counts. Every
problem is collected, so one run lists them all.
"""

import csv
import io
from pathlib import Path

import polars as pl

from ctx.dataset import Dataset
from ctx.gather.own_fastq import table_rows

REPLICATE_KEYS = ("donor", "cell_type", "condition")
# One library: an ENA experiment holds one library, its runs the lanes or reruns of it, and
# ctx reads names a sample's lanes by the sample (D-063).
LIBRARY = "experiment_accession"
# (check, key, detail) of one curation problem.
Problem = tuple[str, str, str]
PROBLEM_COLUMNS = {
    "check": ("none", "The curation step that found the problem."),
    "key": ("none", "The run, header, value or group the problem concerns."),
    "detail": ("none", "What was found, and what was expected."),
}


class CurationError(ValueError):
    """Curation found problems; ``problems`` holds (check, key, detail) for each."""

    def __init__(self, problems: list[Problem]) -> None:
        """Keep the problems and summarise them in the message."""
        self.problems = problems
        lines = [f"{check}: {key}: {detail}" for check, key, detail in problems]
        super().__init__(f"{len(problems)} curation problems:\n  " + "\n  ".join(lines))


def read_table(path: Path) -> tuple[str, str]:
    """A student's own table as curation reads an SDRF: tab-separated text, and its SHA-256.

    A spreadsheet's export is taken as saved, by ``own_fastq.table_rows``, which ``ctx reads``
    reads the same table with: the byte-order mark and Windows line ends are removed, and a
    ``.csv`` is read with commas, quotes honoured, every other ending with tabs.

    Args:
        path: The table, as ``Dataset.table`` resolves it.

    Returns:
        The table as tab-separated text, and the SHA-256 of the file's bytes as saved.
    """
    rows, digest = table_rows(path)
    buffer = io.StringIO()
    csv.writer(buffer, delimiter="\t", lineterminator="\n").writerows(rows)
    return buffer.getvalue(), digest


def normalise(value: str) -> str:
    """Strip and collapse internal whitespace: '2nd  follow up' becomes '2nd follow up'."""
    return " ".join(value.split())


def rank_replicates(
    frame: pl.DataFrame, keys: tuple[str, ...], order: dict[str, tuple[str, ...]]
) -> pl.Series:
    """Number the samples within each key group: by the declared level order, then by run.

    Args:
        frame: One row per run, with ``run_accession``, the key columns and the ordered ones.
        keys: Columns defining a group, e.g. donor, cell type and condition.
        order: Column to its levels in order, e.g. time points; compared in the given order.

    Returns:
        A replicate number per row, from 1, in the frame's row order.
    """
    ranks = [
        pl.col(column).replace_strict(list(levels), list(range(len(levels))), return_dtype=pl.Int64)
        for column, levels in order.items()
    ]
    run = pl.col("run_accession").str.extract(r"(\d+)$").cast(pl.Int64)
    return (
        frame.with_row_index("_row")
        # Then by name: ctx reads' Illumina names all end in _001, so their numbers tie.
        .sort([*ranks, run, pl.col("run_accession")])
        .with_columns(replicate=(pl.int_range(pl.len()).over(list(keys)) + 1).cast(pl.Int32))
        .sort("_row")["replicate"]
    )


def library_replicates(frame: pl.DataFrame, order: dict[str, tuple[str, ...]]) -> pl.Series:
    """``rank_replicates`` over libraries: every run of one library takes the library's number.

    Each library is ranked once, by its first run, within the group its runs share; a library
    whose runs differ in a key or an ordered column is a library per value, ranked apart.

    Args:
        frame: One row per run, with ``run_accession``, ``LIBRARY``, the keys and ``order``.
        order: Column to its levels in order, as ``rank_replicates`` takes it.

    Returns:
        A replicate number per row, from 1, in the frame's row order.
    """
    columns = [LIBRARY, *REPLICATE_KEYS, *(c for c in order if c not in REPLICATE_KEYS)]
    run = pl.col("run_accession").str.extract(r"(\d+)$").cast(pl.Int64)
    first = frame.sort(run, "run_accession", nulls_last=True).unique(
        columns, keep="first", maintain_order=True
    )
    ranked = first.select(*columns, replicate=rank_replicates(first, REPLICATE_KEYS, order))
    return frame.join(ranked, on=columns, how="left", nulls_equal=True, maintain_order="left")[
        "replicate"
    ]


def _header(text: str, spec: Dataset) -> tuple[list[str], list[list[str]], list[Problem]]:
    """The SDRF header and rows, and a problem per configured header missing or repeated."""
    rows = list(csv.reader(io.StringIO(text), delimiter="\t"))
    header, body = rows[0], [r for r in rows[1:] if any(cell.strip() for cell in r)]
    assert spec.sdrf is not None  # noqa: S101 - curate() is called only with an SDRF
    wanted = [spec.sdrf.join[0], *(a for a, _ in spec.sdrf.cross_check)]
    wanted += [c.source for c in spec.columns]
    problems = [
        ("sdrf header", name, f"appears {header.count(name)} times; expected once")
        for name in dict.fromkeys(wanted)
        if header.count(name) != 1
    ]
    return header, body, problems


def _join(runs: pl.DataFrame, table: pl.DataFrame, spec: Dataset) -> list[Problem]:
    """Problems of the join in both directions, and of each cross-checked pair."""
    assert spec.sdrf is not None  # noqa: S101 - see _header
    key, column = spec.sdrf.join
    problems = [
        ("join", k, f"{n} SDRF rows share this {key}")
        for k, n in table.group_by(key).len().filter(pl.col("len") > 1).rows()
    ]
    sdrf_keys, run_keys = set(table[key]), set(runs[column])
    problems += [("join", k, f"SDRF row matches no {column}") for k in sorted(sdrf_keys - run_keys)]
    problems += [("join", k, f"run has no SDRF row ({key})") for k in sorted(run_keys - sdrf_keys)]
    joined = runs.join(table, left_on=column, right_on=key, how="inner")
    for sdrf_header, sheet_column in spec.sdrf.cross_check:
        mismatched = joined.filter(
            pl.col(sheet_column).map_elements(normalise, return_dtype=pl.String)
            != pl.col(sdrf_header)
        )
        problems += [
            ("cross-check", run, f"{sdrf_header} {a!r} != {sheet_column} {b!r}")
            for run, a, b in mismatched.select(column, sdrf_header, sheet_column).rows()
        ]
    return problems


def _map_values(table: pl.DataFrame, spec: Dataset) -> tuple[pl.DataFrame, list[Problem]]:
    """Curated columns from normalised SDRF values; an unmapped value is a problem."""
    problems = []
    for c in spec.columns:
        if c.values is None:
            table = table.with_columns(pl.col(c.source).alias(c.name))
            continue
        for value, n in table.group_by(c.source).len().sort(c.source).rows():
            if value not in c.values:
                problems.append(("value map", f"{c.name}: {value!r}", f"{n} rows; not in values"))
        mapped = pl.col(c.source).replace_strict(c.values, default=None, return_dtype=pl.String)
        table = table.with_columns(mapped.alias(c.name))
    return table, problems


def _batch(frame: pl.DataFrame, spec: Dataset) -> tuple[pl.DataFrame, list[Problem]]:
    """Read length from ENA counts, then batch from the declared bins; never rounded."""
    mates = pl.when(pl.col("library_layout") == "PAIRED").then(2).otherwise(1)
    length = pl.col("base_count").cast(pl.Float64) / pl.col("read_count").cast(pl.Float64) / mates
    frame = frame.with_columns(read_length=length)
    if not spec.batch_bins:
        return frame.with_columns(batch=pl.lit("")), []
    level = {n: name for name, lengths in spec.batch_bins.items() for n in lengths}
    problems = []
    batches = []
    for run, value in frame.select("run_accession", "read_length").rows():
        exact = value is not None and float(value).is_integer()
        batches.append(level.get(int(value)) if exact else None)
        if batches[-1] is None:
            problems.append(("batch", run, f"read length {value} nt is in no declared bin"))
    return frame.with_columns(batch=pl.Series(batches, dtype=pl.String)), problems


def _expect(frame: pl.DataFrame, spec: Dataset) -> list[Problem]:
    """Declared totals and per-value counts against the curated frame, value by value."""
    problems = []
    if spec.expect_runs is not None and frame.height != spec.expect_runs:
        problems.append(("expect", "runs", f"{frame.height} runs; expected {spec.expect_runs}"))
    for column, expected in spec.expect_counts.items():
        values = frame[column]
        if values.dtype == pl.Float64:
            values = values.map_elements(lambda v: f"{v:g}", return_dtype=pl.String)
        found = {str(k): n for k, n in values.value_counts().rows()}
        for value in sorted(set(found) | set(expected)):
            if found.get(value, 0) != expected.get(value, 0):
                detail = f"{found.get(value, 0)} runs; expected {expected.get(value, 0)}"
                problems.append(("expect", f"{column} = {value}", detail))
    return problems


def curate(runs: pl.DataFrame, sdrf_text: str, spec: Dataset) -> pl.DataFrame:
    """Curate the runs of one fetchngs outdir with the dataset's SDRF.

    Args:
        runs: ``fetchngs.read(outdir).runs``.
        sdrf_text: The SDRF, as text.
        spec: The dataset file, with an ``[sdrf]`` table.

    Returns:
        One row per run, whatever column the join names: ``run_accession``, the curated
        columns, ``read_length``, ``batch`` and ``replicate``.

    Raises:
        CurationError: With every problem found; nothing is returned in part.
    """
    assert spec.sdrf is not None  # noqa: S101 - the caller checks the dataset has an SDRF
    header, body, problems = _header(sdrf_text, spec)
    if problems:
        raise CurationError(problems)
    wanted = list(dict.fromkeys([spec.sdrf.join[0], *(a for a, _ in spec.sdrf.cross_check),
                                 *(c.source for c in spec.columns)]))  # fmt: skip
    table = pl.DataFrame(
        {name: [normalise(row[header.index(name)]) for row in body] for name in wanted}
    )
    table, problems = _map_values(table, spec)
    problems += _join(runs, table, spec)
    names = [c.name for c in spec.columns]
    key, column = spec.sdrf.join
    # The run is carried whatever the join: a table of ctx reads' samples joins on sample, and
    # each of a sample's lanes is still its own run.
    carried = list(dict.fromkeys(("run_accession", column, LIBRARY)))
    frame = runs.select(*carried, "library_layout", "base_count", "read_count").join(
        table.select(key, *names).unique(key, keep="first"),
        left_on=column,
        right_on=key,
        how="inner",
        maintain_order="left",
    )
    frame, batch_problems = _batch(frame, spec)
    problems += batch_problems
    for name, levels in spec.order.items():
        stray = sorted(set(frame[name].drop_nulls()) - set(levels))
        problems += [("order", f"{name}: {v!r}", "value missing from [order]") for v in stray]
    problems += _expect(frame, spec)
    if problems:
        raise CurationError(problems)
    replicate = library_replicates(frame, dict(spec.order))
    return frame.select("run_accession", *names, "read_length", "batch").with_columns(
        replicate=replicate
    )
