"""Write a table together with its ``_README.md``: provenance plus a data dictionary.

Provenance lives beside the table, never stamped inside it, so every CSV stays readable by a
plain CSV reader (CONVENTIONS.md, Outputs).
"""

import hashlib
import shlex
import sys
from collections.abc import Mapping, Sequence
from datetime import UTC, datetime
from pathlib import Path

import polars as pl
import polars.selectors as cs

from ctx import __version__
from ctx.io.naming import check
from ctx.runtime import source, versions

FORMATS = ("parquet", "csv")


def write_table(  # noqa: PLR0913 - keyword-only after the third, each a separate choice
    frame: pl.DataFrame,
    directory: Path,
    table_stem: str,
    *,
    columns: Mapping[str, tuple[str, str]],
    data_releases: Mapping[str, str],
    formats: Sequence[str] = FORMATS,
    csv_decimals: int | None = None,
) -> list[Path]:
    """Validate, write the data files, hash them, then write the README last.

    Args:
        frame: The table to publish.
        directory: Where the files go; created if missing, only after validation passes.
        table_stem: A stem from ``ctx.io.naming.stem()``.
        columns: For every column of ``frame``, its (unit, meaning). Use ``"none"`` as the
            unit of a label or identifier.
        data_releases: External resource to release identifier (GENCODE, GDC, ...); pass an
            empty mapping, deliberately, when the table uses none.
        formats: Any of ``"parquet"`` and ``"csv"``, written in this order.
        csv_decimals: The decimals the source carries; the CSV rounds every float column to
            them and writes a whole number without ``.0``, the parquet keeps full precision.
            Digits past the source's precision are arithmetic, not measurement.

    Returns:
        The data files, then the README.
    """
    validate(frame, table_stem, columns, formats)
    directory.mkdir(parents=True, exist_ok=True)
    written = []
    for fmt in formats:
        path = directory / f"{table_stem}.{fmt}"
        if fmt == "parquet":
            frame.write_parquet(path)
        elif csv_decimals is None:
            frame.write_csv(path)
        else:
            frame.with_columns(_shown(csv_decimals)).write_csv(path)
        written.append(path)
    digests = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in written}
    readme = directory / f"{table_stem}_README.md"
    readme.write_text(_readme(frame, table_stem, columns, data_releases, digests, csv_decimals))
    return [*written, readme]


def _shown(decimals: int) -> pl.Expr:
    """Every float column as the CSV writes it: rounded, and ``0`` rather than ``0.0``.

    A count matrix is mostly zeros, so the ``.0`` alone is a third of its bytes.
    """
    return cs.float().round(decimals).cast(pl.String).str.strip_suffix(".0")


def validate(
    frame: pl.DataFrame,
    table_stem: str,
    columns: Mapping[str, tuple[str, str]],
    formats: Sequence[str],
) -> None:
    """Raise before anything touches the disk.

    Args:
        frame: The table to publish.
        table_stem: Its stem.
        columns: (unit, meaning) per column.
        formats: Formats to write.

    Raises:
        ValueError: On an illegal stem, an unknown format, or a column without both a unit
            and a meaning.
    """
    if not check(f"{table_stem}.csv"):
        raise ValueError(f"not a legal output stem: {table_stem!r}")
    if not formats or set(formats) - set(FORMATS):
        raise ValueError(f"formats must be a non-empty subset of {FORMATS}, got {formats!r}")
    if missing := [c for c in frame.columns if c not in columns]:
        raise ValueError(f"no (unit, meaning) for columns {missing}")
    if extra := [c for c in columns if c not in frame.columns]:
        raise ValueError(f"described columns not in the frame: {extra}")
    if blank := [c for c, (unit, meaning) in columns.items() if not unit or not meaning]:
        raise ValueError(f"empty unit or meaning for columns {blank}")


def _read_back_types(path: Path) -> dict[str, str]:
    """Column types as a default reader of the written file infers them."""
    if path.suffix == ".parquet":
        schema = pl.read_parquet_schema(path)
    else:
        schema = dict(pl.scan_csv(path).collect_schema())
    return {name: str(dtype) for name, dtype in schema.items()}


def _readme(
    frame: pl.DataFrame,
    table_stem: str,
    columns: Mapping[str, tuple[str, str]],
    data_releases: Mapping[str, str],
    digests: Mapping[Path, str],
    csv_decimals: int | None = None,
) -> str:
    """Render the companion README."""
    created = datetime.now(UTC).isoformat(timespec="seconds").replace("+00:00", "Z")
    types = {p.suffix.lstrip("."): _read_back_types(p) for p in digests}
    lines = [
        f"# {table_stem}",
        "",
        "Written by `ctx`. Do not edit; rerun the command instead.",
        "",
        "| Field | Value |",
        "| --- | --- |",
        f"| Created (UTC) | {created} |",
        f"| Command | `{shlex.join(sys.argv)}` |",
        f"| ctx | {__version__}, {source()} |",
        f"| Directory | {next(iter(digests)).parent.resolve()} |",
        f"| Rows | {frame.height} |",
        f"| Columns | {frame.width} |",
        f"| Estimated in-memory size (bytes) | {frame.estimated_size()} |",
        *(
            [
                f"| CSV floats | rounded to {csv_decimals} decimals, whole numbers without .0; "
                "the parquet keeps every digit |"
            ]
            if csv_decimals is not None
            else []
        ),
        "",
        "## Files",
        "",
        "| File | Bytes | SHA-256 |",
        "| --- | --- | --- |",
        *(f"| `{p.name}` | {p.stat().st_size} | `{d}` |" for p, d in digests.items()),
        "",
        "## Data releases",
        "",
        *(
            ["| Resource | Release |", "| --- | --- |"]
            + [f"| {k} | {v} |" for k, v in data_releases.items()]
            if data_releases
            else ["None: this table uses no external data release."]
        ),
        "",
        "## Software versions",
        "",
        "| Package | Version |",
        "| --- | --- |",
        *(f"| {k} | {v} |" for k, v in versions().items()),
        "",
        "## Data dictionary",
        "",
        "| Column | Unit | Meaning | " + " | ".join(f"Type ({f})" for f in types) + " |",
        "| --- | --- | --- | " + " | ".join("---" for _ in types) + " |",
        *(
            f"| `{c}` | {columns[c][0]} | {columns[c][1]} | "
            + " | ".join(f"`{t[c]}`" for t in types.values())
            + " |"
            for c in frame.columns
        ),
        "",
    ]
    return "\n".join(lines)
