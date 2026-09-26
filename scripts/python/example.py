#!/usr/bin/env python3
# /// codebase-meta
# name: example
# order: 20
# summary: Keep rows at or under a length cutoff, before anything expensive runs.
# entrypoint: true
# inputs:
#   - path: settings/samples.tsv
#     external: true
#     desc: Delimited table with an identifier column and a length column.
# outputs:
#   - path: outputs/20_example/out_example_kept.parquet
#     desc: One row per surviving record, carrying its identifier and length.
# options:
#   - flag: --max-length
#     desc: Longest value of the length column to keep.
#     default: 150
#   - flag: --length-column
#     desc: Column holding the length to threshold on.
#     default: length
#   - flag: --id-column
#     desc: Column holding the record identifier.
#     default: id
# next: final
# ///
"""Prefilter a table before anything expensive touches it.

Cutting rows here is the cheapest speed-up available anywhere in a pipeline:
every later step pays for whatever this one lets through.

This file is a worked example, not a fixture. Rename it to what it does, and
change `name:` and the output paths in the block above to match — `example` is
not a name, and the run-order table is only as useful as the names in it.

Note the shape of the write: `labcore.io.write_table` refuses a column it has no
description for, so a script cannot pass an unknown schema straight through. The
step therefore selects the columns it means to publish. That constraint is the
point — it is what stops a data dictionary from rotting.

Usage:
    pixi run python scripts/python/example.py -i settings/samples.tsv
    -o outputs/20_example/out_example_kept.parquet --max-length 150
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import polars as pl
from labcore.cli import setup_logging
from labcore.io import write_table

log = setup_logging("cancer_transcriptomics_teaching.example")

COLUMN_DESCRIPTIONS = {
    "id": "Record identifier, carried over unchanged from --id-column.",
    "length": "Value of --length-column for this record, at or below --max-length.",
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Read the command line.

    Args:
        argv: Argument list without the program name; `sys.argv[1:]` when None.

    Returns:
        The parsed arguments.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("-i", "--input", required=True, type=Path, help="Delimited input table.")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Parquet to write.")
    parser.add_argument("--max-length", type=int, default=150, help="Longest length to keep.")
    parser.add_argument("--length-column", default="length", help="Column to threshold on.")
    parser.add_argument("--id-column", default="id", help="Column holding the identifier.")
    parser.add_argument("--separator", default="\t", help="Input field separator.")
    return parser.parse_args(argv)


def keep_short(
    frame: pl.DataFrame, *, id_column: str, length_column: str, limit: int
) -> pl.DataFrame:
    """Select the records at or under a length limit.

    Args:
        frame: Input table.
        id_column: Column holding the record identifier.
        length_column: Column holding the length to threshold on.
        limit: Longest length to keep, inclusive.

    Returns:
        A two-column frame of `id` and `length`, in input order.

    Raises:
        KeyError: If either named column is absent, listing what the file has.
    """
    for column in (id_column, length_column):
        if column not in frame.columns:
            raise KeyError(f"no column '{column}'; the file has {', '.join(frame.columns)}")
    return frame.filter(pl.col(length_column) <= limit).select(
        pl.col(id_column).alias("id"),
        pl.col(length_column).alias("length"),
    )


def main(argv: list[str] | None = None) -> int:
    """Filter the input table and write the survivors.

    Args:
        argv: Argument list without the program name; `sys.argv[1:]` when None.

    Returns:
        0 on success, 1 when the input is missing a required column.
    """
    args = parse_args(argv)
    frame = pl.read_csv(args.input, separator=args.separator)
    try:
        kept = keep_short(
            frame,
            id_column=args.id_column,
            length_column=args.length_column,
            limit=args.max_length,
        )
    except KeyError as exc:
        log.error("%s: %s", args.input, exc)
        return 1

    write_table(kept, args.output, descriptions=COLUMN_DESCRIPTIONS)
    share = 100 * len(kept) / max(len(frame), 1)
    log.info("kept %d of %d rows (%.1f%%) -> %s", len(kept), len(frame), share, args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
