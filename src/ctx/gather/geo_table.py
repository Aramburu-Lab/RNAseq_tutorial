r"""A sample table from a GEO series matrix, one row per sample, for ``ctx gather --counts``.

GEO's ``<GSE>_series_matrix.txt.gz`` lists each sample's accession, title, source and
``key: value`` characteristics and descriptions; an authors' matrix names its columns in a
description, often ``Library name: <column>``, which becomes the ``library_name`` column.
This writes them as a CSV the dataset file's ``[sdrf]`` joins on ``sample``, the column
names GEO's own count tables use, or on ``library_name``. Each value is written as GEO
gives it; the dataset file's ``[columns]`` map them to the course's levels. Standard library
only. Usage::

    python -m ctx.gather.geo_table <series_matrix.txt.gz> <table.csv> [--donor-pattern REGEX]

``--donor-pattern`` takes a regular expression with one group, matched against each title,
whose group is the donor: ``_(\d+)$`` reads 12 from "ANM_12". Without it each sample is its own
donor, which is right only for a study with one sample per individual.
"""

import argparse
import csv
import gzip
import re
import sys
from collections import Counter
from pathlib import Path

# The series matrix's per-sample lines that become columns, and their names in the table.
FIELDS = {"!Sample_geo_accession": "sample", "!Sample_title": "title",
          "!Sample_source_name_ch1": "source_name"}  # fmt: skip
CHARACTERISTICS = "!Sample_characteristics_ch1"
# Read as characteristics are; a line with no "key:" goes to the column DESCRIPTION.
DESCRIPTIONS = "!Sample_description"
DESCRIPTION = "description"
# Values tallied on screen per column; a column with more distinct values is named only.
SHOWN_LEVELS = 8


class TableError(ValueError):
    """A series matrix that cannot give a table; the message names every problem."""


def _key(text: str) -> str:
    """A characteristic's key as a column name: lowercase, words joined by underscores."""
    return re.sub(r"[^a-z0-9]+", "_", text.strip().lower()).strip("_") or "characteristic"


def read_series(path: Path) -> list[dict[str, str]]:
    """One dict per sample: accession, title, source and each characteristic by its key.

    Raises:
        TableError: If the file lists no sample accessions, or the lines differ in length.
    """
    fields: dict[str, list[str]] = {}
    traits: list[list[str]] = []
    with gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open() as h:
        for row in csv.reader(h, delimiter="\t"):
            if row and row[0] in FIELDS:
                fields[FIELDS[row[0]]] = row[1:]
            elif row and row[0] in (CHARACTERISTICS, DESCRIPTIONS):
                traits.append(row[1:])
    samples = fields.get("sample", [])
    if not samples:
        raise TableError(f"{path.name}: no !Sample_geo_accession line; not a series matrix")
    if any(len(v) != len(samples) for v in [*fields.values(), *traits]):
        raise TableError(f"{path.name}: the sample lines differ in length")
    records: list[dict[str, str]] = []
    for i in range(len(samples)):
        record = {name: values[i] for name, values in fields.items()}
        for line in traits:
            key, colon, value = line[i].partition(":")
            if not colon and key.strip():
                given = (record.get(DESCRIPTION), key.strip())
                record[DESCRIPTION] = "; ".join(filter(None, given))
            elif value or key:
                record[_key(key)] = value.strip()
        records.append(record)
    return records


def add_donor(rows: list[dict[str, str]], pattern: str | None) -> list[str]:
    """Set each row's donor from its title by the pattern, or to its own sample.

    Returns:
        The titles the pattern did not match, which leave their donors unset.
    """
    unmatched = []
    for row in rows:
        if pattern is None:
            row["donor"] = row["sample"]
            continue
        found = re.search(pattern, row.get("title", ""))
        if found is None or not found.groups():
            unmatched.append(row.get("title", ""))
            continue
        row["donor"] = found.group(1)
    return unmatched


def write(rows: list[dict[str, str]], out: Path) -> list[str]:
    """Write the table, sample first and donor second, every column any sample has.

    Returns:
        The column names, in the order written.
    """
    columns = ["sample", "donor"]
    for row in rows:
        columns += [k for k in row if k not in columns]
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, restval="")
        writer.writeheader()
        writer.writerows(rows)
    return columns


def tally(rows: list[dict[str, str]], columns: list[str]) -> list[str]:
    """One line per column: its values and how many samples hold each, or how many distinct."""
    lines = []
    for column in columns:
        counted = Counter(row.get(column, "") for row in rows)
        if len(counted) <= SHOWN_LEVELS:
            shown = "; ".join(f"{v!r} {n}" for v, n in sorted(counted.items()))
        else:
            shown = f"{len(counted)} distinct values"
        lines.append(f"{column}: {shown}")
    return lines


def main(argv: list[str] | None = None) -> int:
    """Write the table and print each column's values, for the dataset file's ``[columns]``."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("series", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--donor-pattern", default=None)
    args = parser.parse_args(argv)
    try:
        rows = read_series(args.series)
        if unmatched := add_donor(rows, args.donor_pattern):
            raise TableError(f"--donor-pattern matched no group in {len(unmatched)} titles, "
                             f"e.g. {unmatched[0]!r}")  # fmt: skip
    except (TableError, re.error, OSError) as error:
        print(error, file=sys.stderr)
        return 1
    columns = write(rows, args.out)
    print(f"{args.out}: {len(rows)} samples")
    print("\n".join(tally(rows, columns)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
