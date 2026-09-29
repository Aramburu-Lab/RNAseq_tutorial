#!/usr/bin/env python3
# /// codebase-meta
# exempt: instructor-side, once per cohort; picks the Track A teaching subset from a gather sheet
# ///
r"""Pick the Track A teaching subset: one donor per stratum, every listed cell type from each.

Track A runs FastQC, fastp, STAR and Salmon on a few runs subsampled to 1M reads (TODO.md,
P3.0). The runs must be chosen before any count is seen, from metadata alone, by a rule written
down here. Every cohort value is an argument, so the script holds none (test_repo's
test_no_cohort_identifiers).

The rule, stratum by stratum in the order given:
- a candidate donor has exactly one run of each listed cell type in the stratum's condition
  and read length, after the filters and the minimum read count;
- each cell type's isolation method is the same in every stratum: the most common method among
  the candidates that every stratum offers; if none is shared the script stops and says so;
- a donor chosen for one stratum is not a candidate for a later one;
- candidates are ranked with an empty sex last, because the sex check (M2) cannot be taught on
  a donor with none; then the sex chosen least so far; then by how close their mean read count
  is to the stratum's median; then by donor.

Usage:
    pixi run python scripts/python/pick_teaching_subset.py --sheet <gather samplesheet parquet> \
        --stratum <condition>:<read length> ... --cell-types <a>,<b> --filter <column>=<value> ...
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

import polars as pl


class PickError(ValueError):
    """The rule cannot be met; the message names the stratum or cell type that fails it."""


def candidates(
    sheet: pl.DataFrame, condition: str, read_length: float, cell_types: list[str]
) -> pl.DataFrame:
    """Donors with exactly one run of each cell type in the stratum, one row per run."""
    rows = sheet.filter(
        (pl.col("condition") == condition)
        & (pl.col("read_length") == read_length)
        & pl.col("cell_type").is_in(cell_types)
    )
    per = rows.group_by("donor", "cell_type").len()
    complete = (
        per.filter(pl.col("len") == 1)
        .group_by("donor")
        .agg(pl.col("cell_type").n_unique().alias("n"))
        .filter(pl.col("n") == len(cell_types))["donor"]
    )
    return rows.filter(pl.col("donor").is_in(complete.implode()))


def shared_isolation(pools: dict[str, pl.DataFrame], cell_types: list[str]) -> dict[str, str]:
    """Per cell type, the most common isolation method that every stratum's candidates offer.

    Raises:
        PickError: If a cell type has no method common to every stratum.
    """
    chosen = {}
    for cell in cell_types:
        offered = [set(p.filter(pl.col("cell_type") == cell)["isolation"]) for p in pools.values()]
        common = set.intersection(*offered) if offered else set()
        if not common:
            detail = {s: sorted(o) for s, o in zip(pools, offered, strict=True)}
            raise PickError(f"{cell}: no isolation method in every stratum: {detail}")
        counts = Counter(
            v for p in pools.values() for v in p.filter(pl.col("cell_type") == cell)["isolation"]
        )
        chosen[cell] = max(sorted(common), key=lambda method: counts[method])
    return chosen


def rank(pool: pl.DataFrame, sexes: Counter) -> list[str]:
    """Donors of one stratum in the rule's order: a known sex, the rarer sex, depth, donor."""
    depth = pool.group_by("donor").agg(pl.col("read_count").mean(), pl.col("sex").first())
    median = float(depth["read_count"].median())

    def key(row: dict) -> tuple:
        sex = row["sex"] or ""
        return (sex == "", sexes[sex], abs(row["read_count"] - median), row["donor"])

    return [r["donor"] for r in sorted(depth.iter_rows(named=True), key=key)]


def pick(
    sheet: pl.DataFrame, strata: list[tuple[str, float]], cell_types: list[str]
) -> pl.DataFrame:
    """The chosen runs, one donor per stratum.

    Raises:
        PickError: If a stratum has no donor left who meets the rule.
    """
    pools = {f"{c}:{r:g}": candidates(sheet, c, r, cell_types) for c, r in strata}
    methods = shared_isolation(pools, cell_types)
    chosen, sexes, used = [], Counter(), set()
    for label, pool in pools.items():
        matched = pool.filter(
            pl.struct("cell_type", "isolation").map_elements(
                lambda s: methods[s["cell_type"]] == s["isolation"], return_dtype=pl.Boolean
            )
        )
        whole = matched.group_by("donor").len().filter(pl.col("len") == len(cell_types))["donor"]
        free = pl.col("donor").is_in(whole.implode()) & ~pl.col("donor").is_in(list(used))
        matched = matched.filter(free)
        if matched.is_empty():
            raise PickError(f"{label}: no donor left with one run of each of {cell_types}")
        donor = rank(matched, sexes)[0]
        runs = matched.filter(pl.col("donor") == donor).with_columns(stratum=pl.lit(label))
        used.add(donor)
        sexes[runs["sex"][0] or ""] += 1
        chosen.append(runs)
    return pl.concat(chosen).sort("stratum", "cell_type")


def main(argv: list[str] | None = None) -> int:
    """Print the chosen runs and their accessions, one per line, for the fetchngs ID list."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sheet", type=Path, required=True, help="A gather samplesheet parquet.")
    parser.add_argument("--stratum", action="append", required=True, help="condition:read_length")
    parser.add_argument("--cell-types", required=True, help="Comma-separated cell types.")
    parser.add_argument("--filter", action="append", default=[], help="column=value; empty allowed")
    parser.add_argument("--min-reads", type=int, default=0)
    args = parser.parse_args(argv)
    sheet = pl.read_parquet(args.sheet).filter(pl.col("read_count") >= args.min_reads)
    for rule in args.filter:
        column, value = rule.split("=", 1)
        sheet = sheet.filter(pl.col(column).cast(pl.String) == value)
    strata = [(c, float(r)) for c, r in (s.rsplit(":", 1) for s in args.stratum)]
    try:
        runs = pick(sheet, strata, args.cell_types.split(","))
    except PickError as error:
        print(error, file=sys.stderr)
        return 1
    shown = ["stratum", "donor", "cell_type", "isolation", "sex", "read_count", "run_accession"]
    print(runs.select(shown).write_csv(separator="\t"), end="")
    print(f"fastq_bytes total: {runs['fastq_bytes'].cast(pl.Int64).sum():,}")
    print("\n".join(runs["run_accession"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
