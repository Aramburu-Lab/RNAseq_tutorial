#!/usr/bin/env python3
# /// codebase-meta
# exempt: P2C.3 one-off; writes the Step-1 record that tests/test_sanity_record.py reads
# ///
r"""Write the Step-1 record into deliverables/: the rule verdicts and the metrics behind them.

`ctx report` decides the dataset file's [sanity] rules on the newest M5 and M6 metrics and
writes outputs/09_report/<date>_report_<slug>_sanity_v<N>. The record tests read
deliverables/*_de_<slug>_sanity_v*.csv and *_de_<slug>_metrics_v*.csv and decide every rule
again from that metrics file alone, so it must hold both modules' metrics: the gsea.* ones
come from M6 only. This script copies the hub's verdicts byte for byte and writes the union of
the metrics, after checking that deciding again from the union gives the hub's verdicts. Any
failed check writes nothing (D-026).

Usage:
    pixi run python scripts/python/record_sanity.py --root ~/ctx_Non-cancerous-blood \
        --dataset Non-cancerous-blood --dataset-file conf/datasets/Non-cancerous-blood.toml
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import subprocess
import sys
from datetime import date
from pathlib import Path

import polars as pl

from ctx.dataset import load
from ctx.de.run import newest
from ctx.de.sanity import VERDICT_COLUMNS, decide, parse
from ctx.io.naming import MODULES, stem
from ctx.io.readme_companion import write_table

ROOT = Path(__file__).resolve().parents[2]
MODULE = "de"
UNION_COLUMNS = {
    "module": ("none", "The module whose metrics table holds the metric: de (M5), function (M6)."),
    "metric": ("none", "Name the sanity rules read, e.g. permanova.<pca>.<factor>.r2."),
    "value": ("varies", "The measured value; the unit follows from the metric's name."),
    "detail": ("none", "What the value was computed on, or why it is missing."),
}


def folder(outputs: Path, module: str) -> Path:
    """A module's output directory, numbered as ctx numbers it."""
    return outputs / f"{dict(MODULES)[module]:02d}_{module}"


def sha256(path: Path) -> str:
    """SHA-256 of a file's bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def metrics_union(outputs: Path, dataset: str) -> tuple[pl.DataFrame, dict[str, Path]]:
    """The newest M5 and M6 metrics tables, stacked; a metric in both is refused.

    Raises:
        ValueError: If a metric name appears in both tables.
    """
    paths = {m: newest(folder(outputs, m), dataset, "metrics") for m in ("de", "function")}
    frames = [
        pl.read_parquet(p)
        .select("metric", "value", "detail")
        .insert_column(0, pl.lit(module).alias("module"))
        for module, p in paths.items()
    ]
    union = pl.concat(frames)
    twice = union.filter(pl.col("metric").is_duplicated())["metric"].unique().sort().to_list()
    if twice:
        raise ValueError(f"metrics written by both M5 and M6: {twice}")
    return union, paths


def record_from_csv(text: str) -> dict[str, float]:
    """The metrics as tests/test_sanity_record.py reads them back from the CSV."""
    return {
        r["metric"]: float(r["value"]) if r["value"] else float("nan")
        for r in csv.DictReader(io.StringIO(text))
    }


def verdict_rows(text: str) -> list[tuple[str, str, str, str]]:
    """(rule, condition, verdict, rule_verdict) of every row of a verdict CSV."""
    return [
        (r["rule"], r["condition"], r["verdict"], r["rule_verdict"])
        for r in csv.DictReader(io.StringIO(text))
    ]


def check(
    union: pl.DataFrame, read: dict[str, Path], hub_parquet: Path, dataset_file: Path
) -> list[str]:
    """Every reason the record would not match what the hub decided; empty when it matches."""
    problems = []
    hub_csv = hub_parquet.with_suffix(".csv")
    if pl.read_parquet(hub_parquet).write_csv() != hub_csv.read_text():
        problems.append(f"{hub_csv.name} is not the CSV of {hub_parquet.name}")
    readme = hub_parquet.with_name(f"{hub_parquet.stem}_README.md").read_text()
    for path in read.values():
        if f"| sha256:{path.name} | {sha256(path)} |" not in readme:
            problems.append(f"{hub_parquet.name} was not decided on {path.name} as it is now")
    decided = decide(record_from_csv(union.write_csv()), parse(load(dataset_file).sanity))
    again = [(r["rule"], r["condition"], r["verdict"], r["rule_verdict"]) for r in decided]
    if again != verdict_rows(hub_csv.read_text()):
        problems.append("deciding the rules again from the metrics union differs from the hub")
    return problems


def rule_commit(dataset_file: Path) -> str:
    """The newest commit that registered the dataset file's rules, or 'none found'.

    Followed across a rename (D-114), whose commit adds the line under the new path.
    """
    found = subprocess.run(
        ["git", "-C", str(ROOT), "log", "-1", "--follow", "--format=%H", "-G", "^registered = ",
         "--", str(dataset_file.resolve())],
        capture_output=True, text=True, check=False,
    ).stdout.strip()  # fmt: skip
    return found or "none found"


def write(root: Path, dataset: str, dataset_file: Path, out: Path) -> list[Path]:
    """Check, then write the two record CSVs with their READMEs; nothing on any failure.

    Raises:
        FileExistsError: If out already holds a record for this dataset.
        ValueError: If a check fails; the message lists every problem.
    """
    if existing := sorted(out.glob(f"*_{MODULE}_{dataset}_sanity_v*.csv")):
        raise FileExistsError(f"a record exists: {existing[0]}; a record is written once")
    outputs = root / "outputs"
    union, read = metrics_union(outputs, dataset)
    hub_parquet = newest(folder(outputs, "report"), dataset, "sanity")
    if problems := check(union, read, hub_parquet, dataset_file):
        raise ValueError("no record written:\n  " + "\n  ".join(problems))
    releases = {
        "hub verdicts": hub_parquet.with_suffix(".csv").name,
        f"sha256:{hub_parquet.with_suffix('.csv').name}": sha256(hub_parquet.with_suffix(".csv")),
        **{f"sha256:{p.name}": sha256(p) for p in read.values()},
        f"sha256:{dataset_file.name}": sha256(dataset_file),
        "rule commit": rule_commit(dataset_file),
    }
    day = date.today()
    csv_only = {"data_releases": releases, "formats": ("csv",)}
    written = write_table(union, out, stem(MODULE, f"{dataset}_metrics", 1, day=day),
                          columns=UNION_COLUMNS, **csv_only)  # fmt: skip
    verdicts = pl.read_parquet(hub_parquet)
    written += write_table(verdicts, out, stem(MODULE, f"{dataset}_sanity", 1, day=day),
                           columns=VERDICT_COLUMNS, **csv_only)  # fmt: skip
    if sha256(written[2]) != sha256(hub_parquet.with_suffix(".csv")):
        for path in written:
            path.unlink()
        raise ValueError("no record written: the verdict CSV is not byte-identical to the hub's")
    return written


def main(argv: list[str] | None = None) -> int:
    """Write the record, or print why not and return 1."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--dataset-file", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=ROOT / "deliverables")
    args = parser.parse_args(argv)
    try:
        written = write(args.root, args.dataset, args.dataset_file, args.out)
    except (FileExistsError, ValueError) as error:
        print(error, file=sys.stderr)
        return 1
    for path in written:
        print(path, sha256(path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
