"""wo13-s2 R2.5: rows placed and the left-out share, R2.2's root against this package's.

Read from each root's own gather tables, never typed: the matrix's rows, the rows placed, the
Ensembl rows the annotation lacks, those placed through the ID history, the share of all
counts left out, the ten largest rows left out, and DE calls per contrast.
Usage: python r25_numbers.py <before root> <before slug> <after root> <after slug>
"""

import sys
from pathlib import Path

import polars as pl


def newest(folder: Path, pattern: str) -> Path:
    """The newest file of a pattern in a folder."""
    return sorted(folder.glob(pattern))[-1]


def numbers(root: Path, slug: str) -> None:
    """Print one root's placement and its DE calls."""
    gather = root / "outputs" / "01_gather"
    genes = pl.read_parquet(newest(gather, f"*_gather_{slug}_genes_v1.parquet"))
    left = pl.read_parquet(newest(gather, f"*_gather_{slug}_genes_left_out_v1.parquet"))
    counts = pl.read_parquet(newest(gather, f"*_gather_{slug}_counts_v1.parquet"))
    kept = float(counts.select(pl.exclude("gene_id").cast(pl.Float64)).sum().sum_horizontal()[0])
    lost = float(left["total"].sum())
    rows = genes.height + left.height
    print(
        f"== {root.name}: {genes.height:,} of {rows:,} rows placed, {left.height:,} left out, "
        f"{lost:,.0f} of {kept + lost:,.0f} counts left out = {100 * lost / (kept + lost):.4f} %"
    )
    print("left out by reason:", sorted(left.group_by("kind", "reason").len().rows()))
    print(
        "left out with an ID-history note:",
        left.filter(pl.col("detail").str.starts_with("ID history"))
        .group_by("detail")
        .len()
        .sort("len", descending=True)
        .head(5)
        .rows(),
    )
    page = newest(gather, f"*_gather_{slug}_overview_v1.html").read_text()
    for key in ("Placed on the annotation", "Placed through Ensembl", "Placed through HGNC"):
        start = page.find(key)
        print(key, page[start : start + 260].replace("\n", " ") if start >= 0 else "not on page")
    print("largest left out:")
    for row in (
        left.sort("total", descending=True)
        .head(10)
        .select("source_id", "reason", "detail", "total")
        .rows()
    ):
        print("  ", row)
    for results in sorted((root / "outputs" / "05_de").glob(f"*_de_{slug}_*_results_v1.parquet")):
        frame = pl.read_parquet(results)
        print(
            results.name,
            "tested",
            frame.height,
            "calls",
            sorted(frame.group_by("call").len().rows()),
        )


def main() -> None:
    """Both roots, then the rows the history placed, by the counts they carry."""
    before, after = (Path(sys.argv[i]).expanduser() for i in (1, 3))
    numbers(before, sys.argv[2])
    numbers(after, sys.argv[4])
    genes = pl.read_parquet(
        newest(after / "outputs" / "01_gather", f"*_gather_{sys.argv[4]}_genes_v1.parquet")
    )
    old = pl.read_parquet(
        newest(before / "outputs" / "01_gather", f"*_gather_{sys.argv[2]}_genes_v1.parquet")
    )
    gained = genes.join(old, on="source_id", how="anti")
    lost = old.join(genes, on="source_id", how="anti")
    print(f"rows gained {gained.height:,}, rows lost {lost.height:,}")
    changed = genes.join(old, on="source_id").filter(pl.col("gene_id") != pl.col("gene_id_right"))
    print(f"rows placed on another gene than before: {changed.height}")
    print("gained, first 10 by source ID:", gained.sort("source_id").head(10).rows())
    sod2 = genes.filter(pl.col("source_id").str.starts_with("ENSG00000112096"))
    print("ENSG00000112096 now:", sod2.rows())


if __name__ == "__main__":
    main()
