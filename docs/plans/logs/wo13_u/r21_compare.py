"""wo13-u R2.1: GSE130605 before (R2's root) and after (HGNC placed), from the files written."""

import json
from pathlib import Path

import polars as pl

SLUG = "r2_gse130605"
ROOTS = {"before": Path.home() / "ctx_r2_gse130605", "after": Path.home() / "ctx_u_gse130605"}
NAMED = ["IGKC", "IGHG1", "IGHA1", "MT-CO1", "MT-RNR2", "SNORD3A", "RPPH1"]


def table(root: Path, folder: str, name: str) -> pl.DataFrame:
    """One of the run's parquet tables, the newest day's: the after run crossed midnight."""
    return pl.read_parquet(sorted((root / "outputs" / folder).glob(f"*_{name}_v1.parquet"))[-1])


def numbers(root: Path) -> dict:
    """The gather, DE and sanity numbers a reader compares."""
    counts = table(root, "01_gather", f"gather_{SLUG}_counts")
    left = table(root, "01_gather", f"gather_{SLUG}_genes_left_out")
    kept = float(counts.drop("gene_id").sum().sum_horizontal()[0])
    held = float(left["total"].sum())
    de = table(root, "05_de", f"de_{SLUG}_tumour_normal_results")
    called = de.filter(pl.col("call") != "ns")
    sanity = table(root, "09_report", f"report_{SLUG}_sanity")
    named = (
        de.filter(pl.col("gene_name").is_in(NAMED))
        .select("gene_name", "log2FoldChange", "padj", "call")
        .sort("gene_name")
    )
    return {
        "rows placed": counts.height,
        "rows left out": left.height,
        "left-out share %": round(held / (kept + held) * 100, 4),
        "left out by reason": dict(left.group_by("reason").len().sort("reason").rows()),
        "genes tested": de.height,
        "genes called": called.height,
        "called up": called.filter(pl.col("call") == "up").height,
        "sanity verdicts": dict(sanity.select("rule", "rule_verdict").unique().sort("rule").rows()),
        "named genes": named.rows(),
    }


def main() -> None:
    """Print both runs' numbers, then the genes placed in one run and not the other."""
    found = {k: numbers(r) for k, r in ROOTS.items()}
    print(json.dumps(found, indent=1, default=str))
    before, after = (table(r, "01_gather", f"gather_{SLUG}_genes") for r in ROOTS.values())
    lost = before.join(after, on="source_id", how="anti")
    print("placed before and not after:", lost.height, lost.head(10).rows())
    gained = after.join(before, on="source_id", how="anti")
    print("placed after and not before:", gained.height)
    left = table(ROOTS["after"], "01_gather", f"gather_{SLUG}_genes_left_out")
    print("largest left out after:", left.sort("total", descending=True).head(10).rows())


if __name__ == "__main__":
    main()
