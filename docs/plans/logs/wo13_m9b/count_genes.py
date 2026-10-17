"""Package M9b: rows of each root's newest gather counts table, for the reader's docstring."""

from pathlib import Path

import polars as pl

from ctx.de.run import newest

for root, slug in (("/Users/rytis/ctx_breast", "tcga_brca"), ("/Users/rytis/ctx_tier1", "tier1")):
    path = newest(Path(root) / "outputs" / "01_gather", slug, "counts")
    rows = pl.scan_parquet(path).select(pl.len()).collect().item()
    print(f"{path.name}: {rows:,} genes")
