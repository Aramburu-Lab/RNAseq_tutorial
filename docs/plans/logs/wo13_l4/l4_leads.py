"""L.4: print the class figures' leads and captions from a root's M7 page, and count the genes."""

import html
import re
import sys
from pathlib import Path

import polars as pl

out = Path(sys.argv[1]) / "outputs" / "07_proteo"
page = sorted(out.glob("*_proteo_cptac_ccrcc_v1.html")) or sorted(out.glob("*overview*.html"))
text = html.unescape(re.sub(r"<[^>]+>", " ", page[-1].read_text()))
text = re.sub(r"\s+", " ", text)
for start in ("Median ρ over the tumours", "Median ρ over the normals", "Spearman ρ per gene over"):
    for m in re.finditer(re.escape(start) + r"[^.]*\.[^.]*\.", text):
        print("PAGE:", m.group(0))
genes = pl.read_parquet(sorted(out.glob("*_genes_v1.parquet"))[-1])
for prefix in ("", "normal_"):
    tested = genes.filter(pl.col(f"{prefix}tested") & (pl.col("gene_class") != ""))
    print(
        f"{prefix or 'tumour_'}tested {tested.height}, with a rho "
        f"{tested[f'{prefix}rho'].drop_nulls().len()}, with a padj "
        f"{tested[f'{prefix}rho_padj'].drop_nulls().len()}; BH set (any class) "
        f"{genes[f'{prefix}rho_p'].drop_nulls().len()}"
    )
for f in sorted(out.glob("*rho_by_class*_v1.csv")):
    print(f.name, pl.read_csv(f, comment_prefix="#").height, "rows")
