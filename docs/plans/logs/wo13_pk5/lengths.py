"""PK.5: the ORA bars' lengths at -log10 p and -log10 FDR on Step-1 (read only).

Two cases: the static figure's bars (the build's ORA table of the featured contrast) and
Rytis's 6 genes as M6's live panel tests them (``live.listed``). Each row: family, set,
-log10 p (D-100's length), -log10 FDR (D-102's) and the FDR itself.
"""

import math
import tempfile
from pathlib import Path

import polars as pl

from ctx.function import lite, live

ROOT = "/Users/rytis/ctx_tier1"
SIX = ("MMP8", "NECAB1", "PLAC8", "ZDHHC19", "LAIR1", "CD300E")
results = f"{ROOT}/outputs/05_de/2026-10-16_de_tier1_ndn_gcsf_results_v1.parquet"
table = f"{ROOT}/outputs/06_function/2026-10-20_function_tier1_ndn_gcsf_ora_v1.parquet"
sets = lite.read_gmt(Path(f"{ROOT}/resources/msigdb/h.all.v2026.1.Hs.symbols.gmt"))
with tempfile.TemporaryDirectory() as folder:
    listed = Path(folder) / "six.txt"
    listed.write_text("\n".join(SIX) + "\n")
    cases = {
        "static, featured contrast": pl.read_parquet(table),
        "live, Rytis's 6 genes": live.listed(listed, pl.read_parquet(results), sets),
    }
for name, frame in cases.items():
    drawn = live.bars(frame)
    print(f"{name}: {drawn.height} bars, {int((~drawn['passes']).sum())} grey")
    for r in drawn.iter_rows(named=True):
        p, fdr = -math.log10(r["p"]), -math.log10(r["fdr"])
        print(f"  {r['series']:<13} {r['gene_set']:<45} -log10 p {p:8.4f}  "
              f"-log10 FDR {fdr:.3g}  FDR {r['fdr']!r}")  # fmt: skip
print(f"threshold -log10(0.05) = {-math.log10(live.ALPHA):.4f}")
