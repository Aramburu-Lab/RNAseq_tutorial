"""M8's F3 cut-offs on a real page, at 1,280 and 375 px in both themes, against classify.

python m8_real.py <08_ncrna folder> <dataset> <downloads dir>
"""

import sys
from pathlib import Path

import polars as pl

WT = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(WT / "tests")]

from test_export_browser import _browse  # noqa: E402
from test_ncrna_views_browser import SETTINGS, STRATA  # noqa: E402

from ctx.de.stats import classify  # noqa: E402

out, dataset, downloads = Path(sys.argv[1]).absolute(), sys.argv[2], Path(sys.argv[3]).absolute()
downloads.mkdir(parents=True, exist_ok=True)
(page,) = sorted(out.glob(f"*_ncrna_{dataset}_overview_v1.html"))[-1:]
frame = pl.read_parquet(sorted(out.glob(f"*_{dataset}_results_by_biotype_v1.parquet"))[-1])
SIZES = (("1280_light", 1280, "light"), ("1280_dark", 1280, "dark"),
         ("375_light", 375, "light"), ("375_dark", 375, "dark"))  # fmt: skip
act = STRATA.replace("SETTINGS", repr([list(s) for s in SETTINGS]))
tasks = [{"name": n, "url": page.as_uri(), "width": w, "height": 900, "scheme": s, "act": act}
         for n, w, s in SIZES]  # fmt: skip
seen = _browse(tasks, downloads, timeout=1800)
print(f"page {page.name} {page.stat().st_size:,} B")
for name, _, _ in SIZES:
    acted = seen[name]["acted"]
    equal = 0
    for (basis, cutoff, fc), drawn in zip(SETTINGS, acted["settings"], strict=True):
        calls = frame.select("gene_name", "biotype", call=classify(
            pl.col(basis), pl.col("log2FoldChange"), float(cutoff), float(fc)))  # fmt: skip
        same = all(sorted(drawn["traces"][b]) == sorted(
            calls.filter((pl.col("biotype") == b) & (pl.col("call") != "ns"))["gene_name"]
            .fill_null("").to_list()) for b in drawn["traces"] if b != "ns")  # fmt: skip
        equal += same
        print(f"  {name} {basis} {cutoff} {fc}: {drawn['summary']}")
    print(f"{name}: {equal} of {len(SETTINGS)} settings equal classify gene for gene; box "
          f"{acted['box']} scroll {acted['scroll']}; reset moved {acted['reset']['moved']}; "
          f"errors {seen[name]['errors']}")  # fmt: skip
