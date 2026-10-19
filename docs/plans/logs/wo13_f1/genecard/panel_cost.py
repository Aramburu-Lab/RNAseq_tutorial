"""wo13-f1 M9: the tier1 page's bytes without its Reproduce panels, everything else the same.

The page is laid out as before F.1, from the blocks built before the files (no panel), on the
same scratch root; the run with panels is moved aside first and page.py restored after.
Usage: wo13_f1_genecard_py.sh python docs/plans/logs/wo13_f1/genecard/panel_cost.py
"""

import subprocess
from pathlib import Path

W = Path(__file__).resolve().parents[5]
ROOT = W / "outputs" / "scratch_roots" / "tier1_m9"
ASIDE = W / "outputs" / "scratch_roots" / "tier1_m9_aside" / "11_genecard_panels"
PAGE = W / "src" / "ctx" / "genecard" / "page.py"
OLD, NEW = "built = _module_page(ctx, again,", "built = _module_page(ctx, drawn[:kept],"
GENES = "ISG15,IFI6,STAT1,OAS2,CD274,AIM2,GBP5"

(ROOT / "outputs" / "11_genecard").rename(ASIDE)
original = PAGE.read_text()
assert original.count(OLD) == 1
PAGE.write_text(original.replace(OLD, NEW))
try:
    done = subprocess.run(["python", "-m", "ctx", "--root", str(ROOT), "genecard", "--dataset",
                           "tier1", "--dataset-file", str(W / "conf/datasets/tier1.toml"),
                           "--offline", "--genes", GENES], capture_output=True, text=True,
                          check=False)  # fmt: skip
finally:
    PAGE.write_text(original)
print("exit", done.returncode)
for label, folder in (("without panels", ROOT / "outputs" / "11_genecard"), ("with panels", ASIDE)):
    (page,) = folder.glob("*_overview_v1.html")
    print(f"{label}: {page.stat().st_size:,} B")
