"""Package wo13-l: measure every module page of the scratch roots of one label.

Usage: <main env python> docs/plans/logs/wo13_l/env.py python docs/plans/logs/wo13_l/sweep.py
       <label>
The pages: each root's setup, overview, hub and glossary pages; the result goes to
docs/plans/logs/wo13_l/widths_<label>.json through measure.py.
"""

import subprocess
import sys
from pathlib import Path

WT = Path(__file__).absolute().parents[4]
label = sys.argv[1]
pages = sorted(
    str(p.relative_to(WT))
    for root in sorted((WT / "outputs" / "scratch_roots").glob(f"*_{label}"))
    for p in (root / "outputs").glob("*/*.html")
    if any(k in p.name for k in ("_overview_", "_hub_", "_environment_", "_glossary_"))
)
print(len(pages), "pages", flush=True)
out = WT / "docs" / "plans" / "logs" / "wo13_l" / f"widths_{label}.json"
sys.exit(subprocess.run([sys.executable, str(WT / "docs/plans/logs/wo13_l/measure.py"), str(out),
                         *pages], cwd=WT, check=False).returncode)  # fmt: skip
