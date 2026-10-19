"""wo13-f1: run a real page's Reproduce panel, steps 2 to 5, in a clean folder; compare bytes.

Usage: wo13_f1_py.sh python docs/plans/logs/wo13_f1/redraw.py <page.html> <figure id> <folder>
Step 1 clones the course from GitHub and builds it with pixi; here CTX_REPO is this worktree
and ctx the environment's, the only line changed. Then the redrawn CSV and the PDF and SVG of
both themes are compared byte for byte with the module's files beside the page.
"""

import base64
import html
import re
import subprocess
import sys
from pathlib import Path

W = Path(__file__).resolve().parents[4]
page_path, figure, folder = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
page = page_path.read_text()
found = re.findall(r'<code id="ctx-repro-\d+">(.*?)</code>.*?base64,([^"]+)"', page, re.S)
scripts = [base64.b64decode(s).decode() for p, s in found
           if f"ctx figure {figure} " in html.unescape(p)]  # fmt: skip
print(f"panels on the page: {len(found)}; for {figure}: {len(scripts)}")
script = scripts[0]
head, rest = script.split("\n# 2. ", 1)
mkdir = next(line for line in head.splitlines() if line.startswith("mkdir -p repro_"))
local = ["#!/usr/bin/env bash", "set -euo pipefail", f'CTX_REPO="{W}"', mkdir]
folder.mkdir(parents=True, exist_ok=False)
(folder / "reproduce.sh").write_text("\n".join(local) + "\n\n# 2. " + rest)
done = subprocess.run(["bash", "reproduce.sh"], cwd=folder, capture_output=True, text=True,
                      check=False)  # fmt: skip
print("exit", done.returncode)
print(done.stdout[-1500:], done.stderr[-1500:], sep="\n")
(csv,) = folder.glob("repro_*/*.csv")
stem = csv.name[: -len(".csv")]
module_stem = next(p.name[: -len(".csv")] for p in page_path.parent.glob(f"*{stem[10:]}.csv"))
for end in (".csv", "_white.pdf", "_white.svg", "_black.pdf", "_black.svg"):
    new, old = csv.parent / f"{stem}{end}", page_path.parent / f"{module_stem}{end}"
    same = new.read_bytes() == old.read_bytes()
    print(f"{end:12} {'byte-equal' if same else 'DIFFERS'} {new.stat().st_size:>10,} B")
