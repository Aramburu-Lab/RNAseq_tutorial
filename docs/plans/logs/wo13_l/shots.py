"""Package wo13-l: PNGs of each figure card of a page, for an eye check, through shot.mjs.

Usage: <main env python> docs/plans/logs/wo13_l/shots.py <page.html> <width> <scheme> <tag>
Writes outputs/wo13_l_shots/<tag>_<n>.png, one per figure card.
"""

import json
import re
import subprocess
import sys
from pathlib import Path

BIN = "/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin"
WT = Path(__file__).absolute().parents[4]
page, width, scheme, tag = sys.argv[1:]
html = Path(page).read_text(encoding="utf-8")
n = len(re.findall(r'<figure class="fig"', html))
out = WT / "outputs" / "wo13_l_shots"
out.mkdir(parents=True, exist_ok=True)
shots = [{"url": Path(page).absolute().as_uri(), "width": int(width), "scheme": scheme,
          "act": "new Promise(r => setTimeout(r, 2500))",
          "selector": f"figure.fig:nth-of-type(1) >> {i}", "out": str(out / f"{tag}_{i}.png")}
         for i in range(n)]  # fmt: skip
for s in shots:
    i = int(s["selector"].rsplit(" ", 1)[1])
    s["selector"] = "figure.fig"
    s["act"] = (f"(async () => {{ await new Promise(r => setTimeout(r, 2500)); "
                f"const all = document.querySelectorAll('figure.fig'); "
                f"all.forEach((f, j) => {{ if (j !== {i}) f.classList.remove('shot'); }}); "
                f"all[{i}].classList.add('shot'); }})()")  # fmt: skip
    s["selector"] = "figure.fig.shot"
plan = {"chromium": "/Applications/Chromium.app/Contents/MacOS/Chromium", "shots": shots}
done = subprocess.run([f"{BIN}/node", str(WT / "docs/plans/logs/wo13_l/shot.mjs")],
                      input=json.dumps(plan), capture_output=True, text=True,
                      check=False)  # fmt: skip
print(done.stdout[-2000:], done.stderr[-2000:])
