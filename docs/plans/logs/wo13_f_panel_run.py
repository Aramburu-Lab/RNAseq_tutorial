"""wo13-f: run one Reproduce panel of a real page, steps 2 to 5, in a clean folder.

Usage: wo13_f_py.sh python docs/plans/logs/wo13_f_panel_run.py <page.html> <figure id> <folder>
Step 1 clones the course from GitHub and builds it with pixi; here CTX_REPO is this worktree
and ctx the environment's, which is the only line changed. Prints the script's output.
"""

import base64
import html
import re
import subprocess
import sys
from pathlib import Path

W = Path(__file__).resolve().parents[3]
page, figure, folder = Path(sys.argv[1]).read_text(), sys.argv[2], Path(sys.argv[3])
found = re.findall(r'<code id="ctx-repro-\d+">(.*?)</code>.*?base64,([^"]+)"', page, re.S)
(script,) = [
    base64.b64decode(s).decode() for p, s in found if f"ctx figure {figure} " in html.unescape(p)
]
head, rest = script.split("\n# 2. ", 1)
environment = head.splitlines()
mkdir = next(line for line in environment if line.startswith("mkdir -p repro_"))
local = ["#!/usr/bin/env bash", "set -euo pipefail", f'CTX_REPO="{W}"', mkdir]
folder.mkdir(parents=True, exist_ok=True)
(folder / "reproduce.sh").write_text("\n".join(local) + "\n\n# 2. " + rest)
done = subprocess.run(["bash", "reproduce.sh"], cwd=folder, capture_output=True, text=True,
                      check=False)  # fmt: skip
print(done.returncode, done.stdout[-3000:], done.stderr[-3000:], sep="\n")
