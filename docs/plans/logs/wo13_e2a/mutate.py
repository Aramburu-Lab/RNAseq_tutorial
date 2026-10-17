"""Package wo13-e2a: plant each mutation in place, run its tests, restore, log CAUGHT or SURVIVED.

Run: python docs/plans/logs/wo13_e2a/mutate.py <mutations.json> <log> [label ...]
Each mutation: {"label", "file", "old", "new", "tests"}; each old must occur once in the file.
The source is restored whatever happens, and the run writes no bytecode (lessons.md). The
tests run with this worktree's src as an absolute PYTHONPATH, so a subprocess that changes
directory still imports it.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

WT = Path(__file__).resolve().parents[4]
BIN = Path("/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin")
ENV = {**os.environ, "PATH": f"{BIN}:/usr/bin:/bin", "PYTHONPATH": str(WT / "src"),
       "PYTHONDONTWRITEBYTECODE": "1"}  # fmt: skip

mutations = json.loads(Path(sys.argv[1]).read_text())
log = Path(sys.argv[2])
wanted = set(sys.argv[3:])
for m in mutations:
    if wanted and m["label"] not in wanted:
        continue
    path = WT / m["file"]
    original = path.read_text()
    if original.count(m["old"]) != 1:
        line = f"NOT PLANTED: {m['label']} ({original.count(m['old'])} matches)"
    else:
        try:
            path.write_text(original.replace(m["old"], m["new"]))
            done = subprocess.run(
                [str(BIN / "python"), "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider",
                 *m["tests"]],
                capture_output=True, text=True, check=False, env=ENV, cwd=WT,
            )  # fmt: skip
        finally:
            path.write_text(original)
        tail = (done.stdout.strip().splitlines() or ["(no output)"])[-1]
        line = f"{'CAUGHT' if done.returncode else 'SURVIVED'}: {m['label']} ({tail})"
    print(line, flush=True)
    with log.open("a") as out:
        out.write(line + "\n")
