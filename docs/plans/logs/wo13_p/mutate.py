"""Package P: plant each mutation in place, run its tests, restore, and log CAUGHT or SURVIVED.

Run: python mutate.py <mutations.json> <log>
Each mutation: {"label", "file", "old" and "new" or "edits": [[old, new]], "tests"}; each old
must occur once in the file.
The source is restored whatever happens, and the run writes no bytecode (lessons.md).
"""

import json
import os
import subprocess
import sys
from pathlib import Path

WT = Path(__file__).resolve().parents[4]
PYTEST = WT / "docs" / "plans" / "logs" / "wo13_p" / "pytest.sh"

mutations = json.loads(Path(sys.argv[1]).read_text())
log = Path(sys.argv[2])
for m in mutations:
    path = WT / m["file"]
    original = path.read_text()
    # One edit, or several in the one file as [old, new] pairs under "edits".
    edits = m.get("edits") or [[m["old"], m["new"]]]
    if any(original.count(old) != 1 for old, _ in edits):
        line = f"NOT PLANTED: {m['label']} ({[original.count(o) for o, _ in edits]} matches)"
    else:
        mutated = original
        for old, new in edits:
            mutated = mutated.replace(old, new)
        try:
            path.write_text(mutated)
            done = subprocess.run(["bash", str(PYTEST), "-q", "-x", *m["tests"]],
                                  capture_output=True, text=True, check=False,
                                  env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})  # fmt: skip
        finally:
            path.write_text(original)
        tail = (done.stdout.strip().splitlines() or ["(no output)"])[-1]
        line = f"{'CAUGHT' if done.returncode else 'SURVIVED'}: {m['label']} ({tail})"
    print(line, flush=True)
    with log.open("a") as out:
        out.write(line + "\n")
