"""Package P: plant each mutation in place, run its tests, restore, and log CAUGHT or SURVIVED.

Run: python mutate.py <mutations.json> <log>
Each mutation: {"label", "file", "old", "new", "tests": [pytest args]}; "old" must occur once.
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
    if original.count(m["old"]) != 1:
        line = f"NOT PLANTED: {m['label']} ({original.count(m['old'])} matches)"
    else:
        try:
            path.write_text(original.replace(m["old"], m["new"]))
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
