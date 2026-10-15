"""K.9's real check: a kit from ~/ctx_kidney alone, its LEFT_OUT.json commands tabulated.

Run through K_kidney.sh, which builds the kit under the worktree's .scratch_k folder. Prints,
per left-out file, the --dataset-file the reader is shown and whether README.md names the
builder's root and any dataset file of the builder's own.
"""

import json
import sys
from pathlib import Path

kit = Path(sys.argv[1])
entries = json.loads((kit / "LEFT_OUT.json").read_text(encoding="utf-8"))["files"]
readme = (kit / "README.md").read_text(encoding="utf-8")
for e in entries:
    words = e["rebuild"][-1].split()
    shown = words[words.index("--dataset-file") + 1] if "--dataset-file" in words else None
    print(f"{e['path']}\t{shown}\troot={e['builder_root']}\town={e['builder_dataset_file']}")
mapped = sum("--dataset-file conf/datasets/" in e["rebuild"][-1] for e in entries)
roots = readme.count("is the builder's ctx root")
print(f"{len(entries)} left out; {mapped} name conf/datasets/; README.md names the root {roots}x")
