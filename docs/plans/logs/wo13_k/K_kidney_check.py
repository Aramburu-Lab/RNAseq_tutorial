"""K.9's real check on ~/ctx_kidney, read only: the LEFT_OUT.json entries build_kit writes.

The whole kit is refused on this root today (D-072: outputs/11_genecard and 12_clinical hold
no index.html until M9 and M10 are rebuilt there), so this runs build_kit's own select() and
carry_matrices() on the root, which compute every left-out entry before that refusal, and
writes the GDC manifests and LEFT_OUT.json under the given scratch folder. Prints, per entry,
the --dataset-file the reader is shown, the builder's root and any dataset file of the
builder's own, and the README lines build_kit would write under it.
"""

import importlib.util
import json
import sys
from pathlib import Path

W = Path(__file__).resolve().parents[4]
spec = importlib.util.spec_from_file_location("build_kit", W / "scripts/python/build_kit.py")
assert spec is not None and spec.loader is not None
kit_script = importlib.util.module_from_spec(spec)
spec.loader.exec_module(kit_script)

root, out = Path(sys.argv[1]).expanduser(), Path(sys.argv[2])
slugs = frozenset(kit_script.slugs_in(root))
_, matrices, problems = kit_script.select(root, (), slugs, None, kit_script.codebase_datasets())
print("refusals on this root:", *problems, sep="\n  ")
out.mkdir(parents=True)
entries = kit_script.carry_matrices(out, "kidney", (root / "outputs").resolve(), matrices)
(out / "LEFT_OUT.json").write_text(json.dumps({"files": entries}, indent=1) + "\n")
for e in entries:
    words = e["rebuild"][-1].split()
    shown = words[words.index("--dataset-file") + 1] if "--dataset-file" in words else None
    print(f"{Path(e['path']).name}\t{shown}\troot={e['builder_root']}"
          f"\town={e['builder_dataset_file']}")  # fmt: skip
mapped = sum("--dataset-file conf/datasets/" in e["rebuild"][-1] for e in entries)
print(f"{len(entries)} left out; {mapped} name conf/datasets/")
print("README.md under the first:", kit_script.left_out_text(entries[0]), sep="\n")
