"""wo13-f: a scratch root inside the worktree, de and function on tests/synth.py's synthetic data.

Usage: wo13_f_py.sh python docs/plans/logs/wo13_f_scratch.py <root>
Nothing real is read; the root is for looking at pages, never committed.
"""

import subprocess
import sys
from datetime import date
from pathlib import Path

W = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(W / "tests"))
from synth import DE_GROUPS, design_counts, write_gather_outputs, write_gene_sets  # noqa: E402

from ctx.dataset import load  # noqa: E402

root = Path(sys.argv[1]).resolve()
root.mkdir(parents=True, exist_ok=True)
spec = W / "tests" / "fixtures" / "datasets" / "de_synthetic.toml"
write_gather_outputs(root, "desynth", design_counts(DE_GROUPS, batch_lfc=2.0), date.today())
write_gene_sets(root / load(spec).resource("sets").dest)
for module in ("de", "function"):
    args = ["ctx", "--root", str(root), module, "--dataset", "desynth", "--dataset-file", str(spec)]
    done = subprocess.run(args, cwd=root, check=False, capture_output=True, text=True)
    print(module, done.returncode, done.stderr[-2000:])
