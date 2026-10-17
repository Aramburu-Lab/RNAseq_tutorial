"""S.8 probe: ctx qc twice on one synthetic root; every CSV compared by bytes and as sets.

PYTHONPATH=src:tests python docs/plans/logs/wo13_e2a/s8_probe.py <scratch dir> [runs]
"""

import hashlib
import shutil
import subprocess
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "tests"))
from test_qc_run import CTX, ENV, plant, with_annotation  # noqa: E402

scratch = Path(sys.argv[1]).resolve()
runs = int(sys.argv[2]) if len(sys.argv) > 2 else 2  # noqa: PLR2004
shutil.rmtree(scratch, ignore_errors=True)
base = scratch / "base"
base.mkdir(parents=True)
plant(base)
spec = with_annotation(base, scratch)
digests: dict[str, set[str]] = {}
frames: dict[str, list[pl.DataFrame]] = {}
for i in range(runs):
    root = scratch / f"run{i}"
    shutil.copytree(base, root)
    done = subprocess.run([str(CTX), "--root", str(root), "qc", "--dataset", "desynth",
                           "--dataset-file", str(spec)], env=ENV, capture_output=True, text=True,
                          check=False)  # fmt: skip
    assert done.returncode == 0, done.stderr[-2000:]
    for csv in sorted((root / "outputs" / "02_qc").glob("*.csv")):
        name = csv.name.split("_desynth_", 1)[1]
        digests.setdefault(name, set()).add(hashlib.sha256(csv.read_bytes()).hexdigest())
        frames.setdefault(name, []).append(pl.read_csv(csv))
for name, seen in sorted(digests.items()):
    first = frames[name][0]
    as_sets = all(f.sort(f.columns).equals(first.sort(first.columns)) for f in frames[name])
    print(f"{name}: {len(seen)} distinct byte strings over {runs} runs; equal as sets {as_sets}")
