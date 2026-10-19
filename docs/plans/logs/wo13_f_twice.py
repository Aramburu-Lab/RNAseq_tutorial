"""wo13-f: draw each spec'd figure of a scratch root twice with ctx figure; compare the files.

Usage: wo13_f_py.sh python docs/plans/logs/wo13_f_twice.py <root> <out>
Prints, per figure and file kind, whether the two renders are byte-equal, and what differs.
"""

import shutil
import subprocess
import sys
from pathlib import Path

from ctx.dataset import load
from ctx.figures import record
from ctx.figures.params import params_file

W = Path(__file__).resolve().parents[3]
root, out = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
load(W / "tests" / "fixtures" / "datasets" / "de_synthetic.toml")
for folder in ("05_de", "06_function"):
    for csv in sorted((root / "outputs" / folder).glob("*_v1.csv")):
        spec = record.describe(csv.name[:-4], csv.parent, csv)
        if spec is None:
            continue
        drawn = []
        for n in (1, 2):
            place = out / f"{spec.figure_id}_{n}"
            shutil.rmtree(place, ignore_errors=True)
            place.mkdir(parents=True)
            for item in spec.inputs:
                top = {"root": root, "repo": W}.get(item.base)
                src = top / item.source if top else Path(item.source).expanduser()
                shutil.copyfile(src, place / item.file)
            (place / "params.toml").write_text(params_file(spec))
            subprocess.run(["ctx", "figure", spec.figure_id, "--params", "params.toml"],
                           cwd=place, check=True, capture_output=True)  # fmt: skip
            drawn.append(place)
        for ext in ("_white.pdf", "_black.pdf", "_white.svg", "_black.svg", ".csv"):
            a, b = (next(p.glob(f"*{ext}")).read_bytes() for p in drawn)
            first = next((i for i, (x, y) in enumerate(zip(a, b, strict=False)) if x != y), None)
            where = "" if a == b else a[max(0, (first or 0) - 40) : (first or 0) + 40]
            print(f"{spec.figure_id} {ext} {'EQUAL' if a == b else 'DIFFER'} {len(a)} {where!r}")
