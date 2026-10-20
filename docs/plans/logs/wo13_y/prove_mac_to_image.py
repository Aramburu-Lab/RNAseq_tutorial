"""wo13-y item 6: a figure drawn on this Mac, redrawn in the figures image, compared byte for byte.

    python prove_mac_to_image.py --image ctx-figures:y --root <Mac-drawn root> --work <empty>

containers/prove_figures.py's own steps and comparison, with one change: the Docker tab clones
the course repository from GitHub at the drawing commit, which this unpushed branch's commit
is not on, so the clone line here reads this worktree's repository instead. The clone serves
the one repository input, the dataset file, whose bytes the tab still checks by SHA-256.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "containers"))
from prove_figures import compare, volcano_tab  # noqa: E402

GITHUB = "https://github.com/Aramburu-Lab/RNAseq_tutorial.git"

parser = argparse.ArgumentParser()
parser.add_argument("--image", required=True)
parser.add_argument("--root", required=True, type=Path)
parser.add_argument("--work", required=True, type=Path)
args = parser.parse_args()
script, stem = volcano_tab(args.root)
# A root drawn outside a git checkout (the control, from an export of main's src) names no
# commit, so its tab copies the dataset file from the repository folder and clones nothing.
assert script.count(f"git clone {GITHUB} ") <= 1, "the tab clones twice"
reader = args.work / "reader"
reader.mkdir(parents=True)
(reader / "reproduce.sh").write_text(script.replace(f"git clone {GITHUB} ", f"git clone {REPO} "))
env = {**os.environ, "CTX_FIGURES_IMAGE": args.image}
done = subprocess.run(["bash", "reproduce.sh"], cwd=reader, env=env, check=False)
if done.returncode:
    sys.exit(done.returncode)
(redrawn,) = reader.glob("repro_*")
if not compare(args.root / "outputs" / "05_de", redrawn, stem):
    sys.exit(1)
print("every compared file is byte-equal")
