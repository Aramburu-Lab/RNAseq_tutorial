"""Package wo13-l: run one command with main's environment and this worktree's src.

Usage: <main env python> docs/plans/logs/wo13_l/env.py <command ...>
"""

import os
import subprocess
import sys
from pathlib import Path

BIN = "/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin"
WT = Path(__file__).absolute().parents[4]
env = dict(os.environ)
env["PATH"] = BIN + os.pathsep + env.get("PATH", "")
env["PYTHONPATH"] = env.get("SRC") or str(WT / "src")
env["PYTHONDONTWRITEBYTECODE"] = "1"
env["CTX_CHROMIUM"] = "/Applications/Chromium.app/Contents/MacOS/Chromium"
sys.exit(subprocess.run(sys.argv[1:], cwd=WT, env=env, check=False).returncode)
