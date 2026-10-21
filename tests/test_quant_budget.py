"""M4.11: `ctx quant`'s laptop path on the 8 real teaching runs keeps Track A.

Marked ``budget`` and deselected by default; run with ``-m budget -s``. The runs are M3's
trimmed reads of tier1_teach in the reader's Step-1 root (``CTX_TEACH_ROOT``), linked in with
the gather sheet, the annotation and chr22, so a cold run builds the reference, Salmon's
index with chr22 as decoys and every quantification.
"""

import os
import subprocess
from pathlib import Path

import pytest
from test_budget import CTX, GENOME, TEACH_ROOT, _measured

from ctx.report.builder import unpacked
from ctx.runtime import BUDGETS

LAPTOP_INPUTS = ("outputs/01_gather", "outputs/03_process", "resources", GENOME)


@pytest.mark.budget
def test_quant_laptop_path_keeps_track_a(tmp_path):
    """A cold laptop path on 8 runs keeps 15 minutes and 6 GB; its page stays under 15 MB.

    Memory is bounded as M3's Track A test bounds it: the command plus threads times its
    largest child, Salmon. Measured 2026-10-09 under a load average near 45: wall 253.8 s,
    self 353,206,272 B, largest child 1,269,972,992 B, page 5,850,923 B.
    """
    missing = [p for p in LAPTOP_INPUTS if not (TEACH_ROOT / p).exists()]
    if missing:
        pytest.skip(f"the teaching runs' M3 outputs are not on this machine: {missing}")
    for part in LAPTOP_INPUTS:
        (tmp_path / part).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / part).symlink_to(TEACH_ROOT / part)
    env = {k: v for k, v in os.environ.items() if not k.startswith(("SLURM_", "CTX_"))}
    repo = Path(__file__).resolve().parents[1]
    done = subprocess.run(
        [str(CTX), "--root", str(tmp_path), "quant", "--dataset", "tier1_teach",
         "--dataset-file", str(repo / "conf" / "datasets" / "tier1_teach.toml"),
         "--genome", str(tmp_path / GENOME)],
        cwd=tmp_path, env=env, capture_output=True, text=True, check=False,
    )  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    (page,) = (tmp_path / "outputs" / "04_quant").glob("*_tier1_teach_overview_v1.html")
    wall, self_rss, children = _measured(unpacked(page.read_text()))
    threads = min(4, os.cpu_count() or 1)
    bound = self_rss + threads * children
    print(f"\nquant laptop: wall {wall:.1f} s, self {self_rss:,} B, largest child {children:,} B, "
          f"bound {bound:,} B ({threads} threads), page {page.stat().st_size:,} B")  # fmt: skip
    assert wall < BUDGETS["A"].minutes * 60
    assert bound < BUDGETS["A"].ram_gb * 10**9
    assert page.stat().st_size < 15 * 10**6
