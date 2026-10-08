"""M10 at real size (plan M10.9): gather and clinical on the breast cohort from this Mac's
folder, against the Track B budget, under 10 min and 4 GB each (CONVENTIONS.md, Stack), and the
15 MB page cap.
Opt-in: -m "realdata and budget"; minutes."""

import importlib.util
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

from ctx.report.builder import MAX_BYTES
from ctx.runtime import BUDGETS

ROOT = Path(__file__).resolve().parents[1]
CTX = Path(sys.executable).parent / "ctx"
SPEC = ROOT / "conf" / "datasets" / "tcga_brca.toml"
NAME = SPEC.stem
XLSX = importlib.util.find_spec("openpyxl") is not None


def measured(page: Path) -> tuple[float, int]:
    """Wall seconds and this process's peak RSS in bytes, as the page's Reproduce table prints."""
    html = page.read_text()
    wall = re.search(r'<td class="wall">([\d.]+) s</td>', html)
    rss = re.search(r'<td class="rss-self">[\d.]+ MiB \((\d+) bytes\)</td>', html)
    assert wall and rss
    return float(wall.group(1)), int(rss.group(1))


@pytest.mark.realdata
@pytest.mark.budget
@pytest.mark.skipif(not XLSX, reason="openpyxl ==3.1.5 is pinned in pixi.toml; run pixi install")
def test_tcga_brca_gather_and_clinical_fit_track_b(tmp_path):
    """Catches a run over the Track B budget, and a page over the cap, at real size: the
    lesson that fixture sizes hide both (lessons.md, 2026-09-26)."""
    folder = Path(os.environ.get("CTX_GDC_HOME", Path.home())) / "ctx_breast" / "data/raw" / NAME
    if not folder.is_dir():
        pytest.skip(f"{folder} is not on this machine")
    budget = BUDGETS["B"]
    root = tmp_path / "root"
    for args in (("gather", "--gdc-dir", str(folder)), ("clinical",)):
        done = subprocess.run([str(CTX), "--root", str(root), *args, "--dataset", NAME,
                               "--dataset-file", str(SPEC)], capture_output=True, text=True,
                              check=False)  # fmt: skip
        assert done.returncode == 0, done.stderr[-3000:]
    for module in ("01_gather", "12_clinical"):
        (page,) = (root / "outputs" / module).glob(f"*_{NAME}_overview_v1.html")
        wall, rss = measured(page)
        assert wall < budget.minutes * 60 and rss < budget.ram_gb * 2**30, (module, wall, rss)
        assert page.stat().st_size <= MAX_BYTES
