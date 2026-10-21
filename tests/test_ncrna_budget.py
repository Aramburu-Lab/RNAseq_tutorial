"""M8's TCGA half at real size (plan step 14): ``ctx ncrna`` on the breast cohort from this Mac's
folders, against Track B, under 10 min and 4 GB (CONVENTIONS.md, Stack), and the 15 MB cap.

The root is built as a reader builds it: gather with the miRNA folders, setup registering the
target files and GDC's GENCODE v36 by hash, then M5's paired fit. M5's TCGA page is over the
cap (16,164,875 B on 2026-10-08), so ``ctx de`` exits non-zero after writing its tables; this
test fits M5's contrast itself, with M5's own ``fit``, and writes the results table ``ctx de``
writes. Opt-in: -m "realdata and budget"; about five minutes.
"""

import os
import re
import subprocess
import sys
from datetime import date
from pathlib import Path

import polars as pl
import pytest

from ctx.dataset import load
from ctx.de.deseq import RESULT_COLUMNS, fit
from ctx.gather import gdc
from ctx.io.naming import stem
from ctx.io.readme_companion import write_table
from ctx.report.builder import unpacked
from ctx.report.cap import MAX_BYTES
from ctx.runtime import BUDGETS

REPO = Path(__file__).resolve().parents[1]
CTX = Path(sys.executable).parent / "ctx"
NAME = "tcga_brca"
SPEC = REPO / "conf" / "datasets" / f"{NAME}.toml"
HOME = Path(os.environ.get("CTX_GDC_HOME", Path.home()))
LINKS = {f"data/raw/{NAME}{k}": HOME / "ctx_breast" / "data/raw" / f"{NAME}{k}"
         for k in ("", "_mirna", "_isoform")}  # fmt: skip
LINKS |= {
    "data/raw/mirna_targets": HOME / "ctx_tier1" / "data/raw/mirna_targets",
    "data/raw/gdc_gencode_v36": HOME / "ctx_breast" / "data/raw/gdc_gencode_v36",
}


def _run(root: Path, *args: str) -> None:
    """One ctx command on the root; its failure shows its last lines."""
    done = subprocess.run([str(CTX), "--root", str(root), *args, "--dataset-file", str(SPEC)],
                          capture_output=True, text=True, check=False)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]


def _m5(root: Path) -> None:
    """M5's tumour-normal fit, written as ``ctx de`` writes its results table."""
    spec = load(SPEC)
    out = root / "outputs" / "01_gather"
    sheet, counts, genes = (pl.read_parquet(next(out.glob(f"*_{NAME}_{w}_v1.parquet")))
                            for w in ("samplesheet", "counts", "genes"))  # fmt: skip
    contrast = spec.contrasts[0]
    result = fit(counts, gdc.chosen(sheet), contrast, spec.order, threads=4, shrink=False)
    frame = result.results.join(genes.select("gene_id", "gene_name"), on="gene_id", how="left")
    write_table(frame, root / "outputs" / "05_de",
                stem("de", f"{NAME}_{contrast.name}_results", 1, day=date.today()),
                columns={**RESULT_COLUMNS, "gene_name": ("none", "Gene symbol.")},
                data_releases={})  # fmt: skip


@pytest.mark.realdata
@pytest.mark.budget
def test_tcga_brca_ncrna_keeps_track_b(tmp_path):
    """Catches a run over Track B and a page over the cap at real size: 1,188 miRNA libraries
    and 1,213 STAR files, the sizes fixtures hide (lessons.md, 2026-09-26)."""
    missing = [str(p) for p in LINKS.values() if not p.is_dir()]
    if missing:
        pytest.skip(f"not on this machine: {missing}")
    root = tmp_path / "root"
    for dest, source in LINKS.items():
        (root / dest).parent.mkdir(parents=True, exist_ok=True)
        (root / dest).symlink_to(source)
    raw = root / "data" / "raw"
    _run(
        root,
        "gather",
        "--dataset",
        NAME,
        "--gdc-dir",
        str(raw / NAME),
        "--mirna-dir",
        str(raw / f"{NAME}_mirna"),
        "--isoform-dir",
        str(raw / f"{NAME}_isoform"),
    )
    _run(root, "setup")
    _m5(root)
    _run(root, "ncrna", "--dataset", NAME)
    (page,) = (root / "outputs" / "08_ncrna").glob(f"*_{NAME}_overview_v1.html")
    html = unpacked(page.read_text())
    wall = float(re.search(r'<td class="wall">([\d.]+) s</td>', html).group(1))  # type: ignore[union-attr]
    rss = int(re.search(r'<td class="rss-self">[\d.]+ MiB \((\d+) bytes\)</td>', html).group(1))  # type: ignore[union-attr]
    print(f"\nncrna {NAME}: wall {wall:.1f} s, self {rss:,} B, page {page.stat().st_size:,} B")
    assert wall < BUDGETS["B"].minutes * 60 and rss < BUDGETS["B"].ram_gb * 2**30
    assert page.stat().st_size <= MAX_BYTES
