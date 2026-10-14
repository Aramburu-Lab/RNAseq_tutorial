"""M7 at real size (plan P7.11): ctx proteo on CPTAC breast with F4 and on Step-1 against
Rieckmann, from this Mac's roots, against the Track B budget, under 10 min and 4 GB each
(CONVENTIONS.md, Stack), and the 15 MB page cap. Inputs are linked, never written.

CPTAC's PAM50 and ER calls are joined onto gather's sheet here, by ``pair_samples``, the way the
declared annotation join is to write them; that join is not built (it needs ctx.gather and
ctx.dataset). Opt-in: -m "realdata and budget"; minutes."""

import os
import re
import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

from ctx.gather import gdc
from ctx.proteo import tables
from ctx.proteo.run import pair_samples
from ctx.report.builder import MAX_BYTES
from ctx.runtime import BUDGETS

ROOT = Path(__file__).resolve().parents[1]
HOME = Path(os.environ.get("CTX_GDC_HOME", Path.home()))
FREEZE = "resources/datasets/cptac_brca/S060_Breast_Cancer_data_freeze_GCTfiles_v5.4-public.tar.gz"
MEMBER = "prosp-brca-v5.4-public-proteome-ratio-norm-NArm.gct"
# The freeze's annotation columns and the sheet columns cptac_brca.toml's contrasts name.
JOINED = {
    "pam50": ("PAM50", {}),
    "er_status": ("ER.Updated.Clinical.Status", {"positive": "ER+", "negative": "ER-"}),
}


def measured(page: Path) -> tuple[float, int]:
    """Wall seconds and this process's peak RSS in bytes, as the page's Reproduce table prints."""
    html = page.read_text()
    wall = re.search(r'<td class="wall">([\d.]+) s</td>', html)
    rss = re.search(r'<td class="rss-self">[\d.]+ MiB \((\d+) bytes\)</td>', html)
    assert wall and rss
    return float(wall.group(1)), int(rss.group(1))


def ctx(root: Path, *args: str) -> None:
    """Run ctx in this interpreter, so the subprocess imports the checkout under test."""
    done = subprocess.run([sys.executable, "-m", "ctx", "--root", str(root), *args],
                          capture_output=True, text=True, check=False)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]


def link(source: Path, root: Path, relative: str) -> None:
    """Link one input into the scratch root at its place under a ctx root."""
    target = root / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    target.symlink_to(source / relative)


def joined_sheet(home: Path, root: Path) -> None:
    """Gather's CPTAC sheet with pam50 and er_status from the freeze, written into ``root``."""
    gather = home / "outputs" / "01_gather"
    (path,) = sorted(gather.glob("*_gather_cptac_brca_samplesheet_v1.parquet"))[-1:]
    sheet = pl.read_parquet(path)
    layer, annotation = tables.read_gct(home / FREEZE, MEMBER)
    pairs, _ = pair_samples(layer.samples, gdc.chosen(sheet, "Primary Tumor"), "X")
    added = {}
    for column, (source, spelling) in JOINED.items():
        found = {s: spelling.get(v, v) for j, s in pairs if (v := annotation[source][j])}
        added[column] = pl.col("sample").replace_strict(found, default=None,
                                                        return_dtype=pl.String)  # fmt: skip
    (root / "outputs" / "01_gather").mkdir(parents=True, exist_ok=True)
    sheet.with_columns(**added).write_parquet(root / "outputs" / "01_gather" / path.name)
    for what in ("counts", "genes"):
        link(home, root, str(path.relative_to(home)).replace("samplesheet", what))


def within(root: Path, module: str, slug: str) -> None:
    """The module's page within Track B and the page cap."""
    budget = BUDGETS["B"]
    (page,) = (root / "outputs" / module).glob(f"*_{slug}_overview_v1.html")
    wall, rss = measured(page)
    assert wall < budget.minutes * 60 and rss < budget.ram_gb * 2**30, (module, wall, rss)
    assert page.stat().st_size <= MAX_BYTES, page.stat().st_size


@pytest.mark.realdata
@pytest.mark.budget
def test_cptac_brca_with_f4_fits_track_b(tmp_path):
    """Catches a run over the Track B budget or a page over the cap at real size, with both F4
    contrasts drawn; and K1 to K3 lost. Fixture sizes hide both (lessons.md, 2026-09-26)."""
    home = HOME / "ctx_breast"
    if not (home / FREEZE).is_file():
        pytest.skip(f"{home / FREEZE} is not on this machine")
    root = tmp_path / "root"
    joined_sheet(home, root)
    for name in (FREEZE, "resources/msigdb/c5.go.cc.v2026.1.Hs.symbols.gmt",
                 "resources/msigdb/c2.cp.reactome.v2026.1.Hs.symbols.gmt"):  # fmt: skip
        link(home, root, name)
    spec = ROOT / "conf" / "datasets" / "cptac_brca.toml"
    with_de = tmp_path / "cptac_brca.toml"
    with_de.write_text(spec.read_text() + '\n[de]\nfeatured = "basal_vs_luma"\n')
    ctx(root, "de", "--dataset", "cptac_brca", "--dataset-file", str(with_de))
    ctx(root, "proteo", "--dataset", "cptac_brca", "--dataset-file", str(spec))
    for module in ("05_de", "07_proteo"):
        within(root, module, "cptac_brca")
    out = root / "outputs" / "07_proteo"
    (verdicts,) = out.glob("*_sanity_v1.parquet")
    rules = pl.read_parquet(verdicts).filter(pl.col("rule").is_in(["K1", "K2", "K3"]))
    assert set(rules["rule_verdict"]) == {"PASS"}, rules.rows()
    assert len(list(out.glob("*_contrast_*_v1_white.pdf"))) == 2


@pytest.mark.realdata
@pytest.mark.budget
def test_tier1_against_rieckmann_fits_track_b(tmp_path):
    """Catches the Step-1 page over the budget or the cap at real size."""
    home = HOME / "ctx_tier1"
    results = sorted(
        (home / "outputs" / "05_de").glob("*_de_tier1_ndn_vs_mono_steady_results_v1.parquet")
    )
    if not results:
        pytest.skip(f"{home} has no ndn_vs_mono_steady results")
    root = tmp_path / "root"
    for name in (str(results[-1].relative_to(home)),
                 "resources/datasets/tier1/41590_2017_BFni3693_MOESM13_ESM.xlsx",
                 "resources/msigdb/c5.go.cc.v2026.1.Hs.symbols.gmt"):  # fmt: skip
        link(home, root, name)
    ctx(root, "proteo", "--dataset", "tier1", "--dataset-file",
        str(ROOT / "conf" / "datasets" / "tier1.toml"))  # fmt: skip
    within(root, "07_proteo", "tier1")
