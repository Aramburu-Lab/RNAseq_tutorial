"""``ctx proteo`` end to end on made-up roots (tests/proteo_root.py), both designs, and the
rules' registration order. Each test names the mutation it catches.
"""

import re
import subprocess
from pathlib import Path

import polars as pl
import pytest
from proteo_root import TUMOURS, across, cross

from ctx.report.builder import MAX_BYTES

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def across_root(tmp_path_factory):
    """The across-samples root, ctx proteo run."""
    return across(tmp_path_factory.mktemp("m7_across"))


@pytest.fixture(scope="module")
def cross_root(tmp_path_factory):
    """The cross-dataset root, ctx proteo run."""
    return cross(tmp_path_factory.mktemp("m7_cross"))


def output(root: Path, what: str) -> Path:
    """The one proteo output of this kind."""
    (path,) = (root / "outputs" / "07_proteo").glob(f"*_proteo_*_{what}")
    return path


@pytest.mark.parametrize("which", ["across_root", "cross_root"])
def test_page_offline_and_every_table_described(which, request):
    """Catches a page over the cap or fetching a script, and a table without its README."""
    root = request.getfixturevalue(which)
    html = output(root, "overview_v1.html").read_text()
    assert len(html.encode()) <= MAX_BYTES
    assert not re.search(r"<script[^>]+\bsrc\s*=", html, re.IGNORECASE)
    for table in (root / "outputs" / "07_proteo").glob("*.csv"):
        assert table.with_name(table.stem + "_README.md").is_file(), table.name


@pytest.mark.parametrize("which", ["across_root", "cross_root"])
def test_missingness_first_and_rules_pass(which, request):
    """Catches the class or fold figure drawn before the missingness figures (K3 fails), and
    the planted class shift or fold agreement lost (K2, P1 or P2 fails)."""
    root = request.getfixturevalue(which)
    verdicts = pl.read_parquet(output(root, "sanity_v1.parquet"))
    assert set(verdicts["rule_verdict"]) == {"PASS"}, verdicts.select("rule", "left").rows()
    html = output(root, "overview_v1.html").read_text()
    stems = re.findall(r"_proteo_fake_\w+?_(\w+?)_v1\.csv", html)
    assert stems[:2] == ["quantified", "coverage"], stems


def test_across_missing_never_zero_and_one_row_per_gene(across_root):
    """Catches a missing protein read as 0 (G10 would count every sample), the extra G1 row
    kept as a second gene, and samples joined by position instead of case."""
    genes = pl.read_parquet(output(across_root, "genes_v1.parquet"))
    assert genes["gene"].n_unique() == genes.height
    g10 = genes.filter(pl.col("gene") == "G10").row(0, named=True)
    assert g10["n_quantified"] == TUMOURS // 2 and not g10["tested"]
    g9 = genes.filter(pl.col("gene") == "G9").row(0, named=True)
    assert g9["rho"] > 0.95
    metrics = pl.read_parquet(output(across_root, "metrics_v1.parquet"))
    assert metrics["metric"].n_unique() == metrics.height


def test_cross_untested_below_the_share(cross_root):
    """Catches the 70 % filter skipped: G30 has 5 of 8 values."""
    genes = pl.read_parquet(output(cross_root, "genes_v1.parquet"))
    g30 = genes.filter(pl.col("gene") == "G30").row(0, named=True)
    assert g30["n_samples"] == 5 and not g30["tested"]
    html = output(cross_root, "overview_v1.html").read_text()
    assert "Cross-study comparison: RNA from fake_cells" in html


def git(*args: str) -> str:
    """git in this checkout."""
    return subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True,
                          check=True).stdout.strip()  # fmt: skip


@pytest.mark.parametrize("name", ["tier1", "cptac_brca"])
def test_rules_registered_before_the_protein_reader(name):
    """The registration commit of [proteo.sanity] precedes the commit that added the reader.

    Catches a rule edited after protein values could be read without registering again."""
    if git("rev-parse", "--is-shallow-repository") == "true":
        pytest.skip("shallow clone: the history is absent")
    path = f"conf/datasets/{name}.toml"
    registered = git("log", "--format=%H", "-G", r"^registered_at = \"2026-10", "--", path).split()
    reader = git("log", "--reverse", "--diff-filter=A", "--format=%H", "--",
                 "src/ctx/proteo/tables.py").split()  # fmt: skip
    assert registered and reader
    assert registered[0] != reader[0]
    ancestor = subprocess.run(["git", "-C", str(ROOT), "merge-base", "--is-ancestor",
                               registered[0], reader[0]], check=False)  # fmt: skip
    assert ancestor.returncode == 0
