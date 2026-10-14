"""The made-up M10 preview draws through ctx.clinical, the built module's one definition.

The preview once carried its own Kaplan-Meier, log-rank, Breslow Cox, PH test, C-index and power
(TODO.md, Priority 2): a second definition that could drift from the page it previews.
"""

import ast
import importlib.util
import sys
from datetime import date
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from ctx.clinical import analyses, cutpoint

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "python" / "preview_clinical.py"
DATASET = "made"
DAY = date(2026, 10, 9)
# What the built module defines and the preview must import, never define again.
RETIRED = {"kaplan_meier", "logrank", "cox", "ph_test", "c_index", "power", "rmst", "cramers_v",
           "km_path", "first_below_half", "risk_table"}  # fmt: skip


@pytest.fixture(scope="module")
def preview():
    sys.path.insert(0, str(SCRIPT.parent))
    spec = importlib.util.spec_from_file_location("preview_clinical", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def built(preview, tmp_path_factory):
    """The preview page on a made-up gathered sheet of 24 samples; returns the root."""
    root = tmp_path_factory.mktemp("preview")
    gather = root / "outputs" / "01_gather"
    gather.mkdir(parents=True)
    rng = np.random.default_rng(3)
    n = 24
    pl.DataFrame(
        {
            "sample": [f"s{i}" for i in range(n)],
            "condition": rng.choice(["a", "b"], n),
            "cell_type": rng.choice(["x", "y"], n),
            "batch": rng.choice(["1", "2", "3"], n),
            "sex": rng.choice(["f", "m"], n),
        },
    ).write_parquet(gather / f"{DAY}_gather_{DATASET}_samplesheet_v1.parquet")
    pl.DataFrame({"gene_id": ["g1"], "gene_name": ["GENE1"]}).write_parquet(
        gather / f"{DAY}_gather_{DATASET}_genes_v1.parquet"
    )
    page = preview.page(root, DATASET, DAY)
    assert page.is_file()
    return root


def written(root: Path, name: str, suffix: str = "parquet") -> Path:
    return root / "outputs" / "12_clinical" / f"{DAY}_clinical_{DATASET}_{name}_v1.{suffix}"


def table(root: Path, name: str) -> pl.DataFrame:
    return pl.read_parquet(written(root, name))


def test_the_preview_defines_no_statistic_of_its_own():
    """Mutation caught: a local Cox, Kaplan-Meier or log-rank put back into the preview."""
    defined = {n.name for n in ast.walk(ast.parse(SCRIPT.read_text()))
               if isinstance(n, ast.FunctionDef)}  # fmt: skip
    assert defined.isdisjoint(RETIRED), defined & RETIRED


def test_the_previews_cox_table_is_ctx_clinicals(built, preview):
    """The preview's Cox rows equal ``analyses.cox_model`` on the cohort it wrote. Mutation
    caught: any other fit of the same data, such as the retired Breslow fit."""
    cohort, cox = table(built, "cohort"), table(built, "cox")
    t, e = cohort["months"].to_numpy(), cohort["event"].to_numpy()
    for row, used in zip(cox.iter_rows(named=True), (["z"], ["z", *preview.COVARIATES]),
                         strict=True):  # fmt: skip
        model = analyses.cox_model(t, e, cohort.select(used).to_numpy().astype(float))
        assert row["hr"] == model.hr and row["ph_p"] == model.ph_p
        assert row["c_index"] == model.c_index and row["events"] == model.events


def test_the_previews_km_and_scan_are_ctx_clinicals(built, preview):
    """Mutation caught: the curves estimated elsewhere, or the scan on another generator."""
    cohort = table(built, "cohort")
    t, e = cohort["months"].to_numpy(), cohort["event"].to_numpy()
    names = list(preview.groups())
    expected = analyses.km_frame(t, e, cohort["group"].to_numpy(), names)
    assert table(built, "km_estimates").equals(expected)
    found = cutpoint.scan(t, e, cohort["log2cpm"].to_numpy())
    scan = pl.read_csv(written(built, "cut_scan", "csv"))
    assert np.allclose(scan["chi2"].to_numpy(), found.statistic, rtol=1e-12, atol=0)
