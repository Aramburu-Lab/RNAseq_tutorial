"""scripts/python/pick_teaching_subset.py applies its written rule, on a synthetic sheet."""

import importlib.util
from pathlib import Path

import polars as pl
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "python" / "pick_teaching_subset.py"
CELLS = ["neut", "mono"]


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("pick_teaching_subset", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run(donor, condition, length, cell, **given):
    fields = {"isolation": "kit", "sex": "F", "read_count": 1_000_000, **given}
    return {"donor": donor, "condition": condition, "read_length": float(length),
            "cell_type": cell, "run_accession": f"r_{donor}_{cell}_{condition}",
            "fastq_bytes": "10", **fields}  # fmt: skip


def sheet(rows):
    return pl.DataFrame(rows)


def both(donor, condition, length, **kw):
    return [run(donor, condition, length, c, **kw) for c in CELLS]


def test_one_donor_per_stratum_each_cell_type(script):
    rows = both("d1", "a", 76) + both("d2", "b", 76) + both("d3", "a", 101) + both("d4", "c", 101)
    picked = script.pick(sheet(rows), [("a", 76), ("b", 76), ("a", 101), ("c", 101)], CELLS)
    assert picked.height == 8
    assert picked.group_by("stratum").agg(pl.col("donor").n_unique())["donor"].to_list() == [1] * 4
    assert set(picked["donor"]) == {"d1", "d2", "d3", "d4"}


def test_a_donor_is_used_once(script):
    # d1 qualifies for both strata; the second must take d2.
    rows = both("d1", "a", 76) + both("d1", "b", 76, sex="M") + both("d2", "b", 76)
    picked = script.pick(sheet(rows), [("a", 76), ("b", 76)], CELLS)
    assert picked.filter(pl.col("stratum") == "b:76")["donor"].unique().to_list() == ["d2"]


def test_isolation_is_shared_across_strata(script):
    # d2 is closer in depth but its monocytes were isolated with beads, which stratum a lacks.
    rows = both("d1", "a", 76) + both("d3", "b", 76, read_count=900_000)
    rows += [run("d2", "b", 76, "neut"), run("d2", "b", 76, "mono", isolation="beads")]
    picked = script.pick(sheet(rows), [("a", 76), ("b", 76)], CELLS)
    assert set(picked["isolation"]) == {"kit"}
    assert "d2" not in set(picked["donor"])


def test_empty_sex_last_then_the_rarer_sex(script):
    rows = both("d1", "a", 76, sex="F")
    rows += both("d2", "b", 76, sex="", read_count=1_000_000) + both(
        "d3", "b", 76, sex="F", read_count=5
    )
    rows += both("d4", "b", 76, sex="M", read_count=7)
    picked = script.pick(sheet(rows), [("a", 76), ("b", 76)], CELLS)
    # d2 has the median depth but no sex; F is taken, so the rarer M wins over F.
    assert picked.filter(pl.col("stratum") == "b:76")["donor"].unique().to_list() == ["d4"]


def test_an_impossible_stratum_is_named(script):
    rows = both("d1", "a", 76) + [run("d2", "b", 76, "neut")]
    with pytest.raises(script.PickError, match="b:76"):
        script.pick(sheet(rows), [("a", 76), ("b", 76)], CELLS)


def test_deterministic(script):
    rows = both("d1", "a", 76) + both("d2", "a", 76) + both("d3", "b", 76)
    strata = [("a", 76), ("b", 76)]
    assert script.pick(sheet(rows), strata, CELLS).equals(script.pick(sheet(rows), strata, CELLS))
