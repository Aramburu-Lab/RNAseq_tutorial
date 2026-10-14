"""ctx.clinical's expression tests, Dunn's pairs, grade and the box's two levels on made-up data.

Dunn's test has one definition, ``ctx.proteo.measure.class_test`` (M7); M10 imports it.
"""

import numpy as np
import polars as pl
import pytest
from scipy import stats

from ctx.clinical import ClinicalError, analyses, figures, harmonise
from ctx.proteo import measure


def patients(n: int = 60, seed: int = 1) -> tuple[pl.DataFrame, dict[str, np.ndarray]]:
    """Made-up patients: stage 1 to 4, a rare histology, sex, age; one gene rising with stage."""
    rng = np.random.default_rng(seed)
    stage = np.repeat([1, 2, 3, 4], n // 4)
    histology = np.array(["ductal"] * (n - 25) + ["lobular"] * 20 + ["rare"] * 5)
    table = pl.DataFrame({"stage": stage, "histology": histology,
                          "sex": rng.choice(["female", "male"], n),
                          "age": rng.integers(30, 90, n)})  # fmt: skip
    gene = stage + rng.normal(0, 1.5, n)
    gene[:3] = gene[3]  # ties, which Dunn's z corrects for
    return table, {"G": gene}


def test_plans_follow_each_variable_in_the_declared_order():
    """Catches a variable's tests reordered, or sex tested as ER's positive and negative."""
    assert analyses.plans(("stage", "sex", "age")) == [
        ("stage", "kruskal"),
        ("stage", "spearman"),
        ("sex", "kruskal"),
        ("age", "spearman"),
    ]


def test_dunn_pairs_are_measures_one_definition():
    """The Kruskal-Wallis row and its pairs equal ``measure.class_test`` on the levels with at
    least MIN_GROUP patients ("rare", 5 patients, left out). Mutations caught: Dunn's z without
    the tie correction, BH over every pair of the table in place of each test's, and a level
    under MIN_GROUP kept."""
    table, values = patients()
    tests, pairs = analyses.expression_tests(table, values, [("stage", "kruskal"),
                                                             ("histology", "kruskal")])  # fmt: skip
    gene = values["G"]
    for variable in ("stage", "histology"):
        labels = table[variable].to_numpy()
        kept = [lv for lv in np.unique(labels) if (labels == lv).sum() >= analyses.MIN_GROUP]
        found = measure.class_test({str(lv): gene[labels == lv] for lv in kept})
        row = tests.filter(pl.col("variable") == variable).row(0, named=True)
        assert row["p"] == found.p == stats.kruskal(*(gene[labels == lv] for lv in kept)).pvalue
        assert row["effect"] == found.epsilon2 and row["n"] == sum(found.n.values())
        mine = pairs.filter(pl.col("variable") == variable)
        assert mine.select("level_a", "level_b", "z", "p", "p_adj").rows() == list(found.pairs)
    assert pairs.filter(pl.col("variable") == "histology").height == 1
    assert "rare" not in pairs["level_b"].to_list()


def test_a_whole_number_level_prints_whole_beside_a_missing_value():
    """Stage with one patient missing: Dunn's levels print as 1 to 4, the missing patient in
    none. Mutation caught: the column read through numpy, which makes stage 1 the level 1.0."""
    table, values = patients()
    stage = [None, *table["stage"].to_list()[1:]]
    gapped = table.with_columns(pl.Series("stage", stage, dtype=pl.Int64))
    tests, pairs = analyses.expression_tests(gapped, values, [("stage", "kruskal")])
    assert set(pairs["level_a"]) | set(pairs["level_b"]) == {"1", "2", "3", "4"}
    assert tests["n"].item() == len(stage) - 1


def test_a_kruskal_test_with_one_level_is_refused():
    """Catches a cohort declaring a variable its patients hold one level of: no test exists."""
    table, values = patients()
    one = table.with_columns(pl.lit("only").alias("histology"))
    with pytest.raises(ClinicalError, match="no test"):
        analyses.expression_tests(one, values, [("histology", "kruskal")])


def test_grade_reads_g1_to_g4_and_keeps_gx_as_a_reason():
    """Catches GX read as a grade and a missing token read as GX."""
    assert harmonise.grade("G3") == (3, None)
    assert harmonise.grade("GX") == (None, "GX")
    assert harmonise.grade("[Not Available]") == (None, "not available")


def test_stage_groups_split_two_from_three():
    """Catches stage II put with III to IV, and a patient without stage given a level."""
    table = pl.DataFrame({"stage": [1, 2, 3, 4, None], "er": [None] * 5},
                         schema={"stage": pl.Int64, "er": pl.String})  # fmt: skip
    levels, (early, late) = figures.two_levels(table, "stage_group")
    assert levels.tolist() == [early, early, late, late, None]
    er = pl.DataFrame({"er": ["positive", "negative", None]})
    levels, (negative, positive) = figures.two_levels(er, "er")
    assert levels.tolist() == [positive, negative, None]


def test_competing_causes_come_from_the_os_and_dss_flags():
    """M10.11: a DSS event is a death from the cancer, an OS death with DSS 0 one from other
    causes, an OS death without DSS is left out and counted; 1 - KM of DSS lies at or above the
    cancer's incidence. Mutations caught: the causes swapped, a death of unknown cause kept, and
    a DSS time that differs from OS accepted."""
    rows = pl.DataFrame({
        "OS": [1, 1, 0, 1, 1, 0], "DSS": [1, 0, 0, None, 1, 0],
        "OS_days": [100, 200, 300, 150, 400, 500], "DSS_days": [100, 200, 300, None, 400, 500],
        "OS_months": [1.0, 2.0, 3.0, 1.5, 4.0, 5.0],
    })  # fmt: skip
    found = figures.Competing.of(rows)
    assert found.cause.tolist() == [figures.CANCER, figures.OTHER, 0, figures.CANCER, 0]
    assert found.unknown == 1
    at = found.at(10)
    # By hand: 1/5 at month 1, then 3/5 still event-free times 1/2 at month 4.
    assert at["cancer"] == pytest.approx(100 * (1 / 5 + (3 / 5) * (1 / 2)))
    assert at["other"] == pytest.approx(100 * (4 / 5) * (1 / 4))
    assert at["naive"] >= at["cancer"]
    moved = rows.with_columns(pl.Series("DSS_days", [100, 201, 300, None, 400, 500]))
    with pytest.raises(ValueError, match="DSS time differs"):
        figures.Competing.of(moved)
