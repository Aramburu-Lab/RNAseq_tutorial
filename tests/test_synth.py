"""The synthetic matrix has the declared design and recovers its planted truth."""

import numpy as np
import polars as pl
import pytest
from synth import counts

PLANTED = 2.0


@pytest.fixture(scope="module")
def synthetic():
    return counts()


def test_design_nested(synthetic):
    sheet = synthetic.samples
    assert sheet.height == 24
    assert sheet["sample"].n_unique() == 24
    per_donor = sheet.group_by("donor").agg(
        pl.col("condition").n_unique().alias("conditions"),
        pl.col("cell_type").n_unique().alias("cell_types"),
    )
    assert per_donor.height == 12
    assert (per_donor["conditions"] == 1).all()
    assert (per_donor["cell_types"] == 2).all()
    matrix = synthetic.counts.drop("gene_id")
    assert matrix.columns == sheet["sample"].to_list()
    assert set(matrix.dtypes) == {pl.Int32}
    assert matrix.to_numpy().min() >= 0


def test_seed_determinism():
    assert counts(42).counts.equals(counts(42).counts)
    assert not counts(42).counts.equals(counts(43).counts)


def _log2_ratio(cpm, num, den):
    return np.log2(cpm[:, num].mean(axis=1) + 0.5) - np.log2(cpm[:, den].mean(axis=1) + 0.5)


def test_planted_truth(synthetic):
    matrix = synthetic.counts.drop("gene_id").to_numpy().astype(float)
    cpm = matrix / matrix.sum(axis=0) * 1e6
    sheet, truth = synthetic.samples, synthetic.truth
    is_b = (sheet["cell_type"] == "type_b").to_numpy()
    is_case = (sheet["condition"] == "case").to_numpy()
    contrasts = {
        "lfc_cell_type": _log2_ratio(cpm, is_b & ~is_case, ~is_b & ~is_case),
        "lfc_condition_type_a": _log2_ratio(cpm, ~is_b & is_case, ~is_b & ~is_case),
        "lfc_condition_type_b": _log2_ratio(cpm, is_b & is_case, is_b & ~is_case),
    }
    null = (truth["gene_class"] == "null").to_numpy()
    for column, observed in contrasts.items():
        planted = truth[column].to_numpy()
        hit = np.abs(planted) >= PLANTED
        assert hit.sum() > 0, column
        assert (np.sign(observed[hit]) == np.sign(planted[hit])).all(), column
        assert (np.abs(observed[null]) < PLANTED).all(), column

    log_cpm = np.log2(cpm + 1)
    centred = log_cpm - log_cpm.mean(axis=1, keepdims=True)
    _, _, vt = np.linalg.svd(centred, full_matrices=False)
    side_b, side_a = set(np.sign(vt[0][is_b])), set(np.sign(vt[0][~is_b]))
    assert len(side_b) == len(side_a) == 1
    assert side_b != side_a, "PC1 does not separate the cell types"
