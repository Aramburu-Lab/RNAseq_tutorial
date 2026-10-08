"""One definition of counts per million, shared by M5's filter, M7 and M10 (M7 plan, P7.2)."""

import numpy as np
import pytest

from ctx.clinical import run as clinical_run
from ctx.de import filters

COUNTS = np.array([[10, 0], [30, 50], [60, 50]])


def test_log2cpm_by_hand():
    """Mutations caught: CPM over the gene sums (axis 1); the pseudo-count added before scaling;
    a library size passed in but ignored for the column sums."""
    expected = np.log2(np.array([[1e5, 0], [3e5, 5e5], [6e5, 5e5]]) + 1)
    assert filters.log2cpm(COUNTS) == pytest.approx(expected, rel=1e-15)
    library = np.array([200, 100])
    assert filters.log2cpm(COUNTS, library)[0, 0] == pytest.approx(np.log2(5e4 + 1), rel=1e-15)


def test_one_definition():
    """M10's expression and M5's filter call the shared function, not a copy of it.

    Mutation caught: a second ``log2cpm`` defined in ``ctx.clinical.score`` again."""
    assert clinical_run.log2cpm is filters.log2cpm
    assert not hasattr(clinical_run.score, "log2cpm")
