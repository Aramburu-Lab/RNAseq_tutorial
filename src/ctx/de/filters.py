"""The expression pre-filter: edgeR's ``filterByExpr``, on the contrast's own samples.

Ported from edgeR 4.10.5 ``R/filterByExpr.R`` (``filterByExpr.default`` with ``group`` given
and ``lib.size`` left to the column sums; Bioconductor release tarball SHA-256
083f1ec9ff3999d4df497927c63d2f83d0f8a59b28510a32cdc92aa0834950c7, read 2026-09-27).
"""

import numpy as np

MIN_COUNT = 10
MIN_TOTAL_COUNT = 15
LARGE_N = 10
MIN_PROP = 0.7
# edgeR's tolerance on the two comparisons.
TOLERANCE = 1e-14


def min_sample_size(group_sizes: list[int]) -> float:
    """How many samples must pass the CPM cutoff: the smallest group, eased above ``LARGE_N``."""
    smallest = float(min(n for n in group_sizes if n > 0))
    return LARGE_N + (smallest - LARGE_N) * MIN_PROP if smallest > LARGE_N else smallest


def expressed(counts: np.ndarray, group_sizes: list[int], min_count: int = MIN_COUNT) -> np.ndarray:
    """Which genes to test: enough counts per million in enough samples, and enough in total.

    The CPM cutoff is ``min_count`` reads at the median library size, so it scales with
    depth: a raw-count threshold would keep every gene of a deep library and drop the same
    gene in a shallow one.

    Args:
        counts: Genes by samples, the contrast's samples only.
        group_sizes: Samples per compared group.
        min_count: Reads a gene needs at the median library size.

    Returns:
        A boolean mask over genes.
    """
    library = counts.sum(axis=0)
    cutoff = min_count / np.median(library) * 1e6
    cpm = counts / library * 1e6
    enough = (cpm >= cutoff).sum(axis=1) >= min_sample_size(group_sizes) - TOLERANCE
    total = counts.sum(axis=1) >= MIN_TOTAL_COUNT - TOLERANCE
    return np.asarray(enough & total, dtype=bool)
