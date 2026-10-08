"""The descriptive numbers every reference panel prints: n, median, IQR and tau.

External data are summaries of other studies, so the card describes them and tests nothing
(CONVENTIONS.md, Figures). Standard library and numpy only, so notebook 09 can run it in
Pyodide.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Final, Literal

import numpy as np

# numpy's default quantile, the "linear" rule of polars and of R's type 7.
QUANTILE_METHOD: Final[Literal["linear"]] = "linear"
# HPA's detection limit, 1 nTPM ("nTPM<1 in all immune cells" reads "Not detected"); the
# same 1 TPM floor serves GTEx.
DETECTED = 1.0


@dataclass(frozen=True)
class Spread:
    """n, median and interquartile range of one group's values."""

    n: int
    median: float
    q1: float
    q3: float


def spread(values: Sequence[float] | np.ndarray) -> Spread:
    """n, median, first and third quartile, by numpy's linear rule.

    Raises:
        ValueError: If there are no values: an empty group has no median, never zero.
    """
    x = np.asarray(values, dtype=float)
    if x.size == 0:
        raise ValueError("no values: an empty group has no median")
    q1, median, q3 = np.quantile(x, (0.25, 0.5, 0.75), method=QUANTILE_METHOD)
    return Spread(int(x.size), float(median), float(q1), float(q3))


def tau(values: Sequence[float] | np.ndarray) -> float | None:
    """Tissue specificity index of Yanai et al. 2005 (Bioinformatics 21:650), on log2(x + 1).

    The log comes first, then the maximum: tau on linear values is a different, larger number
    (0.9575 for NPTX2's 51 HPA consensus tissues, against 0.6710 here). 0 is even expression,
    1 expression in one tissue or cell type only.

    Below ``DETECTED`` everywhere the ratios are noise: NPTX2's HPA immune cells, all under
    0.2 nTPM, gave tau 1.00 for a gene HPA calls not detected in immune cells.

    Returns:
        tau, or None where it is undefined: fewer than two values, or none at ``DETECTED``.
    """
    raw = np.asarray(values, dtype=float)
    if raw.size < 2 or raw.max() < DETECTED:  # noqa: PLR2004 - tau needs two types to compare
        return None
    x = np.log2(raw + 1.0)
    return float(np.sum(1.0 - x / x.max()) / (x.size - 1))
