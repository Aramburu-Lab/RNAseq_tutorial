"""The cut-point scan and its corrected p, the lesson that a scanned minimum is small by chance.

Every whole percentile from ``LOW`` to ``HIGH`` splits the patients, at or above it high, and
each split gets its log-rank chi-squared. The smallest p is shown only beside the maximally
selected statistic's permutation p (Lausen and Schumacher 1992): the values are shuffled
against the outcomes ``PERMUTATIONS`` times by mulberry32 from ``SEED`` with Fisher-Yates, as
``ctxSurv.scan`` shuffles them, so the page and the build give the same corrected p exactly.
"""

from dataclasses import dataclass

import numpy as np
from scipy.stats import chi2

from ctx.clinical.survival import Splits, mulberry32

LOW, HIGH = 20, 80
PERMUTATIONS = 1000
SEED = 42


@dataclass(frozen=True)
class Scan:
    """One scan: the cuts, each one's statistic and p, and the correction.

    Attributes:
        cuts: Percentiles scanned, LOW to HIGH.
        thresholds: The values at those percentiles.
        statistic: Log-rank chi-squared per cut.
        p: Its p per cut, uncorrected.
        beaten: Permutations whose largest chi-squared reached the observed largest.
        maxima: The largest chi-squared of each permutation, in draw order.
        corrected: (1 + beaten) / (permutations + 1).
    """

    cuts: np.ndarray
    thresholds: np.ndarray
    statistic: np.ndarray
    p: np.ndarray
    beaten: int
    maxima: np.ndarray
    corrected: float

    @property
    def best(self) -> int:
        """The cut with the largest statistic, the first of equals."""
        return int(self.cuts[int(np.argmax(self.statistic))])


def scan(time: np.ndarray, event: np.ndarray, values: np.ndarray,
         permutations: int = PERMUTATIONS, seed: int = SEED) -> Scan:  # fmt: skip
    """Scan every cut of ``values`` and correct its largest statistic by permutation.

    Args:
        time: Follow-up time per patient.
        event: 1 for an event, 0 for censored.
        values: The score or expression per patient.
        permutations: Shuffles of the values against the outcomes.
        seed: mulberry32's seed.

    Returns:
        The scan.
    """
    values = np.asarray(values, float)
    cuts = np.arange(LOW, HIGH + 1)
    thresholds = np.percentile(values, cuts)
    statistic_of = Splits(time, event)

    def largest(shown: np.ndarray) -> np.ndarray:
        """Chi-squared at every cut of one ordering of the values."""
        return statistic_of(shown[:, None] >= thresholds[None, :])

    observed = largest(values)
    draw = mulberry32(seed)
    maxima = np.array([largest(draw.shuffle(values)).max() for _ in range(permutations)])
    beaten = int(np.sum(maxima >= observed.max()))
    return Scan(cuts, thresholds, observed, chi2.sf(observed, 1), beaten, maxima,
                (1 + beaten) / (permutations + 1))  # fmt: skip
