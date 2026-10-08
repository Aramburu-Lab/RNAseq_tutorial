"""A selection's expression score and the splits it is grouped by, one definition each.

The score's order (GS plan, M10), mirrored by ``ctxSurv`` in the page:

1. log2(CPM + 1) from the protocol's counts, by ``ctx.de.filters.log2cpm``; library size is
   gather's, the ``N_`` rows left out.
2. Stored as int16 of the value times ``SCALE``: the build and the page read the same values.
3. Genes kept with median log2 CPM of at least ``--min-log2cpm`` and a non-zero variance.
4. z per gene over the analysed patients, from the stored values.
5. The mean over the selected genes of sign times z; a list without signs counts each as +1.
6. Standardised to SD 1, so a hazard ratio is per SD of the score.

A split is pre-specified: at or above the median is high, or tertiles by the 1/3 and 2/3
percentiles, numpy's linear percentile in both languages.
"""

import numpy as np

from ctx.clinical import ClinicalError

SCALE = 1000
INT16 = np.iinfo(np.int16)
SPLITS = ("median", "tertiles")
TERTILES = (100 / 3, 200 / 3)


def stored(values: np.ndarray) -> np.ndarray:
    """The values times ``SCALE``, rounded half to even, as int16.

    Raises:
        ClinicalError: If a value leaves int16's range.
    """
    scaled = np.rint(np.asarray(values, float) * SCALE)
    if scaled.min(initial=0) < INT16.min or scaled.max(initial=0) > INT16.max:
        raise ClinicalError("a log2 CPM times 1000 leaves int16")
    return scaled.astype(np.int16)


def expressed(values: np.ndarray, min_log2cpm: float) -> tuple[np.ndarray, np.ndarray]:
    """Which genes the filter keeps, from the stored values: (kept, zero variance), per gene.

    A gene is kept with median log2 CPM of at least ``min_log2cpm`` and a variance above zero.
    """
    real = np.asarray(values, float) / SCALE
    flat = real.var(axis=1) == 0
    return (np.median(real, axis=1) >= min_log2cpm) & ~flat, flat


def zscores(values: np.ndarray) -> np.ndarray:
    """Z per gene over the patients given, from stored values: centred, over the SD (ddof 0)."""
    real = np.asarray(values, float) / SCALE
    sd = real.std(axis=1, keepdims=True)
    return (real - real.mean(axis=1, keepdims=True)) / np.where(sd > 0, sd, np.nan)


def score(z: np.ndarray, signs: np.ndarray) -> np.ndarray:
    """The selection's score per patient: mean of sign times z over its genes, then SD 1.

    Args:
        z: ``zscores`` of the selected genes, genes by patients.
        signs: +1 or -1 per selected gene; +1 for an unsigned list.

    Raises:
        ClinicalError: If the selection is empty or its score does not vary.
    """
    if not len(z):
        raise ClinicalError("no selected gene has data on this page")
    mean = (np.asarray(signs, float)[:, None] * z).mean(axis=0)
    sd = mean.std()
    if not sd > 0:
        raise ClinicalError("the score is the same for every patient")
    return np.asarray((mean - mean.mean()) / sd)


def split(values: np.ndarray, rule: str) -> np.ndarray:
    """Group per patient: 0 low and 1 high at the median, or 0, 1, 2 by tertiles.

    At or above a cut is the higher group, the rule the scan uses at every cut.
    """
    if rule not in SPLITS:
        raise ClinicalError(f"--split {rule!r}: one of {SPLITS}")
    cuts = np.percentile(values, [50.0] if rule == "median" else list(TERTILES))
    return np.asarray((np.asarray(values)[:, None] >= cuts[None, :]).sum(axis=1))
