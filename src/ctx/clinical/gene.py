"""One gene's survival on the ER figures: its patients at or above its median against below.

The page draws it when the reader picks a gene (``ctxSurv.gene`` and ``ctxLive.gene``, M10.G);
this is the definition that twin is held to, and R's survival is held to this
(``tests/test_gene_survival.py``). The split is ``score.split``'s median rule, one definition in
both languages: at or above the median is the upper group, so a tie at the median goes up. It
is set before the data are seen, so its log-rank p needs no correction for a cut chosen.
"""

from dataclasses import dataclass
from typing import Any

import numpy as np
import polars as pl

from ctx.clinical import analyses, figures, score, survival
from ctx.report import text
from ctx.viz.km import KmData

PAGE = "clinical"
# The groups' keys in content/clinical.toml [gene], below first: below takes the reference
# colour and the first slot, above the comparison, as the ER figure's two groups do.
NAMES = ("below", "above")


@dataclass(frozen=True)
class GeneSplit:
    """One gene's split on one endpoint and what the figure draws of it.

    Attributes:
        months: Follow-up of the endpoint's analysed patients, in table order.
        event: 1 for an event, 0 for censored, likewise.
        upper: 1 at or above the gene's median, 0 below, likewise.
        figure: The Kaplan-Meier contract, at-risk times included.
        logrank: Chi-squared, df and p; None when one group is empty, as when more than half
            the patients share the lowest value.
    """

    months: np.ndarray
    event: np.ndarray
    upper: np.ndarray
    figure: KmData
    logrank: tuple[float, int, float] | None


def split(patients: pl.DataFrame, endpoint: str, stored: np.ndarray,
          names: tuple[str, str] = NAMES, horizon: float = 60.0) -> GeneSplit:  # fmt: skip
    """The gene's median split on the endpoint's analysed patients.

    Args:
        patients: The patient table.
        endpoint: OS, DSS, DFI or PFI.
        stored: The gene's stored values (``score.stored``) over every patient, table order.
        names: The two groups' names, below first.
        horizon: Months, for the figure contract alone.

    Returns:
        The split, its curves and its log-rank test.
    """
    analysed = patients[endpoint].is_not_null().to_numpy()
    out = analyses.outcome(patients, endpoint)
    upper = score.split(np.asarray(stored)[analysed], "median")
    labels = np.array(names)[upper]
    frame = analyses.km_frame(out.months, out.event, labels, names)
    found = figures.km(frame, list(names), "", horizon)
    both = len(np.unique(upper)) == len(names)
    test = survival.logrank(out.months, out.event, labels) if both else None
    return GeneSplit(out.months, out.event, upper, found, test)


def group_names() -> tuple[str, str]:
    """The two groups in page words, below first."""
    words = text.raw(PAGE, "gene")
    return words[NAMES[0]], words[NAMES[1]]


def controls(endpoint: str) -> dict[str, Any]:
    """An ER figure's gene box: its endpoint, the groups' names and the templates it fills.

    The expression, the patients and the screen are the survival panel's, read from its
    settings (``ctxLive.settings``): the box adds no expression to the page, only each
    screen row's two p values (``page.screen_rows``).
    """
    words = {k: v for k, v in text.raw(PAGE, "gene").items() if k not in NAMES}
    words |= text.raw(PAGE, "figures.er_gene")
    return {"kind": "gene", "endpoint": endpoint, "names": list(group_names()), "words": words,
            "defaults": {"summary": text.word(PAGE, "gene.start")}}  # fmt: skip
