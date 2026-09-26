#!/usr/bin/env python3
# /// codebase-meta
# exempt: >-
#   One-off Phase 1 feasibility spike (TODO.md P1.5), not a pipeline step. It decides whether
#   the JupyterLite notebooks run pydeseq2 live; its record is a deliverable, read by no step.
# ///
"""Pyodide feasibility spike: can the JupyterLite notebooks run pydeseq2 live?

``decide()`` is the rule, committed before any run so the record cannot shape it.
``agreement()`` turns two fits into the two numbers the rule reads.
"""

import math
from collections.abc import Mapping, Sequence

import numpy as np

Record = Mapping[tuple[str, str, str], str]


def agreement(
    reference_padj: Sequence[float],
    reference_lfc: Sequence[float],
    padj: Sequence[float],
    lfc: Sequence[float],
) -> tuple[int, float]:
    """Compare a fit with the CPython reference fit, gene by gene in the same order.

    A gene disagrees when its padj is NaN in one fit only, or when neither padj lies in the
    closed band [0.04, 0.06] and the two fall on different sides of 0.05. A gene with either
    padj inside the band is exempt, which is what the band is for.

    Args:
        reference_padj: BH-adjusted p-values of the CPython fit.
        reference_lfc: log2 fold changes of the CPython fit.
        padj: BH-adjusted p-values of the fit under test, same gene order.
        lfc: log2 fold changes of the fit under test, same gene order.

    Returns:
        The number of disagreeing genes, and max |dLFC| over genes with an LFC in both fits.
        An LFC present in one fit only makes max |dLFC| infinite.

    Raises:
        ValueError: If the four inputs are not the same length.
    """
    ref_p, new_p, ref_l, new_l = (
        np.asarray(x, dtype=float) for x in (reference_padj, padj, reference_lfc, lfc)
    )
    if not ref_p.shape == new_p.shape == ref_l.shape == new_l.shape:
        raise ValueError("the fits must cover the same genes in the same order")
    # Inside the function, so the source hash recorded with the spike covers the thresholds.
    low, alpha, high = 0.04, 0.05, 0.06
    one_nan = np.isnan(ref_p) != np.isnan(new_p)
    in_band = ((ref_p >= low) & (ref_p <= high)) | ((new_p >= low) & (new_p <= high))
    # NaN < alpha is False on both sides, so a gene NaN in both fits never lands here.
    other_side = ~one_nan & ~in_band & ((ref_p < alpha) != (new_p < alpha))
    if (np.isnan(ref_l) != np.isnan(new_l)).any():
        return int(one_nan.sum() + other_side.sum()), math.inf
    diff = np.abs(ref_l - new_l)
    diff = diff[~np.isnan(diff)]
    return int(one_nan.sum() + other_side.sum()), float(diff.max()) if diff.size else 0.0


def decide(record: Record) -> dict[str, tuple[str, str]]:
    """Apply the P1.5 rule to a spike record; the rule was fixed before the spike ran.

    Live pydeseq2 in JupyterLite at a grid size only if all of these hold: micropip installs
    pydeseq2; ``~ cell_type + condition`` fits; the fit agrees with a CPython fit of the same
    pydeseq2 version (no disagreeing gene as ``agreement()`` counts them, max |dLFC| <= 1e-3);
    wheels plus data are <= 50 MB, taken as 50,000,000 bytes; install plus fit take <= 120 s
    in the browser. Any failure gives the D4 fallback: precomputed grids plus live scipy tests.
    Stage A holding with no stage-B timing gives a provisional verdict. A missing or NaN
    metric fails its criterion, because every comparison is written so NaN reads as False.

    Args:
        record: (stage, component, field) to value, the spike's long table as a mapping.

    Returns:
        For each size in the record's grid, (verdict, why), with verdict ``live``,
        ``fallback`` or ``provisional``. Key ``plotting`` holds (``plotly`` or
        ``matplotlib``, why): the lite notebooks use matplotlib if plotly does not install.
    """

    def num(key: tuple[str, str, str]) -> float:
        value = record.get(key)
        return math.nan if value is None else float(value)

    # Inside the function, so the source hash recorded with the spike covers the thresholds.
    max_dlfc, max_bytes, max_seconds = 1e-3, 50_000_000, 120.0
    a = "stage_a"
    plotly = record.get((a, "plotly", "status")) == "installed"
    verdicts = {
        "plotting": ("plotly", "plotly installed")
        if plotly
        else ("matplotlib", "plotly did not install")
    }
    for size in record[("meta", "grid", "sizes")].split(";"):
        fit = (a, f"fit_{size}")
        wheels_data = num((a, "pydeseq2_wheels", "total_bytes")) + num(
            ("meta", f"input_{size}", "total_bytes")
        )
        browser = ("stage_b", f"browser_{size}", "install_plus_fit_s")
        criteria = (
            (
                record.get((a, "pydeseq2", "status")) == "installed",
                "micropip did not install pydeseq2",
            ),
            (
                record.get((a, "pydeseq2", "version"))
                == record.get(("cpython", "pydeseq2", "version")),
                "the CPython reference is not the same pydeseq2 version",
            ),
            (
                record.get((*fit, "status")) == "ok",
                f"~ cell_type + condition fit: {record.get((*fit, 'status'), 'missing')}",
            ),
            (
                num((*fit, "padj_side_disagreements_n")) == 0,
                f"{num((*fit, 'padj_side_disagreements_n')):g} genes on the other side of 0.05",
            ),
            (
                num((*fit, "max_abs_dlfc_log2")) <= max_dlfc,
                f"max |dLFC| {num((*fit, 'max_abs_dlfc_log2')):g} > {max_dlfc:g}",
            ),
            (wheels_data <= max_bytes, f"wheels plus data {wheels_data:.0f} > {max_bytes} bytes"),
            (
                browser not in record or num(browser) <= max_seconds,
                f"install plus fit {num(browser):g} s in the browser > {max_seconds:g} s",
            ),
        )
        failed = [why for held, why in criteria if not held]
        if failed:
            verdicts[size] = ("fallback", "; ".join(failed))
        elif browser not in record:
            verdicts[size] = ("provisional", "stage A holds; the stage-B browser timing is missing")
        else:
            verdicts[size] = ("live", "every criterion holds")
    return verdicts
