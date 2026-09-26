"""The P1.5 rule, and the Pyodide spike record it was applied to. No network."""

import importlib.util
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SIZES = ("1000x24", "20000x24", "20000x200")


def load_spike():
    spec = importlib.util.spec_from_file_location(
        "spike_pyodide", ROOT / "scripts" / "python" / "spike_pyodide.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules["spike_pyodide"] = module
    spec.loader.exec_module(module)
    return module


spike = load_spike()


def passing_record():
    """Every criterion holds at every size, each metric at its limit where it has one."""
    record = {
        ("meta", "grid", "sizes"): ";".join(SIZES),
        ("stage_a", "plotly", "status"): "installed",
        ("stage_a", "pydeseq2", "status"): "installed",
        ("stage_a", "pydeseq2", "version"): "0.5.4",
        ("cpython", "pydeseq2", "version"): "0.5.4",
        ("stage_a", "pydeseq2_wheels", "total_bytes"): "40000000",
    }
    for size in SIZES:
        record |= {
            ("stage_a", f"fit_{size}", "status"): "ok",
            ("stage_a", f"fit_{size}", "padj_side_disagreements_n"): "0",
            ("stage_a", f"fit_{size}", "max_abs_dlfc_log2"): "0.001",
            ("meta", f"input_{size}", "total_bytes"): "10000000",
            ("stage_b", f"browser_{size}", "install_plus_fit_s"): "120",
        }
    return record


# One entry per criterion: the change that breaks it, and whether it breaks every size.
FAILURES = {
    "not installed": ({("stage_a", "pydeseq2", "status"): "failed"}, True),
    "other version": ({("cpython", "pydeseq2", "version"): "0.5.3"}, True),
    "fit timeout": ({("stage_a", "fit_{size}", "status"): "timeout"}, False),
    "one disagreement": ({("stage_a", "fit_{size}", "padj_side_disagreements_n"): "1"}, False),
    "dlfc over": ({("stage_a", "fit_{size}", "max_abs_dlfc_log2"): "0.0010001"}, False),
    "dlfc nan": ({("stage_a", "fit_{size}", "max_abs_dlfc_log2"): "nan"}, False),
    "bytes over": ({("meta", "input_{size}", "total_bytes"): "10000001"}, False),
    "browser over": ({("stage_b", "browser_{size}", "install_plus_fit_s"): "120.001"}, False),
}


def test_decide_boundaries():
    verdicts = spike.decide(passing_record())
    assert {s: verdicts[s][0] for s in SIZES} == dict.fromkeys(SIZES, "live")
    assert verdicts["plotting"][0] == "plotly"

    for name, (change, global_) in FAILURES.items():
        record = passing_record()
        record |= {(st, c.format(size="20000x24"), f): v for (st, c, f), v in change.items()}
        verdicts = spike.decide(record)
        expect = dict.fromkeys(SIZES, "fallback" if global_ else "live")
        expect["20000x24"] = "fallback"
        assert {s: verdicts[s][0] for s in SIZES} == expect, name

    record = passing_record()
    del record[("stage_a", "fit_1000x24", "max_abs_dlfc_log2")]
    assert spike.decide(record)["1000x24"][0] == "fallback"

    for size in SIZES:
        del record[("stage_b", f"browser_{size}", "install_plus_fit_s")]
    verdicts = spike.decide(record)
    assert verdicts["1000x24"][0] == "fallback", "a stage-A failure is final without stage B"
    assert verdicts["20000x24"][0] == verdicts["20000x200"][0] == "provisional"

    record[("stage_a", "plotly", "status")] = "failed"
    assert spike.decide(record)["plotting"][0] == "matplotlib"


def test_agreement_boundaries():
    nan = math.nan
    ref = [0.039, 0.07, 0.06, 0.061, 0.01, 0.9, nan, 0.02]
    new = [0.061, 0.04, 0.03, 0.039, 0.02, 0.8, nan, nan]
    # Outside the band on opposite sides: g1 and g4. g2 and g3 are exempt, each on an edge
    # of the closed band in one fit only. g7 is NaN in both, g8 in one only, which counts.
    lfc = [0.0] * len(ref)
    assert spike.agreement(ref, lfc, new, lfc) == (3, 0.0)
    assert spike.agreement([0.5, nan], [1.0, 2.0], [0.5, nan], [1.0005, 2.0])[1] == pytest.approx(
        5e-4
    )
    assert spike.agreement([0.5], [1.0], [0.5], [nan])[1] == math.inf
    with pytest.raises(ValueError, match="same genes"):
        spike.agreement([0.5], [1.0], [0.5, 0.5], [1.0, 1.0])
