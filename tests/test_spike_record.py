"""The P1.5 rule, and the Pyodide spike record it was applied to. No network."""

import ast
import csv
import hashlib
import importlib.util
import inspect
import math
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SIZES = ("1000x24", "20000x24", "20000x200")
RECORD_NAME = re.compile(r"(\d{4}-\d{2}-\d{2})_lite_build_pyodide_spike_v(\d+)\.csv")


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


def record_path():
    """The newest record: a rerun supersedes the old one under a new date or version."""
    found = [(RECORD_NAME.fullmatch(p.name), p) for p in (ROOT / "deliverables").iterdir()]
    found = [(m[1], int(m[2]), p) for m, p in found if m]
    assert found, "no spike record in deliverables/"
    return max(found)[2]


@pytest.fixture(scope="module")
def record():
    with record_path().open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert rows and list(rows[0]) == ["stage", "component", "field", "value"]
    keys = [(r["stage"], r["component"], r["field"]) for r in rows]
    assert len(set(keys)) == len(keys), "a (stage, component, field) key repeats"
    return {k: r["value"] for k, r in zip(keys, rows, strict=True)}


def test_decision_follows_rule(record):
    verdicts = spike.decide(record)
    recorded = {
        c: (record[("decision", c, "verdict")], record[("decision", c, "why")])
        for (stage, c, field) in record
        if stage == "decision" and field == "verdict"
    }
    assert verdicts == recorded
    assert set(SIZES) < set(recorded)


def test_rule_source_unchanged(record):
    for name in ("decide", "agreement"):
        source = inspect.getsource(getattr(spike, name))
        recorded = record[("meta", "rule", f"{name}_source_sha256")]
        assert hashlib.sha256(source.encode()).hexdigest() == recorded, name


def git(*args):
    return subprocess.run(
        ["git", "-C", str(ROOT), *args], capture_output=True, text=True, check=True
    ).stdout.strip()


def test_rule_predates_record():
    if git("rev-parse", "--is-shallow-repository") == "true":
        pytest.skip("shallow clone: CI checks out depth 1, so the commit adding decide() is absent")
    script = "scripts/python/spike_pyodide.py"
    rule = git("log", "--reverse", "--format=%H", "-S", "def decide(", "--", script).split()
    path = record_path().relative_to(ROOT).as_posix()
    added = git("log", "--reverse", "--diff-filter=A", "--format=%H", "--", path).split()
    assert rule, "no commit adds decide()"
    assert added, f"{path} is not committed"
    assert rule[0] != added[0], "the rule and the record arrived in one commit"
    ancestor = subprocess.run(
        ["git", "-C", str(ROOT), "merge-base", "--is-ancestor", rule[0], added[0]], check=False
    )
    assert ancestor.returncode == 0, f"{rule[0]} (decide) is not an ancestor of {added[0]}"
    # Order alone is not enough: a rule rewritten between the two commits would pass. The rule
    # as committed first must be the rule in force now, whose hash the record carries.
    at_rule, now = git("show", f"{rule[0]}:{script}"), (ROOT / script).read_text()
    for name in ("decide", "agreement"):
        assert _function_source(at_rule, name) == _function_source(now, name), (
            f"{name}() changed after {rule[0][:7]}, the commit that fixed the rule"
        )


def _function_source(text, name):
    node = next(
        n for n in ast.parse(text).body if isinstance(n, ast.FunctionDef) and n.name == name
    )
    return ast.get_source_segment(text, node)


def test_record_complete(record):
    assert spike.PYODIDE_VERSION == "314.0.6"
    assert record[("meta", "pyodide", "version")] == spike.PYODIDE_VERSION
    assert record[("meta", "pyodide", "npm_version")] == spike.PYODIDE_VERSION
    assert record[("meta", "pyodide", "npm_integrity")].startswith("sha512-")

    for package in ("scipy", "plotly", "pydeseq2"):
        status = record[("stage_a", package, "status")]
        detail = "version" if status == "installed" else "error"
        assert record.get(("stage_a", package, detail)), (package, status)
    for package in ("numpy", "scipy", "scikit-learn", "pydeseq2"):
        assert record.get(("cpython", package, "version")), package
        assert record.get(("stage_a", f"deseq_dist:{package}", "version")), package

    for group in ("pydeseq2_wheels", "fallback_wheels"):
        files = [int(v) for (s, c, f), v in record.items() if c.startswith(f"{group}:")]
        assert files and all(b > 0 for b in files), group
        assert sum(files) == int(record[("stage_a", group, "total_bytes")]), group
        assert len(files) == int(record[("stage_a", group, "files_n")]), group

    assert record[("meta", "grid", "sizes")].split(";") == list(SIZES)
    synth = importlib.util.module_from_spec(
        importlib.util.spec_from_file_location("synth_signature", ROOT / "tests" / "synth.py")
    )
    synth.__spec__.loader.exec_module(synth)
    declared = set(inspect.signature(synth.counts).parameters)
    for size in SIZES:
        call = record[("meta", f"input_{size}", "counts_call")]
        keywords = {k.arg for k in ast.parse(call, mode="eval").body.keywords}
        assert keywords == declared, f"{call} leaves an argument to its default"
        # A numpy that changes counts() fails here by design; the spike is then rerun.
        for kind, blob in zip(("counts", "samples"), spike.input_csv(call), strict=True):
            digest = record[("meta", f"input_{size}", f"{kind}_sha256")]
            assert hashlib.sha256(blob).hexdigest() == digest, (size, kind)


def _results_csv(path, rows):
    header = "gene_id,baseMean,log2FoldChange,pvalue,padj\n"
    path.write_text(header + "".join(",".join(map(str, r)) + "\n" for r in rows))
    return path


def test_compare_counts_what_it_is_fed(tmp_path):
    reference = _results_csv(
        tmp_path / "ref.csv",
        [("g1", 100, 1.0, 0.01, 0.02), ("g2", 100, 1.0, 0.001, 0.01),
         ("g3", 1, 0.2, 0.5, 0.9), ("g4", 100, 2.0, 0.5, 0.8)],
    )  # fmt: skip
    pyodide = _results_csv(
        tmp_path / "new.csv",
        [("g1", 100, 1.0, 0.01, 0.02), ("g2", 100, 1.0, 0.1, 0.2),
         ("g3", 1, 0.2, 0.5, ""), ("g4", 100, 2.5, 0.5, 0.8)],
    )  # fmt: skip
    out = spike.compare("deseq", reference, pyodide)
    assert out["genes_n"] == 4
    assert out["padj_side_disagreements_n"] == 2  # g2 flips side, g3 is NaN in one fit only
    assert out["max_abs_dlfc_log2"] == pytest.approx(0.5)
    assert out["padj_nan_mismatch_n"] == 1
    shuffled = _results_csv(tmp_path / "shuffled.csv", [("g2", 1, 1, 1, 1), ("g1", 1, 1, 1, 1)])
    with pytest.raises(ValueError, match="different genes"):
        spike.compare("deseq", reference, shuffled)


def test_wheels_counts_unique_install_files():
    spec = importlib.util.spec_from_file_location(
        "pyodide_harness", ROOT / "scripts" / "python" / "pyodide_harness.py"
    )
    harness = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(harness)
    base = "https://files.example/packages"
    log = [
        {"phase": "install", "url": f"{base}/a-1.0-py3-none-any.whl", "bytes": 100},
        {"phase": "install", "url": f"{base}/a-1.0-py3-none-any.whl", "bytes": 100},
        {"phase": "install", "url": "https://pypi.org/simple/a/", "bytes": 7},
        {"phase": "install", "url": f"{base}/a-1.0-py3-none-any.whl.metadata", "bytes": 9},
        {"phase": "install", "url": f"{base}/b%2Bc-2.0.zip", "bytes": 50},
        {"phase": "import", "url": f"{base}/late-1.0-py3-none-any.whl", "bytes": 1},
    ]
    # Unique files by name: a wheel fetched twice counts once, so this is not bytes downloaded.
    assert harness.wheels(log) == {"a-1.0-py3-none-any.whl": 100, "b+c-2.0.zip": 50}
