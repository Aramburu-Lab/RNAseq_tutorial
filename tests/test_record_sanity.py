"""scripts/python/record_sanity.py writes the record the Step-1 record tests read (D-026).

The hub's own writer (frozen `ctx.report.hub.sanity`) decides the synthetic rules S1-S3; the
script must copy those verdicts byte for byte, write the union of the M5 and M6 metrics, and
write nothing when any check fails.
"""

import hashlib
import importlib.util
from datetime import date
from pathlib import Path

import polars as pl
import pytest

from ctx import runtime
from ctx.dataset import load
from ctx.de.sanity import decide, parse
from ctx.io.naming import stem
from ctx.io.readme_companion import write_table
from ctx.report.hub import sanity

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "python" / "record_sanity.py"
DE_SPEC = ROOT / "tests" / "fixtures" / "datasets" / "de_synthetic.toml"
COLUMNS = {"metric": ("none", "m"), "value": ("varies", "v"), "detail": ("none", "d")}
DE = [("calls.a_null.significant", 3.0, "a"), ("calls.a_null.tested", 900.0, "a")]
FUNCTION = [("gsea.a_case.UP_SET.nes", 2.1, "f"), ("gsea.a_case.UP_SET.fdr", 0.001, "f")]


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("record_sanity", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def metrics(root: Path, module: str, rows: list) -> Path:
    number = {"de": "05", "function": "06"}[module]
    frame = pl.DataFrame(rows, schema=["metric", "value", "detail"], orient="row")
    table = stem(module, "desynth_metrics", 1, day=date(2026, 9, 30))
    return write_table(frame, root / "outputs" / f"{number}_{module}", table,
                       columns=COLUMNS, data_releases={})[0]  # fmt: skip


@pytest.fixture
def root(tmp_path):
    metrics(tmp_path, "de", DE)
    metrics(tmp_path, "function", FUNCTION)
    sanity(runtime.Settings("laptop", tmp_path, 1, None, 0.0), "desynth", DE_SPEC)
    return tmp_path


def hub_csv(root: Path) -> Path:
    (path,) = (root / "outputs" / "09_report").glob("*_report_desynth_sanity_v1.csv")
    return path


def test_record_is_the_hub_verdicts_and_the_metrics_union(script, root, tmp_path):
    out = tmp_path / "deliverables"
    script.write(root, "desynth", DE_SPEC, out)
    (verdicts,) = out.glob("*_de_desynth_sanity_v1.csv")
    (union,) = out.glob("*_de_desynth_metrics_v1.csv")
    assert verdicts.read_bytes() == hub_csv(root).read_bytes()
    # The record tests decide again from the metrics CSV alone: S2 passes only because the
    # M6 gsea.* rows are in it, and S3 fails because no module writes its metric.
    record = script.record_from_csv(union.read_text())
    decided = decide(record, parse(load(DE_SPEC).sanity))
    assert [(r["rule"], r["verdict"]) for r in decided] == [
        ("S1", "PASS"), ("S2", "PASS"), ("S2", "PASS"), ("S3", "FAIL"),
    ]  # fmt: skip
    assert script.verdict_rows(verdicts.read_text()) == [
        (r["rule"], r["condition"], r["verdict"], r["rule_verdict"]) for r in decided
    ]
    readme = verdicts.with_name(f"{verdicts.stem}_README.md").read_text()
    assert hashlib.sha256(hub_csv(root).read_bytes()).hexdigest() in readme


def test_refuses_a_second_record(script, root, tmp_path):
    out = tmp_path / "deliverables"
    script.write(root, "desynth", DE_SPEC, out)
    with pytest.raises(FileExistsError, match="a record is written once"):
        script.write(root, "desynth", DE_SPEC, out)


def test_refuses_an_edited_hub_verdict(script, root, tmp_path):
    csv_path = hub_csv(root)
    csv_path.write_text(csv_path.read_text().replace("FAIL", "PASS"))
    out = tmp_path / "deliverables"
    with pytest.raises(ValueError, match="no record written"):
        script.write(root, "desynth", DE_SPEC, out)
    assert not out.exists() or not list(out.iterdir())


def test_refuses_metrics_changed_after_the_hub(script, root, tmp_path):
    metrics(root, "de", [(m, v + 1, d) for m, v, d in DE])
    with pytest.raises(ValueError, match="was not decided on"):
        script.write(root, "desynth", DE_SPEC, tmp_path / "deliverables")


def test_refuses_a_metric_in_both_modules(script, root, tmp_path):
    metrics(root, "function", [*FUNCTION, DE[0]])
    with pytest.raises(ValueError, match="metrics written by both"):
        script.write(root, "desynth", DE_SPEC, tmp_path / "deliverables")


def test_refuses_without_the_function_metrics(script, tmp_path):
    metrics(tmp_path, "de", DE)
    with pytest.raises(Exception, match="no desynth metrics table"):
        script.write(tmp_path, "desynth", DE_SPEC, tmp_path / "deliverables")
