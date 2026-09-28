"""The sanity rules: decide() at its boundaries, the rule schema, and (with P2C) the record."""

import csv
import math
import subprocess
from pathlib import Path

import pytest

from ctx.dataset import load
from ctx.de.sanity import RuleError, decide, parse

EPS = 1e-9


def rules(*conditions, gate=True):
    return parse(
        {"rules": [{"id": "R", "source": "test", "gate": gate, "conditions": list(conditions)}]}
    )


def verdict(record, *conditions, gate=True) -> str:
    rows = decide(record, rules(*conditions, gate=gate))
    return rows[0]["rule_verdict"]


@pytest.mark.parametrize(
    ("op", "below", "at", "above"),
    [
        ("<", "PASS", "FAIL", "FAIL"),
        ("<=", "PASS", "PASS", "FAIL"),
        (">", "FAIL", "FAIL", "PASS"),
        (">=", "FAIL", "PASS", "PASS"),
    ],
)
def test_decide_boundaries(op, below, at, above):
    condition = {"metric": "m", "op": op, "value": 0.05}
    assert verdict({"m": 0.05 - EPS}, condition) == below
    assert verdict({"m": 0.05}, condition) == at
    assert verdict({"m": 0.05 + EPS}, condition) == above


def test_times_other_metric():
    # R2's shape: the low bin at least twice the flat level.
    condition = {"metric": "low", "op": ">=", "other": "flat", "times": 2}
    assert verdict({"low": 20.0, "flat": 10.0}, condition) == "PASS"
    assert verdict({"low": 20.0 - EPS, "flat": 10.0}, condition) == "FAIL"


@pytest.mark.parametrize("record", [{}, {"m": math.nan}, {"m": 1.0}])
def test_missing_metric_fails(record):
    # An absent or NaN measurement is not a pass, on either side of the comparison.
    condition = {"metric": "m", "op": ">", "other": "n"}
    assert verdict(record, condition) == "FAIL"


def test_all_conditions_must_hold():
    both = ({"metric": "a", "op": ">", "value": 0}, {"metric": "b", "op": "<", "value": 0.05})
    assert verdict({"a": 1.0, "b": 0.01}, *both) == "PASS"
    assert verdict({"a": 1.0, "b": 0.5}, *both) == "FAIL"


def test_recorded_rule_never_gates():
    condition = {"metric": "m", "op": ">", "value": 0}
    assert verdict({"m": -1.0}, condition, gate=False) == "RECORDED"
    assert verdict({}, condition, gate=False) == "RECORDED"


def test_schema_errors_listed_together():
    raw = {
        "rules": [
            {"id": "A", "source": "s", "conditions": [{"metric": "m", "op": "~", "value": 1}]},
            {"id": "B", "conditions": [{"metric": "m", "op": ">", "value": 1, "other": "n"}]},
            {"id": "A", "source": "s", "conditions": []},
        ]
    }
    with pytest.raises(RuleError) as raised:
        parse(raw)
    message = str(raised.value)
    assert "A condition 0" in message and "B condition 0" in message
    assert "B: no source" in message and "rule id 'A' appears twice" in message


# --- the Step-1 record (P2C.3), against the rules as committed before any count (P2C.1) ------


ROOT = Path(__file__).resolve().parents[1]
SPEC = ROOT / "conf" / "datasets" / "tier1.toml"
# Everything a verdict depends on: the rules, how samples are curated and counted, the
# models, the enrichment, the PERMANOVA, and the step that decides. None may change between
# the rule commit and the record commit (critique M4).
FROZEN = (
    "conf/datasets/tier1.toml", "src/ctx/dataset.py", "src/ctx/gather/curation.py",
    "src/ctx/gather/rnaseq.py", "src/ctx/de", "src/ctx/function", "src/ctx/viz/annotate.py",
    "src/ctx/report/hub.py",
)  # fmt: skip


def git(*args: str) -> str:
    return subprocess.run(
        ["git", "-C", str(ROOT), *args], capture_output=True, text=True, check=True
    ).stdout.strip()


def records() -> tuple[Path | None, Path | None]:
    verdicts = sorted((ROOT / "deliverables").glob("*_de_tier1_sanity_v*.csv"))
    metrics = sorted((ROOT / "deliverables").glob("*_de_tier1_metrics_v*.csv"))
    return (verdicts[-1] if verdicts else None), (metrics[-1] if metrics else None)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def rule_commit() -> str:
    """The registration: the newest commit that changed the ``[sanity] registered`` line.

    A re-registration (D-013) moves it forward; it must still precede the record.
    """
    found = git("log", "--format=%H", "-G", r"^registered = ", "--", str(SPEC))
    assert found, "no commit registers the Step-1 [sanity] rules"
    return found.split()[0]


def test_registration_is_the_committed_date():
    if git("rev-parse", "--is-shallow-repository") == "true":
        pytest.skip("shallow clone: CI checks out depth 1, so the registration is absent")
    registered = load(SPEC).sanity["registered"]
    at_commit = git("show", f"{rule_commit()}:{SPEC.relative_to(ROOT).as_posix()}")
    assert f'registered = "{registered}"' in at_commit


def test_frozen_since_registration():
    """The freeze, checked on every run from the registration on, not only at the record.

    ``test_rule_unchanged`` can only look back once the record exists: a frozen file changed
    after the counts arrived would surface then, when registering again is no longer honest.
    This test fails at the change itself, committed or not.
    """
    if git("rev-parse", "--is-shallow-repository") == "true":
        pytest.skip("shallow clone: CI checks out depth 1, so the registration is absent")
    if records()[0] is not None:
        pytest.skip("the record exists; test_rule_unchanged holds the freeze from here")
    # Against the working tree, so an uncommitted or untracked change counts too.
    changed = git("diff", "--name-only", rule_commit(), "--", *FROZEN).split()
    changed += git("ls-files", "--others", "--exclude-standard", "--", *FROZEN).split()
    assert not changed, (
        f"changed since the registration {rule_commit()[:7]}: {changed}. Revert the change, "
        "or, only while no Step-1 count exists, register again: set [sanity] registered in "
        "conf/datasets/tier1.toml and record why in knowledge/decisions.md (D-013)."
    )


def need_record():
    if git("rev-parse", "--is-shallow-repository") == "true":
        pytest.skip("shallow clone: CI checks out depth 1, so the rule commit is absent")
    verdicts, metrics = records()
    if verdicts is None or metrics is None:
        pytest.skip("no Step-1 record in deliverables/ yet; it is written in P2C.3")
    return verdicts, metrics


def test_rule_predates_record():
    verdicts, _ = need_record()
    added = git("log", "--reverse", "--diff-filter=A", "--format=%H", "--",
                str(verdicts.relative_to(ROOT))).split()  # fmt: skip
    assert added, f"{verdicts.name} is not committed"
    assert rule_commit() != added[0], "the rules and the record arrived in one commit"
    ancestor = subprocess.run(["git", "-C", str(ROOT), "merge-base", "--is-ancestor",
                               rule_commit(), added[0]], check=False)  # fmt: skip
    assert ancestor.returncode == 0


def test_rule_unchanged():
    verdicts, _ = need_record()
    added = git("log", "--reverse", "--diff-filter=A", "--format=%H", "--",
                str(verdicts.relative_to(ROOT))).split()  # fmt: skip
    changed = git("diff", "--name-only", rule_commit(), added[0], "--", *FROZEN)
    assert not changed, f"changed between the rule commit and the record: {changed}"


def test_decision_follows_rule():
    verdicts, metrics = need_record()
    record = {r["metric"]: float(r["value"]) if r["value"] else float("nan") for r in rows(metrics)}
    decided = decide(record, parse(load(SPEC).sanity))
    recorded = [
        (r["rule"], r["condition"], r["verdict"], r["rule_verdict"]) for r in rows(verdicts)
    ]
    assert [
        (r["rule"], r["condition"], r["verdict"], r["rule_verdict"]) for r in decided
    ] == recorded


def test_record_complete():
    verdicts, _ = need_record()
    recorded = rows(verdicts)
    assert {r.id for r in parse(load(SPEC).sanity)} == {r["rule"] for r in recorded}
    for r in recorded:
        assert r["condition"] and r["verdict"] and r["source"], r
        assert r["left"] or r["verdict"] in {"FAIL", "RECORDED"}, r
