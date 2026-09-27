"""The sanity rules: decide() at its boundaries, the rule schema, and (with P2C) the record."""

import math

import pytest

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
