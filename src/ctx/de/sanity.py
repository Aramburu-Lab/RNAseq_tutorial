"""Pre-registered sanity rules: conditions on the metrics M5 and M6 record, and one verdict each.

The rules live in the dataset file's ``[sanity]`` table and are committed before any real
count exists; ``decide`` is pure, so a verdict can be recomputed from the record alone. A
metric that is missing or not a number is NOT MEASURED, with the reason the run gives: an absent
measurement is never a pass (D-005 records the spike rule that let one through), and it is not
a FAIL either, which would say the data contradicted the rule.
"""

import math
import operator
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import Any

OPS: dict[str, Callable[[float, float], bool]] = {
    "<": operator.lt,
    "<=": operator.le,
    ">": operator.gt,
    ">=": operator.ge,
    "==": operator.eq,
}
VERDICT_COLUMNS = {
    "rule": ("none", "Rule id from the dataset file's [sanity] table."),
    "text": ("none", "What the rule checks, in words."),
    "condition": ("none", "The condition, as metric op right-hand side."),
    "left": ("varies", "The metric's recorded value; empty when missing."),
    "right": ("varies", "The threshold, or the multiple of another metric; empty when missing."),
    "verdict": ("none", "PASS, FAIL or NOT MEASURED; RECORDED for a rule that does not gate."),
    "rule_verdict": ("none", "FAIL if any condition fails, else NOT MEASURED if any is."),
    "reason": ("none", "Why a condition is NOT MEASURED, as the run states it; else empty."),
    "source": ("none", "Where the rule comes from: the paper figure or the plan line."),
}


class RuleError(ValueError):
    """A [sanity] table that cannot be read; the message lists every problem."""


@dataclass(frozen=True)
class Condition:
    """``metric op value`` or ``metric op times * other``."""

    metric: str
    op: str
    value: float | None = None
    other: str | None = None
    times: float = 1.0

    def text(self) -> str:
        """The condition as it reads."""
        right = (
            f"{self.value:g}"
            if self.other is None
            else (self.other if self.times == 1 else f"{self.times:g} * {self.other}")
        )
        return f"{self.metric} {self.op} {right}"


@dataclass(frozen=True)
class Rule:
    """One pre-registered rule: all its conditions must hold."""

    id: str
    text: str
    source: str
    conditions: tuple[Condition, ...]
    gate: bool = True


def parse(raw: Mapping[str, Any]) -> tuple[Rule, ...]:
    """Rules from a ``[sanity]`` table with ``[[sanity.rules]]`` entries.

    Raises:
        RuleError: Listing every malformed rule or condition.
    """
    rules, problems = [], []
    for i, entry in enumerate(raw.get("rules", [])):
        rid = entry.get("id", f"rules[{i}]")
        conditions = []
        for j, c in enumerate(entry.get("conditions", [])):
            if c.get("op") not in OPS or "metric" not in c or ("value" in c) == ("other" in c):
                problems.append(
                    f"{rid} condition {j}: needs metric, an op in {list(OPS)}, "
                    "and exactly one of value or other"
                )
                continue
            conditions.append(
                Condition(
                    c["metric"],
                    c["op"],
                    float(c["value"]) if "value" in c else None,
                    c.get("other"),
                    float(c.get("times", 1.0)),
                )
            )
        if not conditions:
            problems.append(f"{rid}: no conditions")
        if not entry.get("source"):
            problems.append(f"{rid}: no source")
        rules.append(
            Rule(
                rid,
                entry.get("text", ""),
                entry.get("source", ""),
                tuple(conditions),
                bool(entry.get("gate", True)),
            )
        )
    ids = [r.id for r in rules]
    problems += [f"rule id {i!r} appears twice" for i in sorted(set(ids)) if ids.count(i) > 1]
    if problems:
        raise RuleError("\n  ".join(["[sanity]:", *problems]))
    return tuple(rules)


def _value(record: Mapping[str, float], name: str) -> float | None:
    """A recorded metric, or None when absent or not a number."""
    value = record.get(name)
    return None if value is None or math.isnan(value) else float(value)


NOT_MEASURED = "NOT MEASURED"
# The reason a condition gets when the run names none for its absent metric.
NO_VALUE = "no value recorded"


def _verdict(holds: bool | None, gate: bool) -> str:
    """PASS, FAIL or NOT MEASURED (``holds`` None) for a gating rule; RECORDED otherwise."""
    if not gate:
        return "RECORDED"
    return NOT_MEASURED if holds is None else ("PASS" if holds else "FAIL")


def decide(record: Mapping[str, float], rules: tuple[Rule, ...],
           reasons: Mapping[str, str] | None = None) -> list[dict[str, Any]]:  # fmt: skip
    """One verdict row per condition, and the rule's: FAIL, else NOT MEASURED, else PASS.

    A rule fails if any condition fails, else is NOT MEASURED if any condition's metric is
    absent, else passes.

    Args:
        record: Metric name to value, from the M5 and M6 metrics tables.
        rules: The pre-registered rules.
        reasons: Why a metric is absent, by metric name, as the run states it.

    Returns:
        Rows with the columns of ``VERDICT_COLUMNS``.
    """
    rows, reasons = [], reasons or {}
    for rule in rules:
        found = []
        for c in rule.conditions:
            left = _value(record, c.metric)
            other = _value(record, c.other) if c.other else None
            right = c.value if c.other is None else (None if other is None else c.times * other)
            absent = [m for m, v in ((c.metric, left), (c.other, right)) if m and v is None]
            holds = None if absent else OPS[c.op](left, right)  # type: ignore[arg-type]
            why = "; ".join(dict.fromkeys(reasons.get(m, NO_VALUE) for m in absent))
            found.append(
                {
                    "rule": rule.id,
                    "text": rule.text,
                    "condition": c.text(),
                    "left": left,
                    "right": right,
                    "verdict": _verdict(holds, rule.gate),
                    "source": rule.source,
                    "reason": why,
                }
            )
        verdicts = {r["verdict"] for r in found}
        overall = next((v for v in ("RECORDED", "FAIL", NOT_MEASURED) if v in verdicts), "PASS")
        rows += [{**r, "rule_verdict": overall} for r in found]
    return rows
