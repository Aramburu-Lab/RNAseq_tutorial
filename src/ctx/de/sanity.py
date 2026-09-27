"""Pre-registered sanity rules: conditions on the metrics M5 and M6 record, and one verdict each.

The rules live in the dataset file's ``[sanity]`` table and are committed before any real
count exists; ``decide`` is pure, so a verdict can be recomputed from the record alone. A
metric that is missing or not a number fails a gating rule: an absent measurement is not a
pass (D-005 records the spike rule that let one through).
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
    "verdict": ("none", "PASS or FAIL; RECORDED for a rule that does not gate."),
    "rule_verdict": ("none", "FAIL if any condition of the rule fails."),
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


def decide(record: Mapping[str, float], rules: tuple[Rule, ...]) -> list[dict[str, Any]]:
    """One verdict row per condition; the rule fails if any of its conditions fails.

    Args:
        record: Metric name to value, from the M5 and M6 metrics tables.
        rules: The pre-registered rules.

    Returns:
        Rows with the columns of ``VERDICT_COLUMNS``.
    """
    rows = []
    for rule in rules:
        found = []
        for c in rule.conditions:
            left = _value(record, c.metric)
            other = _value(record, c.other) if c.other else None
            right = c.value if c.other is None else (None if other is None else c.times * other)
            holds = left is not None and right is not None and OPS[c.op](left, right)
            verdict = ("PASS" if holds else "FAIL") if rule.gate else "RECORDED"
            found.append(
                {
                    "rule": rule.id,
                    "text": rule.text,
                    "condition": c.text(),
                    "left": left,
                    "right": right,
                    "verdict": verdict,
                    "source": rule.source,
                }
            )
        failed = any(r["verdict"] == "FAIL" for r in found)
        overall = ("FAIL" if failed else "PASS") if rule.gate else "RECORDED"
        rows += [{**r, "rule_verdict": overall} for r in found]
    return rows
