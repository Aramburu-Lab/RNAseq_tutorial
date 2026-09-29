"""The QC registry, `conf/qc_metrics.toml`: every metric `ctx qc` judges, its bands, their sources.

One loader in the pattern of `ctx.dataset`: tomllib, frozen dataclasses, and every problem in
the file reported at once, so a registry is fixed in one pass. The registry also writes the
MultiQC config that colours MultiQC's General Statistics table by the same bands, so the course's
page and MultiQC's report cannot disagree about a value.
"""

import tomllib
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ctx.viz.style import BAND_COLOURS

SOURCE_TYPES = frozenset({"tool_code", "tool_docs", "consortium_standard", "pipeline_default",
                          "literature", "course_convention"})  # fmt: skip
STAGES = ("reads", "trim", "quant", "counts")
KINDS = ("numeric", "status", "strandedness", "category")
DIRECTIONS = ("higher", "lower", "none")
ACTIONS = ("keep", "flag", "fix_rerun", "exclude")
SCALES = ("linear", "log10")
# FastQC's own verdicts, as the three bands.
STATUS_BANDS = {"pass": "good", "warn": "usable", "fail": "bad"}
# General Statistics shows read and base counts in millions by default, in its TSV as well, and
# colours the scaled value: 100,000 reads read as 0.1 against a 10,000 band. Multipliers of 1
# keep every count raw, in the table, the TSV and the rules alike (D-029).
RAW_COUNTS = {
    "read_count_multiplier": 1, "read_count_prefix": "", "read_count_desc": "",
    "base_count_multiplier": 1, "base_count_prefix": "", "base_count_desc": "",
}  # fmt: skip
REQUIRED = (
    "stage", "title", "meaning", "computed", "unit", "table", "column", "direction", "causes",
    "action", "matters_for_de", "read_length_tied", "source", "source_type",
)  # fmt: skip


class RegistryError(ValueError):
    """A registry file that cannot be used; the message lists every problem."""


@dataclass(frozen=True)
class Metric:
    """One QC metric: what it measures, where its value comes from, and how it is judged."""

    id: str
    stage: str
    title: str
    meaning: str
    computed: str
    unit: str
    table: str
    column: str
    module: str | None
    kind: str
    scale: str
    direction: str
    good: float | None
    usable: float | None
    rationale: str
    causes: tuple[str, ...]
    action: str
    matters_for_de: bool
    read_length_tied: bool
    source: str
    source_type: str
    stranded: float | None = None
    unstranded: float | None = None

    @property
    def banded(self) -> bool:
        """Whether the metric has absolute bands, not only the dataset's own median ± k MAD."""
        return self.good is not None

    def band(self, value: Any) -> str | None:
        """The metric's absolute band for a value: good, usable or bad; None without bands.

        Numeric bands are closed on the good side: a value exactly at `good` is good. A status
        metric maps FastQC's pass, warn and fail. Missing values have no band.
        """
        if value is None:
            return None
        if self.kind == "status":
            return STATUS_BANDS.get(str(value))
        if not self.banded or self.good is None or self.usable is None:
            return None
        v = float(value)
        if self.direction == "higher":
            return "good" if v >= self.good else "usable" if v >= self.usable else "bad"
        return "good" if v <= self.good else "usable" if v <= self.usable else "bad"

    def strandedness(self, bias: float) -> str:
        """The library type a strand bias SF / (SF + SR) shows, by nf-core/rnaseq's thresholds.

        Forward or reverse when one strand holds at least `stranded` of the reads, unstranded
        when the two differ by less than `unstranded`, and undetermined in between.
        """
        assert self.stranded is not None and self.unstranded is not None  # noqa: S101 - load() checked
        if bias >= self.stranded:
            return "forward"
        if 1 - bias >= self.stranded:
            return "reverse"
        return "unstranded" if abs(2 * bias - 1) < self.unstranded else "undetermined"


@dataclass(frozen=True)
class Registry:
    """The registry file: the dataset-relative rule and every metric."""

    path: Path
    mad_multiplier: float
    mad_scaled: bool
    min_samples: int
    source: str
    rationale: str
    metrics: Mapping[str, Metric]

    def stage(self, stage: str) -> list[Metric]:
        """The metrics of one stage, in file order."""
        return [m for m in self.metrics.values() if m.stage == stage]


def _metric(key: str, raw: Mapping[str, Any]) -> tuple[Metric | None, list[str]]:
    """One [metrics.<id>] table, or every problem it has."""
    problems = [f"[metrics.{key}] needs {f}" for f in REQUIRED if f not in raw]
    kind, scale = raw.get("kind", "numeric"), raw.get("scale", "linear")
    checks = {"stage": STAGES, "direction": DIRECTIONS, "action": ACTIONS}
    problems += [f"[metrics.{key}] {f} = {raw[f]!r}, not one of {v}"
                 for f, v in checks.items() if f in raw and raw[f] not in v]  # fmt: skip
    problems += [f"[metrics.{key}] kind = {kind!r}, not one of {KINDS}"] * (kind not in KINDS)
    problems += [f"[metrics.{key}] scale = {scale!r}, not one of {SCALES}"] * (scale not in SCALES)
    if raw.get("source_type") not in SOURCE_TYPES:
        problems.append(f"[metrics.{key}] source_type = {raw.get('source_type')!r}, not one of "
                        f"{sorted(SOURCE_TYPES)}")  # fmt: skip
    problems += _band_problems(key, raw, kind)
    if problems:
        return None, problems
    return Metric(
        id=key, stage=raw["stage"], title=raw["title"], meaning=raw["meaning"],
        computed=raw["computed"], unit=raw["unit"], table=raw["table"], column=raw["column"],
        module=raw.get("module"), kind=kind, scale=scale, direction=raw["direction"],
        good=raw.get("good"), usable=raw.get("usable"), rationale=raw.get("rationale", ""),
        causes=tuple(raw["causes"]), action=raw["action"], matters_for_de=raw["matters_for_de"],
        read_length_tied=raw["read_length_tied"], source=raw["source"],
        source_type=raw["source_type"], stranded=raw.get("stranded"),
        unstranded=raw.get("unstranded"),
    ), []  # fmt: skip


def _band_problems(key: str, raw: Mapping[str, Any], kind: str) -> list[str]:
    """What is wrong with a metric's bands, rationale and thresholds."""
    problems = []
    if raw.get("source_type") == "course_convention" and not raw.get("rationale"):
        problems.append(f"[metrics.{key}] a course_convention needs a rationale")
    if kind == "numeric":
        banded = "good" in raw or "usable" in raw
        if banded == (raw.get("bands") == "dataset_only"):
            problems.append(f'[metrics.{key}] needs good and usable, or bands = "dataset_only"')
        if banded and not {"good", "usable"} <= set(raw):
            problems.append(f"[metrics.{key}] needs both good and usable")
        if raw.get("bands") == "dataset_only" and not raw.get("rationale"):
            problems.append(f"[metrics.{key}] dataset_only needs a rationale")
        order = {"higher": 1, "lower": -1}.get(raw.get("direction", ""), 0)
        if banded and {"good", "usable"} <= set(raw) and order * (raw["good"] - raw["usable"]) < 0:
            problems.append(f"[metrics.{key}] good must lie on the good side of usable")
    if kind == "strandedness" and not {"stranded", "unstranded"} <= set(raw):
        problems.append(f"[metrics.{key}] strandedness needs stranded and unstranded")
    return problems


def load(path: Path) -> Registry:
    """Read and check a registry file.

    Raises:
        RegistryError: If the file has any problem; the message lists all of them.
    """
    raw = tomllib.loads(path.read_text())
    head = raw.get("registry", {})
    fields = ("mad_multiplier", "mad_scaled", "min_samples", "source", "rationale")
    problems = [f"[registry] needs {f}" for f in fields if f not in head]
    metrics = {}
    for key, entry in raw.get("metrics", {}).items():
        metric, found = _metric(key, entry)
        problems += found
        if metric:
            metrics[key] = metric
    if not metrics and not problems:
        problems.append("no [metrics.*] tables")
    if problems:
        raise RegistryError(f"{path}:\n  " + "\n  ".join(problems))
    return Registry(path, float(head["mad_multiplier"]), bool(head["mad_scaled"]),
                    int(head["min_samples"]), head["source"], head["rationale"],
                    metrics)  # fmt: skip


def cond_rules(metric: Metric) -> dict[str, list[dict[str, float]]]:
    """MultiQC's conditional-format rules that colour a value by the metric's band.

    MultiQC applies the classes in the order pass, warn, fail, and the last whose condition
    matches wins (checked on MultiQC 1.35, D-029); conditions in one list are alternatives.
    """
    assert metric.good is not None and metric.usable is not None  # noqa: S101 - banded only
    if metric.direction == "higher":
        return {"pass": [{"gt": metric.good}, {"eq": metric.good}],
                "warn": [{"lt": metric.good}], "fail": [{"lt": metric.usable}]}  # fmt: skip
    return {"pass": [{"lt": metric.good}, {"eq": metric.good}],
            "warn": [{"gt": metric.good}], "fail": [{"gt": metric.usable}]}  # fmt: skip


def multiqc_config(registry: Registry, comment: str) -> dict[str, Any]:
    """The MultiQC config the registry implies: band colours, rules, and a comment per module.

    Args:
        registry: The loaded registry.
        comment: The report comment, from the page text.

    Returns:
        The config, ready to write as JSON, which MultiQC reads as YAML.
    """
    banded = [m for m in registry.metrics.values()
              if m.table == "general_stats" and m.kind == "numeric" and m.banded]  # fmt: skip
    sections: dict[str, list[str]] = {}
    for m in registry.metrics.values():
        if m.module:
            sections.setdefault(m.module, []).append(f"{m.title}: {m.meaning}")
    return {
        **RAW_COUNTS,
        "report_comment": comment,
        "table_cond_formatting_colours": [
            {"pass": BAND_COLOURS["good"]}, {"warn": BAND_COLOURS["usable"]},
            {"fail": BAND_COLOURS["bad"]},
        ],
        "table_cond_formatting_rules": {m.column: cond_rules(m) for m in banded},
        "section_comments": {module: " ".join(lines) for module, lines in sections.items()},
    }  # fmt: skip
