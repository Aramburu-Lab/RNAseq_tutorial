"""The QC registry, `conf/qc_metrics.toml`: every metric `ctx qc` judges, its bands, their sources.

One loader in the pattern of `ctx.dataset`: tomllib, frozen dataclasses, and every problem in
the file reported at once, so a registry is fixed in one pass. The registry also writes the
MultiQC config that colours MultiQC's General Statistics table by the same bands, so the course's
page and MultiQC's report cannot disagree about a value.
"""

import math
import os
import re
import tomllib
from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any

from ctx.viz.style import BAND_COLOURS

SOURCE_TYPES = frozenset({"tool_code", "tool_docs", "consortium_standard", "pipeline_default",
                          "literature", "course_convention"})  # fmt: skip
STAGES = ("reads", "trim", "align", "quant", "counts")
KINDS = ("numeric", "status", "strandedness", "category")
DIRECTIONS = ("higher", "lower", "none")
ACTIONS = ("keep", "flag", "fix_rerun", "exclude")
SCALES = ("linear", "log10")
# FastQC's own verdicts, as the three bands.
STATUS_BANDS = {"pass": "good", "warn": "usable", "fail": "bad"}
# A tool's own plot shades these, best first: good from the first value, usable from the second.
PLOT_BANDS = ("good", "usable")
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
# Keys a metric may carry besides REQUIRED. Any other is refused: a misspelt key loaded
# silently and drew nothing (the proposal's plot_line, before any code read it; D-052).
OPTIONAL = (
    "module", "kind", "scale", "good", "usable", "bands", "rationale", "stranded", "unstranded",
    "plot_bands", "plot_line", "plot_source", "plot_expected_bias_positions",
    "plot_unstranded_tolerance", "judged_within", "exclude_gene_ids", "exclude_source",
    "flagged",
)  # fmt: skip
# An unversioned Ensembl gene ID, as exclude_gene_ids names one: the version a release bumps
# must not stop the match.
GENE_ID = re.compile(r"ENSG\d{11}")
# Groups a metric may be judged within instead of the --strata: read_length within its
# read-length group, never across a cohort that mixes lengths by design (Rytis 2026-10-09).
# ctx.qc.checkpoint.length_groups is the one definition of the group.
JUDGED_WITHIN = ("read_length_group",)


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
    plot_bands: tuple[float, float] | None = None
    plot_source: str = ""
    plot_line: float | None = None
    plot_expected_bias_positions: int | None = None
    plot_unstranded_tolerance: float | None = None
    judged_within: str | None = None
    exclude_gene_ids: tuple[str, ...] = ()
    exclude_source: str = ""
    # A category metric's value that flags a sample; None is ``checkpoint.FLAGGED_CATEGORY``.
    flagged: str | None = None

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

    def spans(self) -> list[tuple[str, float | None, float | None]]:
        """The absolute bands as ranges of values, (band, low, high), None at an open end.

        The ranges ``band()`` assigns; a usable band with no width, where a source gives one
        cut-off, is left out.
        """
        if not self.banded or self.good is None or self.usable is None:
            return []
        good, usable = self.good, self.usable
        if self.direction == "higher":
            ranges = [("good", good, None), ("usable", usable, good), ("bad", None, usable)]
        else:
            ranges = [("good", None, good), ("usable", good, usable), ("bad", usable, None)]
        return [r for r in ranges if r[1] is None or r[2] is None or r[1] < r[2]]

    def plot_spans(self) -> list[tuple[str, float | None, float | None]]:
        """The bands a figure shades: ``plot_bands`` where the registry has them, else ``spans()``.

        A pair's side is read from the pair: good above usable shades the high side good.
        Shading never judges; ``band()`` reads ``good`` and ``usable`` alone.
        """
        if self.plot_bands is None:
            return self.spans()
        good, usable = self.plot_bands
        if good > usable:
            return [("good", good, None), ("usable", usable, good), ("bad", None, usable)]
        return [("good", None, good), ("usable", good, usable), ("bad", usable, None)]

    def call_regions(self) -> list[tuple[str, float, float]]:
        """A strandedness metric's call regions on the forward share, (call, low, high).

        Reverse at or below 1 - stranded, forward at or above stranded, unstranded where the
        two shares differ by less than ``unstranded``; ``strandedness()`` decides the edges.
        """
        assert self.stranded is not None and self.unstranded is not None  # noqa: S101 - load() checked
        half = self.unstranded / 2
        return [("reverse", 0.0, 1 - self.stranded), ("unstranded", 0.5 - half, 0.5 + half),
                ("forward", self.stranded, 1.0)]  # fmt: skip

    def strandedness(self, bias: float | Fraction) -> str:
        """The library type a forward-strand share shows, by nf-core/rnaseq's thresholds.

        The share is Salmon's strand_mapping_bias, or featureCounts' -s 1 over -s 1 plus -s 2.

        Forward or reverse when one strand holds at least `stranded` of the reads, unstranded
        when the two differ by less than `unstranded`, and undetermined in between. Decided
        in exact arithmetic on the share as given, a Fraction of two counts as it is and a
        float as it prints, so a share and its mirror get mirrored calls; in floats, 9 of 20
        fell inside the unstranded band and 11 of 20 outside it.
        """
        assert self.stranded is not None and self.unstranded is not None  # noqa: S101 - load() checked
        if isinstance(bias, float) and math.isnan(bias):
            return "undetermined"
        share, stranded = Fraction(str(bias)), Fraction(str(self.stranded))
        if share >= stranded:
            return "forward"
        if 1 - share >= stranded:
            return "reverse"
        within = abs(2 * share - 1) < Fraction(str(self.unstranded))
        return "unstranded" if within else "undetermined"


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
    theoretical_gc: str = ""

    def stage(self, stage: str) -> list[Metric]:
        """The metrics of one stage, in file order."""
        return [m for m in self.metrics.values() if m.stage == stage]


def _metric(key: str, raw: Mapping[str, Any]) -> tuple[Metric | None, list[str]]:
    """One [metrics.<id>] table, or every problem it has."""
    problems = [f"[metrics.{key}] needs {f}" for f in REQUIRED if f not in raw]
    known = (*REQUIRED, *OPTIONAL)
    problems += [f"[metrics.{key}] unknown key {f!r}" for f in raw if f not in known]
    kind, scale = raw.get("kind", "numeric"), raw.get("scale", "linear")
    checks = {"stage": STAGES, "direction": DIRECTIONS, "action": ACTIONS}
    problems += [f"[metrics.{key}] {f} = {raw[f]!r}, not one of {v}"
                 for f, v in checks.items() if f in raw and raw[f] not in v]  # fmt: skip
    problems += [f"[metrics.{key}] kind = {kind!r}, not one of {KINDS}"] * (kind not in KINDS)
    problems += [f"[metrics.{key}] scale = {scale!r}, not one of {SCALES}"] * (scale not in SCALES)
    if raw.get("source_type") not in SOURCE_TYPES:
        problems.append(f"[metrics.{key}] source_type = {raw.get('source_type')!r}, not one of "
                        f"{sorted(SOURCE_TYPES)}")  # fmt: skip
    problems += _band_problems(key, raw, kind) + _exclusion_problems(key, raw)
    if "flagged" in raw and (kind != "category" or not isinstance(raw["flagged"], str)):
        problems.append(f"[metrics.{key}] flagged names the value of a category metric that flags")
    if raw.get("judged_within", JUDGED_WITHIN[0]) not in JUDGED_WITHIN:
        problems.append(f"[metrics.{key}] judged_within = {raw['judged_within']!r}, not one of "
                        f"{JUDGED_WITHIN}")  # fmt: skip
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
        plot_bands=tuple(raw["plot_bands"]) if "plot_bands" in raw else None,
        plot_source=raw.get("plot_source", ""), plot_line=raw.get("plot_line"),
        plot_expected_bias_positions=raw.get("plot_expected_bias_positions"),
        plot_unstranded_tolerance=raw.get("plot_unstranded_tolerance"),
        judged_within=raw.get("judged_within"),
        exclude_gene_ids=tuple(raw.get("exclude_gene_ids", ())),
        exclude_source=raw.get("exclude_source", ""), flagged=raw.get("flagged"),
    ), []  # fmt: skip


def _exclusion_problems(key: str, raw: Mapping[str, Any]) -> list[str]:
    """What is wrong with a metric's gene exclusion list: each an unversioned ID, sourced."""
    ids = raw.get("exclude_gene_ids")
    if ids is None:
        orphan = "exclude_source" in raw
        return [f"[metrics.{key}] exclude_source needs exclude_gene_ids"] * orphan
    problems = [f"[metrics.{key}] exclude_gene_ids needs an exclude_source"] * (
        not raw.get("exclude_source")
    )
    if not isinstance(ids, list) or not all(isinstance(i, str) and GENE_ID.fullmatch(i)
                                            for i in ids):  # fmt: skip
        problems.append(f"[metrics.{key}] exclude_gene_ids must list unversioned Ensembl gene "
                        f"IDs (ENSG and 11 digits), not {ids!r}")  # fmt: skip
    return problems


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
    problems += _plot_problems(key, raw, kind)
    if kind == "strandedness" and not {"stranded", "unstranded"} <= set(raw):
        problems.append(f"[metrics.{key}] strandedness needs stranded and unstranded")
    return problems


def _plot_problems(key: str, raw: Mapping[str, Any], kind: str) -> list[str]:
    """What is wrong with the shading a metric's figures draw: plot_bands, plot_line and the rest.

    A pair's order must agree with ``direction``: good above usable where higher is better,
    below it where lower is; with direction none the pair's own order says which side is good.
    """
    problems = []
    plotted = [f for f in OPTIONAL if f.startswith("plot_") and f != "plot_source" and f in raw]
    if plotted and not raw.get("plot_source"):
        problems.append(f"[metrics.{key}] {', '.join(plotted)} needs a plot_source")
    bands = raw.get("plot_bands")
    if bands is not None:
        side = {"higher": 1, "lower": -1}.get(raw.get("direction", ""), 0)
        if len(bands) != len(PLOT_BANDS) or len(set(bands)) != len(PLOT_BANDS):
            problems.append(f"[metrics.{key}] plot_bands needs two different values, "
                            "[good, usable]")  # fmt: skip
        elif side * (bands[0] - bands[1]) < 0:
            problems.append(f"[metrics.{key}] plot_bands {bands}: good must lie on the "
                            f"{raw['direction']} side of usable")  # fmt: skip
    if "plot_line" in raw and not isinstance(raw["plot_line"], int | float):
        problems.append(f"[metrics.{key}] plot_line must be one number")
    positions = raw.get("plot_expected_bias_positions")
    if positions is not None and (kind != "status" or not isinstance(positions, int)
                                  or positions < 1):  # fmt: skip
        problems.append(f"[metrics.{key}] plot_expected_bias_positions needs a status metric "
                        "and a whole number of positions from 1")  # fmt: skip
    tolerance = raw.get("plot_unstranded_tolerance")
    if tolerance is not None and (kind != "strandedness" or not 0 < tolerance < 0.5):  # noqa: PLR2004
        problems.append(f"[metrics.{key}] plot_unstranded_tolerance needs a strandedness metric "
                        "and a value between 0 and 0.5")  # fmt: skip
    return problems


REGISTRY_ENV = "CTX_QC_REGISTRY"


def resolve(given: Path | None) -> Path:
    """The registry file: the one given, else ``$CTX_QC_REGISTRY``, else the repo's own.

    Raises:
        RegistryError: If none of them exists.
    """
    repo = Path(__file__).resolve().parents[3] / "conf" / "qc_metrics.toml"
    env = os.environ.get(REGISTRY_ENV)
    path = given or (Path(env) if env else repo)
    if not path.is_file():
        raise RegistryError(f"no QC registry at {path}; pass --registry or set {REGISTRY_ENV}")
    return path


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
                    metrics, head.get("fastqc_theoretical_gc", ""))  # fmt: skip


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
        **({"fastqc_config": {"fastqc_theoretical_gc": registry.theoretical_gc}}
           if registry.theoretical_gc else {}),
    }  # fmt: skip
