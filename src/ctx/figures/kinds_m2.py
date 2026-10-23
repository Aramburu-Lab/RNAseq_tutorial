"""M2's figure kinds: drawn from the tables ``ctx qc`` wrote, gather's sheet and the registry.

Every kind reads the run's options and the registry it read, so a redraw judges and colours
as the run did.
"""

from dataclasses import dataclass
from typing import Any

from ctx.dataset import Dataset
from ctx.figures.kind import DATASET_FILE, QC_REGISTRY_RUN, Inputs, Kind, Plain, once, unknown
from ctx.qc import fastqc
from ctx.qc import registry as reg

# FastQC's plots by stem, the registry metric each shows (``blocks.fastqc_figures``).
FASTQC = {"base_quality": "fastqc_per_base_quality", "fastqc_adapter": "fastqc_adapter",
          "fastqc_n_content": "fastqc_n_content", "gc_distribution": "fastqc_gc_distribution",
          "overrepresented": fastqc.OVERREPRESENTED[0]}  # fmt: skip


@dataclass(frozen=True)
class OneMetric:
    """One registry metric's strip."""

    dataset: str
    metric: str


@dataclass(frozen=True)
class OnePlot:
    """One FastQC plot, by its figure's name."""

    dataset: str
    plot: str


def _run(_: Any, __: Dataset) -> Inputs:
    """What every M2 figure reads: the sheet, the dataset file, the options, the registry."""
    return {"samplesheet": ("gather", "samplesheet"), "dataset_file": DATASET_FILE,
            "figure_options": ("qc", "figure_options"), "registry": QC_REGISTRY_RUN}  # fmt: skip


def _with(*tables: str) -> Any:
    """``_run``'s inputs and these M2 tables, each by its name."""
    return lambda p, s: {**_run(p, s), **{t: ("qc", t) for t in tables}}


def _strips(spec: Dataset) -> list[OneMetric]:
    """Every metric of the course registry; the stem names the one a run drew."""
    return [OneMetric(spec.slug, m) for m in reg.load(reg.resolve(None)).metrics]


def _plots(spec: Dataset) -> list[OnePlot]:
    """Every FastQC plot M2 draws."""
    return [OnePlot(spec.slug, p) for p in FASTQC]


def _plot_check(p: OnePlot, spec: Dataset) -> list[str]:
    """The plot is one M2 draws."""
    return unknown("plot", p.plot, FASTQC, "a FastQC plot M2 draws for", spec.slug)


DRAW = "ctx.figures.draw_m2"
KINDS = (
    Kind("qc", "band_map", Plain, "ctx.qc.figures:band_map", f"{DRAW}:band_map", once,
         lambda _: "band_map", _with("metrics")),
    Kind("qc", "linked_pca", Plain, "ctx.qc.figures:linked_pca", f"{DRAW}:linked_pca", once,
         lambda _: "linked_pca", _with("pca_scores", "pca_share", "status")),
    Kind("qc", "size_genes", Plain, "ctx.qc.figures:density", f"{DRAW}:size_genes", once,
         lambda _: "size_genes", _with("metrics", "flags")),
    Kind("qc", "composition", Plain, "ctx.qc.figures:composition", f"{DRAW}:composition", once,
         lambda _: "composition", _with("metrics", "flags")),
    Kind("qc", "sex_check", Plain, "ctx.qc.figures:sex_scatter", f"{DRAW}:sex_check", once,
         lambda _: "sex_check", _with("sex_check")),
    Kind("qc", "strip", OneMetric, "ctx.qc.figures:strip", f"{DRAW}:strip", _strips,
         lambda p: f"strip_{p.metric}", _with("metrics", "bands", "flags")),
    Kind("qc", "read_fates", Plain, "ctx.qc.figures:fates_bubble", f"{DRAW}:read_fates", once,
         lambda _: "read_fates", _with("fates", "metrics")),
    Kind("qc", "read_flow", Plain, "ctx.qc.figures:fates_flow", f"{DRAW}:read_flow", once,
         lambda _: "read_flow", _with("fates")),
    Kind("qc", "fastqc", OnePlot, "ctx.qc.fastqc:curve_figure", f"{DRAW}:fastqc", _plots,
         lambda p: p.plot,
         _with("fastqc_curves", "fastqc_theoretical", "fastqc_overrepresented"), _plot_check),
)  # fmt: skip
