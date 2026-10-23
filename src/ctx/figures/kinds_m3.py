"""M3's figure kinds: drawn from the tables ``ctx process`` wrote and the registry it read.

Every kind reads what the page states beside the tables (``figure_options``), the sheet groups
the bars filter by and the registry, whose strand bands the strand figures draw. The laptop
route's kinds and nf-core's route's (``--rnaseq-dir``, D-116) have stems of their own; the two
counter figures serve both, each route drawing its own comparisons.
"""

from dataclasses import dataclass
from typing import Any

from ctx.dataset import Dataset
from ctx.figures.kind import DATASET_FILE, Kind, Plain, once, registry_of, unknown


@dataclass(frozen=True)
class OneComparison:
    """One comparison of two counters, by its name (``tables.COMPARISONS``)."""

    dataset: str
    comparison: str


def _with(*tables: str) -> Any:
    """These M3 tables, the options, the groups, the registry and the dataset file."""
    return lambda _, __: {
        **{t: ("process", t) for t in (*tables, "figure_options", "groups")},
        "registry": registry_of("process"), "dataset_file": DATASET_FILE,
    }  # fmt: skip


def _names() -> tuple[str, ...]:
    """M3's comparisons on both routes, read when asked: their modules import the readers."""
    from ctx.process.nfcore_page import COMPARISON  # noqa: PLC0415
    from ctx.process.tables import COMPARISONS  # noqa: PLC0415

    return (*COMPARISONS, COMPARISON)


def _comparisons(spec: Dataset) -> list[OneComparison]:
    """Every comparison M3 draws; a STAR run's alone holds the tool comparison."""
    return [OneComparison(spec.slug, c) for c in _names()]


def _comparison_check(p: OneComparison, spec: Dataset) -> list[str]:
    """The comparison is one M3 draws."""
    return unknown("comparison", p.comparison, _names(), "a comparison M3 draws for",
                   spec.slug)  # fmt: skip


DRAW = "ctx.figures.draw_m3"
NFCORE = "ctx.figures.draw_m3_nfcore"
FIGURES = "ctx.process.nfcore_figures"
AGREE = _with("concordance", "sums", "gene_comparison", "outcomes")


def _nfcore(figure: str, plotting: str, *tables: str) -> Kind:
    """One figure of nf-core's route, drawn once per dataset from these M3 tables."""
    return Kind("process", figure, Plain, f"{FIGURES}:{plotting}", f"{NFCORE}:{figure}", once,
                lambda _: figure, _with(*tables))  # fmt: skip


KINDS = (
    Kind("process", "mapping", Plain, "ctx.process.figures:outcome_bars", f"{DRAW}:mapping",
         once, lambda _: "mapping", _with("outcomes", "alignment")),
    Kind("process", "featured_settings", Plain, "ctx.process.figures:outcome_bars",
         f"{DRAW}:featured_settings", once, lambda _: "featured_settings",
         _with("setting_outcomes")),
    Kind("process", "strand_assigned", Plain, "ctx.process.figures:strand_bars",
         f"{DRAW}:strand_assigned", once, lambda _: "strand_assigned",
         _with("strandedness", "alignment")),
    Kind("process", "strand_share", Plain, "ctx.process.strand:share_plot",
         f"{DRAW}:strand_share", once, lambda _: "strand_share", _with("strandedness")),
    Kind("process", "counters", OneComparison, "ctx.process.figures:counter_scatter",
         f"{DRAW}:counters", _comparisons, lambda p: f"counters_{p.comparison}", AGREE,
         _comparison_check),
    Kind("process", "agreement", OneComparison, "ctx.process.figures:agreement_plot",
         f"{DRAW}:agreement", _comparisons, lambda p: f"agreement_{p.comparison}", AGREE,
         _comparison_check),
    _nfcore("aligned_mapping", "mapping_bars", "alignment", "outcomes"),
    _nfcore("rseqc_strand", "strand_bars", "strand_check"),
    _nfcore("rseqc_share", "strand_check_bars", "strand_check"),
    _nfcore("read_distribution", "distribution_bars", "distribution_tags"),
    _nfcore("junction_annotation", "junction_bars", "junction_counts"),
    _nfcore("gene_coverage", "curves", "coverage_profile"),
    _nfcore("duplication", "curves", "duplication_curve"),
    _nfcore("library_complexity", "curves", "complexity_curve", "complexity_perfect"),
    _nfcore("biotypes", "biotype_bars", "biotype_reads"),
)  # fmt: skip
