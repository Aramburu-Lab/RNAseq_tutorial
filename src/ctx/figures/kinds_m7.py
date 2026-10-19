"""M7's figure kinds: drawn from the tables ``ctx proteo`` writes before its figures.

Every kind reads the run's options and its genes. F1 and F2 read the protein layer's presence,
its samples' groups and the RNA abundance; the correlation figures (rho by class, the worked
genes, the paired heatmap, each contrast's fold changes) read the drawn genes' matrices and the
paired samples' annotation; the cross-study figures read the protein samples' groups for their n.
"""

from dataclasses import dataclass

from ctx.dataset import Dataset
from ctx.figures.kind import DATASET_FILE, Inputs, Kind, OneContrast, Plain, unknown

MODULE = "proteo"
CROSS = "cross_dataset"
# Each input name and the M7 table it is (``<dataset>_<what>``), which
# ``ctx.proteo.figure_tables`` writes and reads back; M1's and M2's inputs already take the
# names "genes" and "figure_options".
TABLES = {"proteo_options": "figure_options", "proteo_genes": "genes",
          "abundance": "abundance", "protein_samples": "protein_samples",
          "protein_present": "protein_present", "paired_annotation": "paired_annotation",
          "drawn_rows": "drawn_rows", "worked_samples": "worked_samples"}  # fmt: skip


def folds_what(contrast: str) -> str:
    """The table of one contrast's genes tested at both layers.

    Named ``f4_``, not ``contrast_``: a reader of F4's files finds them by that stem.
    """
    return f"f4_{contrast}_folds"


def summary_what(contrast: str) -> str:
    """The table of one contrast's protein samples per group and fold-change correlation."""
    return f"f4_{contrast}_summary"


@dataclass(frozen=True)
class OneGene:
    """One worked gene's RNA against protein."""

    dataset: str
    gene: str


def _tables(*names: str) -> Inputs:
    """M7's tables by their input names."""
    return {n: (MODULE, TABLES[n]) for n in names}


def _missing(_: Plain, __: Dataset) -> Inputs:
    """F1 and F2: the genes, the layer's presence and groups, the abundance, the options."""
    return {**_tables("proteo_options", "proteo_genes", "protein_samples", "protein_present",
                      "abundance"), "dataset_file": DATASET_FILE}  # fmt: skip


def _correlation(_: object, spec: Dataset) -> Inputs:
    """The correlation figures: the genes, the drawn matrices, the annotation, the options."""
    names = ["proteo_options", "proteo_genes", "paired_annotation", "drawn_rows"]
    if spec.proteo is not None and spec.proteo.design == "tumour_normal":
        names.append("worked_samples")
    return {**_tables(*names), "dataset_file": DATASET_FILE}


def _contrast(p: OneContrast, spec: Dataset) -> Inputs:
    """One contrast's F4: the correlation inputs, its genes at both layers and its numbers."""
    return {**_correlation(p, spec),
            "contrast_folds": (MODULE, folds_what(p.contrast)),
            "contrast_summary": (MODULE, summary_what(p.contrast))}  # fmt: skip


def _cross(_: Plain, __: Dataset) -> Inputs:
    """The cross-study figures: the genes, the protein samples' groups, the options."""
    return {**_tables("proteo_options", "proteo_genes", "protein_samples"),
            "dataset_file": DATASET_FILE}  # fmt: skip


def _design(spec: Dataset) -> str | None:
    """The ``[proteo]`` design, or None for a dataset file without one."""
    return spec.proteo.design if spec.proteo else None


def _always(spec: Dataset) -> list[Plain]:
    """F1 and F2, drawn by every design."""
    return [Plain(spec.slug)] if spec.proteo else []


def _across(spec: Dataset) -> list[Plain]:
    """A figure of the designs that correlate across samples."""
    return [Plain(spec.slug)] if _design(spec) not in (None, CROSS) else []


def _studies(spec: Dataset) -> list[Plain]:
    """A figure of the cross-study design."""
    return [Plain(spec.slug)] if _design(spec) == CROSS else []


def _worked(spec: Dataset) -> list[OneGene]:
    """Every worked gene the dataset file names."""
    p = spec.proteo
    return [OneGene(spec.slug, g) for g in p.genes] if p and _across(spec) else []


def contrasts(spec: Dataset) -> list[str]:
    """The contrasts M7 draws F4 for: every one across samples, the paired one of a pair."""
    p = spec.proteo
    if p is None or p.design == CROSS:
        return []
    return [p.contrast] if p.design == "tumour_normal" and p.contrast else [
        c.name for c in spec.contrasts]  # fmt: skip


def _contrasts(spec: Dataset) -> list[OneContrast]:
    """F4 of each contrast M7 draws."""
    return [OneContrast(spec.slug, c) for c in contrasts(spec)]


def _gene_check(p: OneGene, spec: Dataset) -> list[str]:
    """The gene is one of ``[proteo] genes``."""
    named = spec.proteo.genes if spec.proteo else ()
    return unknown("gene", p.gene, named, "a [proteo] gene", spec.slug)


FIGURES = "ctx.proteo.figures"
DRAW = "ctx.figures.draw_m7"
KINDS = (
    Kind(MODULE, "quantified", Plain, f"{FIGURES}:quantified", f"{DRAW}:quantified", _always,
         lambda _: "quantified", _missing),
    Kind(MODULE, "coverage", Plain, f"{FIGURES}:coverage", f"{DRAW}:coverage", _always,
         lambda _: "coverage", _missing),
    Kind(MODULE, "rho_by_class", Plain, f"{FIGURES}:rho_by_class", f"{DRAW}:rho_by_class",
         _across, lambda _: "rho_by_class", _correlation),
    Kind(MODULE, "worked", OneGene, f"{FIGURES}:worked", f"{DRAW}:worked", _worked,
         lambda p: f"worked_{p.gene.lower()}", _correlation, _gene_check),
    Kind(MODULE, "paired_zscores", Plain, f"{FIGURES}:paired", f"{DRAW}:paired_zscores", _across,
         lambda _: "paired_zscores", _correlation),
    Kind(MODULE, "contrast", OneContrast, "ctx.proteo.densities:contrast_density",
         f"{DRAW}:contrast", _contrasts, lambda p: f"contrast_{p.contrast}", _contrast),
    Kind(MODULE, "fold_change", Plain, "ctx.proteo.densities:fold", f"{DRAW}:fold_change",
         _studies, lambda _: "fold_change", _cross),
    Kind(MODULE, "shift_by_class", Plain, f"{FIGURES}:shift", f"{DRAW}:shift_by_class", _studies,
         lambda _: "shift_by_class", _cross),
)  # fmt: skip
