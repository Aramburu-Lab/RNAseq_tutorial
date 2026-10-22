"""M6's figure kinds: drawn from M6's GSEA and ORA tables, or M5's results and the gene sets."""

from dataclasses import dataclass

from ctx.dataset import Dataset
from ctx.figures.kind import (
    DATASET_FILE,
    Inputs,
    Kind,
    OneContrast,
    featured,
    m5_results,
    resource,
    slug,
)


@dataclass(frozen=True)
class OneSet:
    """One gene set in every contrast."""

    dataset: str
    gene_set: str


@dataclass(frozen=True)
class SetInContrast:
    """One gene set's running score in one contrast."""

    dataset: str
    contrast: str
    gene_set: str


def _gene_sets(spec: Dataset) -> Inputs:
    """The gene-set file ``ctx setup`` fetched for ``[function]``."""
    return {"gene_sets": resource(spec.function.gene_sets)} if spec.function else {}


def _gsea_inputs(p: OneContrast, _: Dataset) -> Inputs:
    """M6's GSEA table of the contrast and the dataset file."""
    return {"gsea": ("function", f"{p.contrast}_gsea"), "dataset_file": DATASET_FILE}


def _ora_inputs(p: OneContrast, _: Dataset) -> Inputs:
    """M6's ORA table, M5's results that give its background, and the dataset file."""
    return {"ora": ("function", f"{p.contrast}_ora"), "results": m5_results(p.contrast),
            "dataset_file": DATASET_FILE}  # fmt: skip


def _set_inputs(_: OneSet, spec: Dataset) -> Inputs:
    """Every contrast's GSEA table, in the dataset file's order, and the dataset file."""
    found = {f"gsea_{c.name}": ("function", f"{c.name}_gsea") for c in spec.contrasts}
    return {**found, "dataset_file": DATASET_FILE}


def _ranked_inputs(p: OneContrast | SetInContrast, spec: Dataset) -> Inputs:
    """M5's results of the contrast, the gene sets and the dataset file: M6 ranks them again."""
    return {"results": m5_results(p.contrast), **_gene_sets(spec),
            "dataset_file": DATASET_FILE}  # fmt: skip


def _featured(spec: Dataset) -> list[OneContrast]:
    """The featured contrast, where ``ctx function`` runs: it needs ``[de]`` and ``[function]``."""
    return featured(spec) if spec.function else []


def _sets(spec: Dataset) -> list[OneSet]:
    """The featured gene set of ``[function]``."""
    return [OneSet(spec.slug, spec.function.featured_set)] if spec.function else []


def _running(spec: Dataset) -> list[SetInContrast]:
    """The running score ``[function] running`` names."""
    return [SetInContrast(spec.slug, *spec.function.running)] if spec.function else []


def set_slug(name: str) -> str:
    """A gene set's name as M6's output names carry it (``ctx.function.run._slug``)."""
    return slug(name)


KINDS = (
    Kind("function", "gsea", OneContrast, "ctx.function.run:_gsea_bars",
         "ctx.figures.draw_m6:gsea", _featured, lambda p: f"gsea_{p.contrast}", _gsea_inputs),
    Kind("function", "ora", OneContrast, "ctx.function.run:_ora_bars", "ctx.figures.draw_m6:ora",
         _featured, lambda p: f"ora_{p.contrast}", _ora_inputs),
    Kind("function", "nes", OneSet, "ctx.function.run:_featured_bars",
         "ctx.figures.draw_m6:nes", _sets, lambda p: f"nes_{set_slug(p.gene_set)}", _set_inputs),
    Kind("function", "running", SetInContrast, "ctx.function.running:scatter",
         "ctx.figures.draw_m6:running", _running, lambda p: f"running_{p.contrast}",
         _ranked_inputs),
    Kind("function", "background", OneContrast, "ctx.function.background:bars",
         "ctx.figures.draw_m6:background", _featured, lambda p: f"background_{p.contrast}",
         _ranked_inputs),
)  # fmt: skip
