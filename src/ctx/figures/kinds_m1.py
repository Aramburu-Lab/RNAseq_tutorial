"""M1's figure kinds: drawn from the curated sample sheet ``ctx gather`` wrote."""

from typing import Any

from ctx.dataset import Dataset
from ctx.figures.kind import DATASET_FILE, QC_REGISTRY, Inputs, Kind, Plain


def _sheet(_: Any, __: Dataset) -> Inputs:
    """M1 reads its own sheet and the dataset file."""
    return {"samplesheet": ("gather", "samplesheet"), "dataset_file": DATASET_FILE}


def _curated_reads(spec: Dataset) -> list[Plain]:
    """The one figure of a kind gather draws, from a curated fetchngs sheet alone.

    A GDC sheet and a count matrix's sheet carry no ``cell_type`` or ``read_count``, and
    their pages (``ctx.gather.gdc_run``, ``ctx.gather.counts_run``) draw no figure.
    """
    return [Plain(spec.slug)] if spec.sdrf and not spec.gdc and not spec.counts else []


KINDS = (
    Kind("gather", "composition", Plain, "ctx.gather.harmonise:_composition",
         "ctx.figures.draw_m1:composition", _curated_reads, lambda _: "composition", _sheet),
    Kind("gather", "library_sizes", Plain, "ctx.gather.harmonise:_library_sizes",
         "ctx.figures.draw_m1:library_sizes", _curated_reads, lambda _: "library_sizes",
         lambda p, s: {**_sheet(p, s), "qc_registry": QC_REGISTRY}),
)  # fmt: skip
