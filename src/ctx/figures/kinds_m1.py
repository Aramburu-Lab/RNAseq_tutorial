"""M1's figure kinds: drawn from the curated sample sheet ``ctx gather`` wrote."""

from typing import Any

from ctx.dataset import Dataset
from ctx.figures.kind import DATASET_FILE, QC_REGISTRY, Inputs, Kind, Plain, once


def _sheet(_: Any, __: Dataset) -> Inputs:
    """M1 reads its own sheet and the dataset file."""
    return {"samplesheet": ("gather", "samplesheet"), "dataset_file": DATASET_FILE}


KINDS = (
    Kind("gather", "composition", Plain, "ctx.gather.harmonise:_composition",
         "ctx.figures.draw_m1:composition", once, lambda _: "composition", _sheet),
    Kind("gather", "library_sizes", Plain, "ctx.gather.harmonise:_library_sizes",
         "ctx.figures.draw_m1:library_sizes", once, lambda _: "library_sizes",
         lambda p, s: {**_sheet(p, s), "qc_registry": QC_REGISTRY}),
)  # fmt: skip
