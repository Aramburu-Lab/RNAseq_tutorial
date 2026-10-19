"""M1's figures redrawn from the sheet ``ctx gather`` wrote, by gather's own functions."""

import os
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import polars as pl

from ctx.dataset import Dataset
from ctx.gather import harmonise
from ctx.qc import registry
from ctx.report.builder import FigureText


@contextmanager
def _registry(path: Path) -> Iterator[None]:
    """Gather reads its depth standard from ``$CTX_QC_REGISTRY``: the copy given, for one call."""
    before = os.environ.get(registry.REGISTRY_ENV)
    os.environ[registry.REGISTRY_ENV] = str(path)
    try:
        yield
    finally:
        if before is None:
            del os.environ[registry.REGISTRY_ENV]
        else:
            os.environ[registry.REGISTRY_ENV] = before


def composition(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """Samples per cell type and condition, from the sheet, as gather draws them."""
    data = harmonise._composition(pl.read_parquet(inputs["samplesheet"]), spec)
    return data, harmonise._composition_text(data)


def library_sizes(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """Reads per run over ENCODE's depths, read from the QC registry copy given."""
    with _registry(inputs["qc_registry"]):
        bands, _sources = harmonise._depth_standard()
    data = harmonise._library_sizes(pl.read_parquet(inputs["samplesheet"]), spec, bands)
    return data, harmonise._library_text(data, bands[0])
