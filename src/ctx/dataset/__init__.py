"""The one reader of a dataset file (``conf/datasets/<slug>.toml``).

A dataset file says where a cohort's metadata comes from and how to read it: the SDRF join,
which columns become curation columns and how their values are spelled, the batch rule, and
the counts the curated sheet must reproduce. Every schema error is reported at once.

``load`` is in ``loader``; each table family's dataclasses and parser in a module of its own.
This module re-exports the names the rest of ctx and its tests import from ``ctx.dataset``.
"""

from ctx.dataset.analyses import Contrast, PcaSpec, Split
from ctx.dataset.checks import DatasetError
from ctx.dataset.counts import COUNT_KINDS, CountsSpec
from ctx.dataset.gdc import ENDPOINTS, ClinicalSpec, GdcSpec
from ctx.dataset.loader import Dataset, load, loaded
from ctx.dataset.proteo import ProteoSpec, ProteoTable
from ctx.dataset.sources import Column, Resource

__all__ = [
    "COUNT_KINDS", "ENDPOINTS", "ClinicalSpec", "Column", "Contrast", "CountsSpec", "Dataset",
    "DatasetError", "GdcSpec", "PcaSpec", "ProteoSpec", "ProteoTable", "Resource", "Split",
    "load", "loaded",
]  # fmt: skip
