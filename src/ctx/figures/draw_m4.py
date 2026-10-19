"""M4's figures redrawn from the tables ``ctx quant`` wrote, by M4's own functions.

Nothing is quantified again: the fates are Salmon's as M4 tabled them, and the fold changes
are M4's own, from both counts of copy Q1, with the contrast's samples as M5 picks them.
"""

from collections.abc import Mapping
from pathlib import Path
from typing import Any

import polars as pl

from ctx.dataset import Dataset
from ctx.figures.kind import OneContrast
from ctx.figures.recorder import Recorder
from ctx.quant import blocks as blk
from ctx.quant.page import shown
from ctx.report.builder import FigureText

FOLDER = "04_quant"


def fates(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """Every sample's fragments by fate, as Salmon reported them."""
    recorder = Recorder(FOLDER, spec.slug)
    found = shown(pl.read_parquet(inputs["figure_sheet"]))
    blk.fates_figure(recorder, pl.read_parquet(inputs["fates"]), found)
    return recorder.one("fates")


def estimates(inputs: Mapping[str, Path], p: OneContrast, _: Dataset) -> tuple[Any, FigureText]:
    """The contrast's fold change from estimated against length-scaled counts, every gene."""
    used = pl.read_parquet(inputs["contrast_samples"])
    tested = int((used["group"] == "tested").sum())
    data, lead = blk._estimates(pl.read_parquet(inputs["fold_changes"]), p.contrast,
                                (used.height - tested, tested), None)  # fmt: skip
    return data, blk.figure_text("estimates", lead)
