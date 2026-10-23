"""M3's figures on nf-core's alignment, redrawn from the tables ``ctx process --rnaseq-dir`` wrote.

The run's frames are read back as the page held them when it drew, and each drawer calls the
block that drew the figure (``ctx.process.nfcore_blocks``) with a writer that writes nothing
(``Recorder``). Nothing of nf-core's outdir is read again.
"""

from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any

import polars as pl

from ctx.dataset import Dataset
from ctx.figures.recorder import Recorder
from ctx.process import figure_tables
from ctx.process import nfcore_blocks as blocks
from ctx.process.nfcore_page import STRAND_SHARE
from ctx.qc import registry as reg
from ctx.report.builder import FigureText

FOLDER = "03_process"
# Each table the page wrote and the frame of ``Found`` it fills.
FRAMES = {"alignment": "alignment", "strand_check": "strand", "distribution_tags": "distribution",
          "junction_counts": "junctions", "coverage_profile": "coverage",
          "duplication_curve": "duplication", "complexity_curve": "complexity",
          "complexity_perfect": "perfect", "biotype_reads": "biotypes"}  # fmt: skip


def _read(inputs: Mapping[str, Path]) -> tuple[blocks.Found, blocks.Shown]:
    """The frames the kind reads and what the page showed, from the run's tables."""
    options = pl.read_parquet(inputs["figure_options"]).row(0, named=True)
    frames: dict[str, Any] = {attr: pl.read_parquet(inputs[t]) for t, attr in FRAMES.items()
                              if t in inputs}  # fmt: skip
    if "outcomes" in inputs:
        frames["outcomes"] = figure_tables.outcomes_from(pl.read_parquet(inputs["outcomes"]))
    metric = reg.load(inputs["registry"]).metrics[STRAND_SHARE]
    shown = blocks.Shown(pl.read_parquet(inputs["groups"]), options["colour"], metric)
    return blocks.Found(**frames), shown


def _drawer(block: Callable[..., Any], what: str) -> Callable[..., tuple[Any, FigureText]]:
    """A drawer of one figure: ``block`` on the frames read back, its figure ``what``."""

    def draw(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
        """The figure the block hands to the writer."""
        found, shown = _read(inputs)
        recorder = Recorder(FOLDER, spec.slug)
        block(recorder, found, shown)
        return recorder.one(what)

    draw.__doc__ = f"``{block.__name__}`` redrawn: the figure ``{what}``."
    return draw


aligned_mapping = _drawer(blocks.mapping_figure, "aligned_mapping")
rseqc_strand = _drawer(blocks.strand_figure, "rseqc_strand")
rseqc_share = _drawer(blocks.strand_check_figure, "rseqc_share")
read_distribution = _drawer(blocks.distribution_figure, "read_distribution")
junction_annotation = _drawer(blocks.junction_figure, "junction_annotation")
gene_coverage = _drawer(blocks.coverage_figure, "gene_coverage")
duplication = _drawer(blocks.duplication_figure, "duplication")
library_complexity = _drawer(blocks.complexity_figure, "library_complexity")
biotypes = _drawer(blocks.biotype_figure, "biotypes")
