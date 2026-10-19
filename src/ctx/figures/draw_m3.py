"""M3's figures redrawn from the tables ``ctx process`` wrote, by M3's own block functions.

The run's frames are read back as the page held them when it drew: each bar's reads per
outcome, the strand table, the counters' agreement and gene sums, the sheet groups and what
the page states beside them. Nothing is aligned or counted again. Each drawer calls the block
function that drew the figure with a writer that writes nothing (``Recorder``).
"""

from collections.abc import Mapping
from pathlib import Path
from typing import Any

import polars as pl

from ctx.dataset import Dataset
from ctx.figures.kinds_m3 import OneComparison
from ctx.figures.recorder import Recorder
from ctx.process import blocks, figure_tables
from ctx.process import tables as tab
from ctx.process.page import Frames, Shown
from ctx.qc import registry as reg
from ctx.report.builder import FigureText

FOLDER = "03_process"
EMPTY = pl.DataFrame()


def _table(inputs: Mapping[str, Path], name: str) -> pl.DataFrame:
    """One of the run's tables, read back; an empty frame for one the kind does not read."""
    return pl.read_parquet(inputs[name]) if name in inputs else EMPTY


def _read(inputs: Mapping[str, Path]) -> tuple[Frames, Shown]:
    """The frames and what the page showed, from the tables the kind reads."""
    options = pl.read_parquet(inputs["figure_options"])
    row = options.row(0, named=True)
    shown = Shown(figure_tables.aligner_of(options), row["featured"],
                  _table(inputs, "groups"), row["genome_name"],
                  reg.load(inputs["registry"]).metrics["gene_strand_bias"], row["genome_genes"],
                  partial=row["partial"])  # fmt: skip
    outcomes = _table(inputs, "outcomes")
    settings = _table(inputs, "setting_outcomes")
    found = Frames(
        alignment=_table(inputs, "alignment"),
        outcomes=figure_tables.outcomes_from(outcomes) if outcomes.height else {},
        settings=figure_tables.outcomes_from(settings) if settings.height else {},
        strand=_table(inputs, "strandedness"), counting=EMPTY, gene_counts=EMPTY,
        comparison=_table(inputs, "gene_comparison"),
        agreements=tab.Agreements(_table(inputs, "concordance"), EMPTY, _table(inputs, "sums")),
        reference=EMPTY, incomplete=None, context=None,
    )  # fmt: skip
    return found, shown


def _drawn(
    inputs: Mapping[str, Path], spec: Dataset, block: Any, what: str
) -> tuple[Any, FigureText]:
    """The figure ``what`` the block function hands to the writer."""
    found, shown = _read(inputs)
    recorder = Recorder(FOLDER, spec.slug)
    block(recorder, found, shown)
    return recorder.one(what)


def mapping(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """Every sample as run, one stacked bar each."""
    return _drawn(inputs, spec, blocks.mapping_figure, "mapping")


def featured_settings(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """The featured sample as run and under each setting."""
    return _drawn(inputs, spec, blocks.settings_figure, "featured_settings")


def strand_assigned(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """What featureCounts assigns at each -s, per sample."""
    return _drawn(inputs, spec, blocks.strand_figure, "strand_assigned")


def strand_share(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """Each sample's forward share and call."""
    return _drawn(inputs, spec, blocks.share_figure, "strand_share")


def counters(inputs: Mapping[str, Path], p: OneComparison, spec: Dataset) -> tuple[Any, FigureText]:
    """One comparison's gene sums, one counter against the other."""
    return _drawn(inputs, spec, blocks.agreement_figures, f"counters_{p.comparison}")


def agreement(
    inputs: Mapping[str, Path], p: OneComparison, spec: Dataset
) -> tuple[Any, FigureText]:
    """One comparison's difference against mean, with its limits of agreement."""
    return _drawn(inputs, spec, blocks.agreement_figures, f"agreement_{p.comparison}")
