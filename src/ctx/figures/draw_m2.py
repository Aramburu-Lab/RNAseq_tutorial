"""M2's figures redrawn from the tables ``ctx qc`` wrote, by M2's own block functions.

The run is read back as ``ctx qc`` held it when it drew: gather's sheet through the GDC door
with each run's read-length group, the registry it read with its MAD multiplier, the options
it was given, and its tables. Each drawer calls the block function that drew the figure on the
page, with a writer that writes nothing (``Recorder``), so the figure and its words are M2's
own. Nothing is measured again: the PCA, the FastQC curves and the fastp grid are tables.
"""

from collections.abc import Mapping
from dataclasses import replace
from pathlib import Path
from typing import Any

import polars as pl

from ctx.dataset import Dataset
from ctx.figures.kinds_m2 import OneMetric, OnePlot
from ctx.figures.recorder import Recorder
from ctx.gather import gdc
from ctx.qc import Fastp, blocks, checkpoint, figure_tables
from ctx.qc import figures as m2_figures
from ctx.qc import registry as reg
from ctx.qc.blocks import Reads, Stage
from ctx.qc.fastqc import Plots
from ctx.qc.figure_tables import SEPARATOR
from ctx.qc.run import Options
from ctx.report.builder import FigureText

FOLDER = "02_qc"


def _context(inputs: Mapping[str, Path], spec: Dataset) -> tuple[pl.DataFrame, reg.Registry,
                                                                  Options]:  # fmt: skip
    """The sheet, the registry and the options of the run, as ``ctx qc`` held them."""
    sheet = checkpoint.length_groups(gdc.chosen(pl.read_parquet(inputs["samplesheet"])),
                                     bool(spec.batch_bins))  # fmt: skip
    row = pl.read_parquet(inputs["figure_options"]).row(0, named=True)
    registry = replace(reg.load(inputs["registry"]), mad_multiplier=row["mad"])
    fastp = Fastp(row["fastp_q"], row["fastp_u"], row["fastp_l"], row["fastp_cut_tail"],
                  row["fastp_window"], row["fastp_mean"])  # fmt: skip
    options = Options(dataset_file=inputs["dataset_file"],
                      strata=tuple(c for c in row["strata"].split(SEPARATOR) if c),
                      peers=tuple(c for c in row["peers"].split(SEPARATOR) if c),
                      mad=row["mad"], fastp=fastp, grid_run=row["grid_run"])  # fmt: skip
    return sheet, registry, options


def _table(inputs: Mapping[str, Path], name: str) -> pl.DataFrame:
    """One of the run's tables, read back."""
    return pl.read_parquet(inputs[name])


def _stage(inputs: Mapping[str, Path]) -> Stage:
    """The judgement of the samples: values, and bands and flags where the kind reads them."""
    bands = _table(inputs, "bands") if "bands" in inputs else pl.DataFrame()
    flags = _table(inputs, "flags") if "flags" in inputs else pl.DataFrame()
    return Stage(_table(inputs, "metrics"), bands, flags)


def _drawn(recorder: Recorder, what: str) -> tuple[Any, FigureText]:
    """The figure the block function handed to the writer."""
    return recorder.one(what)


def band_map(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """Every sample's absolute band on every banded metric."""
    _, registry, _ = _context(inputs, spec)
    recorder = Recorder(FOLDER, spec.slug)
    data = m2_figures.band_map(_table(inputs, "metrics"), registry, "sample")
    if data is not None:
        blocks.band_map(recorder, data, "band_map", "sample")
    return _drawn(recorder, "band_map")


def linked_pca(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """The blind VST's PCA beside the samples' labels and QC status, from its tables."""
    sheet, _, options = _context(inputs, spec)
    share = _table(inputs, "pca_share")["share"].to_list()
    dims = (*dict.fromkeys([*options.strata, *options.peers, "batch"]),)
    data = m2_figures.linked_pca(_table(inputs, "pca_scores"), (share[0], share[1]), sheet,
                                 _table(inputs, "status"), dims)  # fmt: skip
    recorder = Recorder(FOLDER, spec.slug)
    blocks.figure(recorder, data, "linked_pca", "linked")
    return _drawn(recorder, "linked_pca")


def _identity(inputs: Mapping[str, Path], spec: Dataset, what: str) -> tuple[Any, FigureText]:
    """One of the identity figures, drawn by ``blocks.identity_figures``."""
    sheet, registry, options = _context(inputs, spec)
    recorder = Recorder(FOLDER, spec.slug)
    if "sex_check" in inputs:
        blocks.sex_figure(recorder, _table(inputs, "sex_check"), sheet)
    else:
        # The PCA is no input of these two; identity_figures reads only values and flags.
        stage = _stage(inputs)
        measured = checkpoint.Measured(stage.values, None, pl.DataFrame(), (0.0, 0.0))
        blocks.identity_figures(recorder, stage, measured, sheet, registry, options)
    return _drawn(recorder, what)


def size_genes(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """Library size against genes detected, flagged samples marked."""
    return _identity(inputs, spec, "size_genes")


def composition(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """Mitochondrial against ribosomal shares, flagged samples marked."""
    return _identity(inputs, spec, "composition")


def sex_check(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """Recorded sex against expression, from the sex check M2 wrote."""
    return _identity(inputs, spec, "sex_check")


def strip(inputs: Mapping[str, Path], p: OneMetric, spec: Dataset) -> tuple[Any, FigureText]:
    """One metric's values by stratum, its bands and flags as judged."""
    sheet, registry, options = _context(inputs, spec)
    units = sheet.unique("sample", keep="first", maintain_order=True).select(
        "sample", *blocks.grouping_columns(sheet, options), checkpoint.LENGTH_GROUP
    )
    recorder = Recorder(FOLDER, spec.slug)
    blocks.strips(recorder, _stage(inputs), units, registry, options, "sample")
    return _drawn(recorder, f"strip_{p.metric}")


def _fates(inputs: Mapping[str, Path], spec: Dataset, what: str) -> tuple[Any, FigureText]:
    """The read-fates figures, drawn by ``blocks.fate_figures``."""
    sheet, registry, options = _context(inputs, spec)
    values = _table(inputs, "metrics") if "metrics" in inputs else None
    recorder = Recorder(FOLDER, spec.slug)
    fate = _table(inputs, "fates")
    if values is None:
        dims = (*dict.fromkeys([*options.strata, *options.peers, "batch"]),)
        blocks.flow_figure(recorder, fate, sheet, dims)
    else:
        blocks.fate_figures(recorder, fate, values, sheet, registry, options)
    return _drawn(recorder, what)


def read_fates(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """Raw against useful fragments per sample."""
    return _fates(inputs, spec, "read_fates")


def read_flow(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """Where every read went, summed over the samples."""
    return _fates(inputs, spec, "read_flow")


def fastqc(inputs: Mapping[str, Path], p: OnePlot, spec: Dataset) -> tuple[Any, FigureText]:
    """One of FastQC's plots behind its verdicts, from the plot data M2 wrote."""
    sheet, registry, options = _context(inputs, spec)
    plots = figure_tables.plots_from(_table(inputs, "fastqc_curves"),
                                     _table(inputs, "fastqc_theoretical"),
                                     _table(inputs, "fastqc_overrepresented"))  # fmt: skip
    units = sheet.unique("sample", keep="first", maintain_order=True)
    recorder = Recorder(FOLDER, spec.slug)
    blocks.fastqc_figures(recorder, plots, units, registry, options)
    return _drawn(recorder, p.plot)


def _reads(inputs: Mapping[str, Path], spec: Dataset, what: str) -> tuple[Any, FigureText]:
    """The reads' grid and band map, drawn by ``blocks.read_figures``.

    FastQC's plots of the reads are no table, so the block is handed none of them.
    """
    sheet, registry, options = _context(inputs, spec)
    grid = _table(inputs, "grid")
    run = options.grid_run or ""
    top = pl.DataFrame(schema={"run_accession": pl.String, "top_pct": pl.Float64})
    none = Plots("run_accession", {}, None, top)
    skipped = pl.DataFrame(schema={"run_accession": pl.String, "status": pl.String})
    seen = Reads(_table(inputs, "read_metrics"), skipped, none, grid, run)
    recorder = Recorder(FOLDER, spec.slug)
    blocks.read_figures(recorder, seen, sheet, registry, options)
    return _drawn(recorder, what)


def fastp_grid(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """Fastp's -q by -l grid on one run."""
    return _reads(inputs, spec, "fastp_grid")


def run_band_map(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """Every run's absolute band on every banded read metric."""
    return _reads(inputs, spec, "run_band_map")
