"""M8's figures redrawn from the tables ``ctx ncrna`` and the modules before it wrote.

The run is read back as ``ctx ncrna`` held it: gather's sheet once per sample through the GDC
door, the contrast's samples as M5 picks them, the counts in the sheet's order, the gene
records and options M8 wrote, M5's results. Each drawer calls the block function that drew
the figure on the page with a writer that writes nothing (``Recorder``), so the figure and its
words are M8's own. The miRNA volcano reads M8's miRNA fit back, its samples as M5's fit picks
them; nothing is fitted again. A block's other tables and figures are dropped.
"""

from collections.abc import Mapping
from pathlib import Path
from time import perf_counter
from typing import Any

import numpy as np
import polars as pl

from ctx import runtime
from ctx.dataset import Dataset
from ctx.de.deseq import Fit, select_samples
from ctx.de.filters import MIN_COUNT
from ctx.figures.recorder import Recorder
from ctx.gather import gdc
from ctx.ncrna import biotypes, figure_tables, libraries_page, mirna, mirna_page, page, targets
from ctx.ncrna.run import Inputs, Tcga
from ctx.report.builder import FigureText

FOLDER = "08_ncrna"
EMPTY = pl.DataFrame()
NO_TARGETS = targets.Targets(EMPTY, EMPTY, {}, frozenset())


def _table(inputs: Mapping[str, Path], name: str) -> pl.DataFrame:
    """One input table read back; an empty frame for one the kind does not read."""
    return pl.read_parquet(inputs[name]) if name in inputs else EMPTY


def _sheet(inputs: Mapping[str, Path]) -> pl.DataFrame:
    """Gather's sheet as ``run.resolve`` reads it: one row per sample, GDC's choice alone."""
    sheet = _table(inputs, "samplesheet")
    return gdc.chosen(sheet.unique("sample", keep="first", maintain_order=True)) if sheet.height \
        else sheet  # fmt: skip


def _context(inputs: Mapping[str, Path], spec: Dataset) -> page.Context:
    """The page's context, from the inputs the kind reads; the rest left empty."""
    assert spec.ncrna is not None  # noqa: S101 - every M8 kind needs [ncrna]
    contrast = next(c for c in spec.contrasts if c.name == spec.ncrna.featured)
    sheet = _sheet(inputs)
    used = select_samples(sheet, contrast, spec.order)[0] if sheet.height else EMPTY
    counts = _table(inputs, "counts")
    matrix = counts.select(sheet["sample"].to_list()).to_numpy() if counts.height \
        and sheet.height else np.empty((0, 0))  # fmt: skip
    options = _table(inputs, "figure_options")
    # A kind without the options draws a block that reads neither the groups nor the count.
    chosen: tuple[str, ...] = biotypes.TESTED
    min_count, gtf = MIN_COUNT, Path()
    if options.height:
        chosen = figure_tables.biotypes_of(options)
        min_count, gtf = int(options["min_count"][0]), Path(options["gtf"][0])
    found = Inputs(spec, contrast, sheet, used, _table(inputs, "gene_biotypes"), matrix,
                   _table(inputs, "results"), gtf, (), _tcga(inputs, spec, counts))  # fmt: skip
    # No block a drawer calls reads the run's settings: the miRNA fit alone took its threads.
    settings = runtime.Settings("", Path.cwd(), 1, None, perf_counter())
    return page.Context(found, chosen, min_count, settings, "", ())


def _tcga(inputs: Mapping[str, Path], spec: Dataset, counts: pl.DataFrame) -> Tcga | None:
    """The miRNA half's inputs the kind reads, as ``run._tcga`` reads them."""
    if spec.gdc is None or spec.ncrna is None or spec.ncrna.mirna is None:
        return None
    named = {n: inputs[n] for n in targets.RESOURCES if n in inputs}
    mirnas = _table(inputs, "mirna_samplesheet")
    return Tcga(spec.ncrna.mirna, gdc.chosen(mirnas) if mirnas.height else mirnas,
                _table(inputs, "samplesheet"), _table(inputs, "mirna_mature_counts"), EMPTY,
                EMPTY, counts, EMPTY, targets.read(named) if named else NO_TARGETS)  # fmt: skip


def _drawn(recorder: Recorder, what: str) -> tuple[Any, FigureText]:
    """The figure the block function handed to the writer."""
    return recorder.one(what)


def _library(inputs: Mapping[str, Path], spec: Dataset, what: str) -> tuple[Any, FigureText]:
    """F1 or F2, drawn by ``page.library``."""
    recorder = Recorder(FOLDER, spec.slug)
    page.library(recorder, _context(inputs, spec))
    return _drawn(recorder, what)


def biotype_shares(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """Each sample's counted reads by biotype group."""
    return _library(inputs, spec, "biotype_shares")


def filter_by_biotype(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """The value M5's filter tests, by biotype, at the run's minimum count."""
    return _library(inputs, spec, "filter_by_biotype")


def volcano_by_biotype(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """The featured contrast by biotype, BH within each."""
    recorder = Recorder(FOLDER, spec.slug)
    page.de(recorder, _context(inputs, spec))
    return _drawn(recorder, "volcano_by_biotype")


def lncrna_neighbour(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """A lncRNA against its nearest protein-coding gene, the pairs computed again."""
    recorder = Recorder(FOLDER, spec.slug)
    page.lncrna(recorder, _context(inputs, spec))
    return _drawn(recorder, "lncrna_neighbour")


def mirna_volcano(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """The mature miRNAs' volcano, from the fit M8 wrote and its samples as M5's fit picks them."""
    assert spec.ncrna is not None and spec.ncrna.mirna is not None  # noqa: S101 - the kind's
    contrast = spec.ncrna.mirna
    table = _table(inputs, "mirna_results")
    # ``mirna.fit_mature`` hands M5's fit no level orders, so its samples are picked without.
    used, dropped = select_samples(gdc.chosen(_table(inputs, "mirna_samplesheet")), contrast, {})
    fitted = Fit(contrast, table.drop("gene_name"), used, dropped, EMPTY, EMPTY)
    recorder = Recorder(FOLDER, spec.slug)
    mirna_page.mirna_volcano(recorder, fitted, table.select("gene_id", "gene_name"))
    return _drawn(recorder, "mirna_volcano")


def _tcga_of(ctx: page.Context) -> Tcga:
    """The context's miRNA half, which every miRNA kind's dataset file declares."""
    assert ctx.inputs.tcga is not None  # noqa: S101 - _mirna makes these kinds for GDC alone
    return ctx.inputs.tcga


def target_shift(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """The first featured miRNA's targets' fold changes against the non-targets'."""
    ctx = _context(inputs, spec)
    recorder = Recorder(FOLDER, spec.slug)
    mirna_page.shift(recorder, ctx, _tcga_of(ctx), _table(inputs, "mirna_results"))
    return _drawn(recorder, "target_shift")


def mirna_target_pairs(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """Each featured miRNA against its targets across tumours, both layers of each sample."""
    ctx = _context(inputs, spec)
    tcga = _tcga_of(ctx)
    recorder = Recorder(FOLDER, spec.slug)
    mirna_page.pairs(recorder, ctx, tcga, mirna.layers(tcga.sheet, tcga.mature, tcga.counts))
    return _drawn(recorder, "mirna_target_pairs")


def library_types(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """Each tumour's rRNA-depleted library against its other one, every file of the sheet."""
    assert spec.gdc is not None  # noqa: S101 - _mirna makes this kind for GDC alone
    ctx = _context(inputs, spec)
    recorder = Recorder(FOLDER, spec.slug)
    libraries_page.group(recorder, _tcga_of(ctx), ctx.inputs.genes, spec.gdc.tumour)
    return _drawn(recorder, "library_types")
