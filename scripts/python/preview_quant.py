#!/usr/bin/env python3
# /// codebase-meta
# exempt: preview builder; draws the M4 page through ctx.quant.page on made-up Salmon outputs
# ///
"""The course preview of M4 Quantification (``ctx quant``), drawn on made-up data.

The made-up dataset has no reads, so no Salmon run. This page is drawn by the real page's own
code and words (``ctx.quant.page``, ``content/quant.toml``) from made-up Salmon outputs: a
``quant.sf`` per sample, written from the dataset's own gather counts split over made-up
isoforms, and the numbers MultiQC reads from Salmon's logs. The ``quant.sf`` files go through
ctx's tximport port as real ones do, so every figure's statistics are computed. Every
caption starts "Made-up data.", and the build refuses a figure file that does not say so,
since a figure is downloaded and shown apart from the page (D-050). ``build_preview.py``
calls ``page``.
"""

from __future__ import annotations

import re
import shutil
import tempfile
import time
from dataclasses import replace
from datetime import date
from pathlib import Path

import numpy as np
import polars as pl
import typer

from ctx import runtime
from ctx.cli import app
from ctx.dataset import load
from ctx.de.run import Writer, newest
from ctx.io.naming import MODULES, stem
from ctx.qc.checkpoint import gene_types
from ctx.quant import CountsFrom, QuantError, laptop, outdir, salmon, tables, tximport
from ctx.quant import page as pages
from ctx.quant import run as quant_run
from ctx.quant.run import Inputs, Measured, Options, featured_contrast
from ctx.report import text
from ctx.report.builder import write
from ctx.report.cards import parameter_cards
from ctx.viz.save import Figure

MODULE = PAGE = "quant"
SEED = 4
# The dataset file build_preview writes into the preview root.
SPEC = "desynth.toml"
# Made-up but plausible for bulk RNA-seq: 1 to 5 isoforms per gene, lengths around 1,800 nt,
# Salmon mapping 80 to 90 % of fragments with a few in decoys and below the score.
ISOFORMS, LENGTH, FLOOR = 5, (1_800.0, 0.6), 300
FRAGMENT = {0: 250.0, 1: 250.0}  # single-end runs take Salmon's prior, 250 nt
RATE, DECOY, BELOW = (0.80, 0.90), (0.03, 0.08), (0.01, 0.03)
# Every fifth gene is made lncRNA, so the isoform figure has two types to compare.
LNCRNA_EVERY = 5
MADE_UP = re.compile(r"Made-up", re.I)


class MadeUpWriter(Writer):
    """A Writer whose figures each start their caption "Made-up data."."""

    def marked(self, data: Figure) -> Figure:
        """The figure with the made-up mark in front of its caption, once."""
        mark = text.word(PAGE, "preview.caption")
        return data if data.caption.startswith(mark) else replace(
            data, caption=f"{mark} {data.caption}")  # fmt: skip


def made_up_quants(rng: np.random.Generator, counts: pl.DataFrame, folder: Path) -> None:
    """A quant.sf per sample: each gene's gather count split over its made-up isoforms."""
    samples = counts.columns[1:]
    names, genes, lengths = [], [], []
    for gene in counts["gene_id"]:
        for k in range(int(rng.integers(1, ISOFORMS + 1))):
            names.append(f"{gene}.t{k + 1}")
            genes.append(gene)
            lengths.append(int(max(FLOOR, rng.lognormal(np.log(LENGTH[0]), LENGTH[1]))))
    gene_index = {g: i for i, g in enumerate(counts["gene_id"])}
    rows = np.array([gene_index[g] for g in genes])
    shares = rng.dirichlet(np.ones(ISOFORMS), size=counts.height)
    first = np.searchsorted(rows, np.arange(counts.height))
    weight = shares[rows, np.arange(len(rows)) - first[rows]]
    weight /= np.bincount(rows, weights=weight)[rows]
    pl.DataFrame({"transcript_id": names, "gene_id": genes, "gene_name": genes}).write_csv(
        folder / outdir.TX2GENE, separator="\t")  # fmt: skip
    for j, sample in enumerate(samples):
        reads = counts[sample].to_numpy()[rows] * weight
        effective = np.maximum(np.array(lengths) - FRAGMENT[j % 2], 1.0)
        rate = reads / effective
        (folder / sample).mkdir()
        pl.DataFrame({"Name": names, "Length": lengths, "EffectiveLength": effective.round(3),
                      "TPM": (1e6 * rate / rate.sum()).round(6),
                      "NumReads": reads.round(3)}).write_csv(folder / sample / "quant.sf",
                                                             separator="\t")  # fmt: skip
    salmon.transcript_tpm(folder, samples)


def made_up_salmon_table(rng: np.random.Generator, port: tximport.GeneTables) -> pl.DataFrame:
    """MultiQC's Salmon table as ``outdir.salmon_table`` gives it, around each sample's reads."""
    mapped = np.rint(port.counts.sum(axis=0)).astype(np.int64)
    n = mapped.size
    processed = np.rint(mapped / rng.uniform(*RATE, n)).astype(np.int64)
    return pl.DataFrame({
        "sample": list(port.samples), "num_processed": processed, "num_mapped": mapped,
        "num_decoy_fragments": np.rint(processed * rng.uniform(*DECOY, n)).astype(np.int64),
        "num_fragments_filtered_vm": np.rint(processed * rng.uniform(*BELOW, n)).astype(np.int64),
        "library_types": ["['U']"] * n, "frag_length_mean": ["250.0"] * n,
        "frag_length_sd": ["25.0"] * n, "frag_length_source": ["prior"] * n,
        "gc_bias_correct": ["False"] * n, "seq_bias_correct": ["False"] * n,
        "opt_type": ["vb"] * n, "salmon_version": ["2.7.0"] * n,
        "num_valid_targets": [str(len(port.rows))] * n, "num_decoy_targets": ["0"] * n,
        "strand_bias": rng.uniform(0.49, 0.51, n),
    })  # fmt: skip


def check_made_up(files: list[Path]) -> None:
    """Refuse a preview figure whose SVG does not say its data are made up.

    Raises:
        RuntimeError: Naming each figure that does not.
    """
    if bare := [f.name for f in files if f.suffix == ".svg" and not MADE_UP.search(f.read_text())]:
        raise RuntimeError(f"preview figures without the made-up mark: {bare}")


def page(root: Path, dataset: str, day: date) -> Path:
    """Write the M4 preview page and its figures under ``root/outputs/04_quant``.

    Args:
        root: A ctx root holding the made-up dataset's gather outputs and dataset file.
        dataset: The dataset slug.
        day: The date every output name carries.

    Returns:
        The page.
    """
    started = time.perf_counter()
    rng = np.random.default_rng(SEED)
    gathered = root / "outputs" / "01_gather"
    # The module's own reader, so the preview takes the runs ctx quant would: a GDC sheet refused.
    path, read, problems = quant_run._sheet(root, dataset)
    if read is None or path is None:
        raise QuantError("; ".join(problems))
    inputs_read = [path, newest(gathered, dataset, "counts")]
    sheet = read.unique("sample", keep="first", maintain_order=True)
    counts = pl.read_parquet(inputs_read[1]).select("gene_id", *sheet["sample"])
    spec = load(root / SPEC)
    gtf = root / spec.resource("annotation_gtf").dest
    genes = gene_types(gtf).with_columns(gene_type=pl.when(
        pl.int_range(pl.len()) % LNCRNA_EVERY == 0).then(pl.lit("lncRNA"))
        .otherwise(pl.col("gene_type")))  # fmt: skip
    folder = Path(tempfile.mkdtemp(prefix="ctx_preview_quant_"))
    try:
        made_up_quants(rng, counts, folder)
        samples = sheet["sample"].to_list()
        port = tximport.summarize(outdir.quants(folder, samples), outdir.tx2gene(folder))
        table = made_up_salmon_table(rng, port)
        fates = tables.fates(table)
        scaled = port.counts_from(CountsFrom.LENGTH_SCALED_TPM).sum(axis=0)
        sums = tables.sums(fates, dict(zip(port.samples, map(float, scaled), strict=True)))
        matrix = laptop.counts(port, {}, {"salmon": "made up"})
        contrast, _ = featured_contrast(spec)
        inputs = Inputs(spec, inputs_read[0], sheet, None, gtf, genes, None, None, contrast)
        measured = Measured(table, fates, sums, port, matrix, folder)
        releases = {"input:samplesheet": inputs_read[0].name, "input:counts": inputs_read[1].name,
                    "made_up": f"numpy default_rng({SEED}), {Path(__file__).name}"}  # fmt: skip
        writer = MadeUpWriter(root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}",
                              dataset, day, releases)  # fmt: skip
        group = typer.main.get_command(app)
        cards = parameter_cards(group, group.commands[MODULE])  # type: ignore[attr-defined]
        settings = runtime.Settings("laptop", root, 1, None, started)
        command = f"ctx quant --dataset {dataset} --dataset-file {SPEC} --rnaseq-dir rnaseq"
        built = pages.build(pages.Context(writer, inputs, Options(root / SPEC, None), measured,
                                          cards, command, settings, made_up=True))  # fmt: skip
    finally:
        # The made-up quant.sf files are this build's scratch; the figures keep their numbers.
        shutil.rmtree(folder)
    check_made_up(writer.files)
    return write(built, writer.out / f"{stem(MODULE, f'{dataset}_overview', 1, day=day)}.html")
