"""M4's laptop path, ahead of Salmon: M3's trimmed reads, and the gene counts it gives M4.

The reads are M3's own: the newest ``folders`` table ``ctx process`` wrote for the dataset
names its trim folder, whose manifest lists every trimmed file at its size. A sample's runs
are quantified together, as M3 aligns them, single-end only, since Salmon's ``-r`` takes one
mate. Nothing here writes; a missing table, folder or file is a problem listed with the rest.
"""

import json
from pathlib import Path

import polars as pl

from ctx.dataset import DatasetError
from ctx.de.results import newest
from ctx.gather import rnaseq
from ctx.io.naming import MODULES
from ctx.process import cache
from ctx.qc.reads import trimmed_names
from ctx.quant import CountsFrom
from ctx.quant.tximport import GeneTables

TRIM_STEP = "trim"
Reads = dict[str, tuple[Path, ...]]


def _process(root: Path) -> Path:
    """M3's output folder."""
    return root / "outputs" / f"{dict(MODULES)['process']:02d}_process"


def trim_folder(root: Path, dataset: str) -> tuple[Path | None, list[str]]:
    """The trim folder M3's newest ``folders`` table names for the dataset, or why not."""
    try:
        table = pl.read_parquet(newest(_process(root), dataset, "folders"))
    except DatasetError as error:
        return None, [f"no M3 run to read trimmed reads from; run ctx process first: {error}"]
    rows = table.filter(pl.col("step") == TRIM_STEP)["folder"].to_list()
    if len(rows) != 1:
        return None, [f"M3's folders table for {dataset} names {len(rows)} trim folders, not 1"]
    folder = _process(root) / "tools" / dataset / rows[0]
    if not (folder / cache.MANIFEST).is_file():
        return None, [f"{folder}: M3's trim folder is missing or unfinished"]
    return folder, []


def trimmed_reads(root: Path, dataset: str, sheet: pl.DataFrame) -> tuple[Reads, Path | None,
                                                                         list[str]]:  # fmt: skip
    """Each sample's trimmed single-end reads, its runs in sheet order, or every problem.

    Args:
        root: The ctx root.
        dataset: The dataset's slug.
        sheet: The gather sheet, one row per run, with ``run_accession`` and ``sample``.

    Returns:
        Sample to read files, the trim folder, and the problems found.
    """
    folder, problems = trim_folder(root, dataset)
    if folder is None:
        return {}, None, problems
    listed = json.loads((folder / cache.MANIFEST).read_text())["files"]
    reads: Reads = {}
    paired_end = pl.col("fastq_2").fill_null("") != ""
    for run, sample, paired in sheet.select("run_accession", "sample", paired_end).rows():
        if paired:
            problems.append(f"{run} is paired-end; the laptop path quantifies single-end runs")
            continue
        name = f"reads/{trimmed_names(run, 1)[0]}"
        path = folder / name
        if name not in listed or not path.is_file() or path.stat().st_size != listed[name]:
            problems.append(f"{path} is not the file M3's manifest lists")
            continue
        reads[sample] = (*reads.get(sample, ()), path)
    return reads, folder, problems


def counts(port: GeneTables, names: dict[str, str], provenance: dict[str, str]) -> rnaseq.Counts:
    """The port's length-scaled counts as gather's integer matrix: what M4's figures read."""
    scaled = port.counts_from(CountsFrom.LENGTH_SCALED_TPM)
    integers, fractional = rnaseq.to_counts(scaled)
    matrix = pl.DataFrame(integers, schema=list(port.samples), orient="row")
    genes = pl.Series("gene_id", list(port.genes))
    named = pl.DataFrame({"gene_id": genes, "gene_name": genes.replace(names)})
    return rnaseq.Counts(matrix.insert_column(0, genes), named, fractional, provenance)
