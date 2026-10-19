"""The tables M10 writes so ``ctx figure`` can redraw its figures from tables alone (F.1).

The run's options as its figures read them, and the stored values of the genes the page
carries, the build list's among them: the box, the score, its scan and power read the build
genes, and the live survival panel reads every carried gene. Both are written before any
figure, so a figure's spec finds them; ``read_options`` and ``read_values`` read them back.
"""

from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import polars as pl

from ctx.clinical import ClinicalError, figures
from ctx.report import text

if TYPE_CHECKING:
    from ctx.clinical.page import Context
    from ctx.de.run import Writer

PAGE = "clinical"
# A list option as one cell: the patient-table columns of --covariates or --cox-strata.
SEPARATOR = ","
GENE_COLUMNS = ("gene_id", "gene_name", "build")


def options_row(ctx: "Context") -> pl.DataFrame:
    """The run's options a figure reads: the endpoint as resolved, the split, the models."""
    o = ctx.options
    return pl.DataFrame({
        "endpoint": [ctx.endpoint], "split": [o.split], "horizon": [float(o.horizon)],
        "covariates": [SEPARATOR.join(o.covariates)], "strata": [SEPARATOR.join(o.strata)],
        "source": [ctx.source],
    })  # fmt: skip


def values_table(ctx: "Context") -> pl.DataFrame:
    """The carried genes' stored values, one row per gene in carried order, a column per sample.

    The samples are the patient table's, in its order, as the values' columns are.
    """
    rows = [int(r) for r in ctx.carried.rows]
    ids, symbols = ctx.gene_table["gene_id"].to_list(), ctx.gene_table["gene_name"].to_list()
    build = set(ctx.chosen)
    stored = ctx.values[rows]
    genes = pl.DataFrame({"gene_id": [ids[r] for r in rows],
                          "gene_name": [symbols[r] for r in rows],
                          "build": [r in build for r in rows]})  # fmt: skip
    samples = ctx.clinical.patients["sample"].to_list()
    columns = pl.DataFrame({s: stored[:, i] for i, s in enumerate(samples)},
                           schema=dict.fromkeys(samples, pl.Int16))  # fmt: skip
    return genes.hstack(columns)


def write(writer: "Writer", ctx: "Context") -> None:
    """Write both tables, each with its README."""
    writer.table(options_row(ctx), "figure_options", figures.columns("columns.figure_options"))
    frame = values_table(ctx)
    unit, meaning = text.raw(PAGE, "panel.sample_column")
    samples = dict.fromkeys(frame.columns[len(GENE_COLUMNS) :], (unit, meaning))
    writer.table(frame, "panel_values", {**figures.columns("columns.panel_values"), **samples})


def read_options(path: Path) -> dict[str, Any]:
    """The options row read back, its lists as tuples."""
    row = pl.read_parquet(path).row(0, named=True)
    split = {k: tuple(c for c in row[k].split(SEPARATOR) if c) for k in ("covariates", "strata")}
    return {**row, **split}


def read_values(path: Path, samples: list[str]) -> tuple[pl.DataFrame, np.ndarray, list[int]]:
    """The carried genes, their stored values (genes by patients) and the build list's rows.

    Raises:
        ClinicalError: If the values' samples are not the patient table's, in its order.
    """
    frame = pl.read_parquet(path)
    if (found := frame.columns[len(GENE_COLUMNS) :]) != samples:
        raise ClinicalError(f"panel_values: {len(found)} sample columns, not the patient "
                            f"table's {len(samples)} in its order")  # fmt: skip
    values = frame.select(samples).to_numpy().astype(np.int16)
    chosen = [int(i) for i in np.flatnonzero(frame["build"].to_numpy())]
    return frame.select("gene_id", "gene_name"), values, chosen
