"""The tables M2 writes so ``ctx figure`` can redraw its figures from tables alone (F.1).

Each holds what a figure drew and no other table holds: the run's options, as ``ctx qc`` was
given them; the blind VST's PCA scores and variance shares; FastQC's plot data, as MultiQC
wrote it. They are written with the others, before any figure, so a figure's spec finds them.
"""

from typing import TYPE_CHECKING

import polars as pl

from ctx.qc import checkpoint, fastqc
from ctx.qc import registry as reg

if TYPE_CHECKING:
    from ctx.qc.run import Options

# A list option as one cell: the sheet columns of --strata or --peers, in order.
SEPARATOR = ","


def options_row(registry: reg.Registry, options: "Options") -> pl.DataFrame:
    """The run's options a figure reads: groupings, k and the registry."""
    return pl.DataFrame({
        "strata": [SEPARATOR.join(options.strata)], "peers": [SEPARATOR.join(options.peers)],
        "mad": [float(registry.mad_multiplier)], "registry": [str(registry.path)],
    })  # fmt: skip


def pca_tables(measured: checkpoint.Measured) -> dict[str, pl.DataFrame]:
    """The PCA scores per sample and the two variance shares the axes print."""
    share = pl.DataFrame({"component": ["pc1", "pc2"], "share": list(measured.share)})
    return {"pca_scores": measured.scores.select("sample", "pc1", "pc2"), "pca_share": share}


def plot_tables(plots: fastqc.Plots) -> dict[str, pl.DataFrame]:
    """FastQC's curves, one row per metric, unit and point; the GC curve; the top shares.

    The theoretical GC table is written empty where MultiQC wrote no such curve, which
    ``plots_from`` reads back as none.
    """
    key = plots.key
    curves = [f.select(pl.lit(m).alias("metric"), key, "x", "y") for m, f in plots.curves.items()]
    empty = pl.DataFrame(schema={"metric": pl.String, key: pl.String, "x": pl.Float64,
                                 "y": pl.Float64})  # fmt: skip
    theoretical = plots.theoretical
    if theoretical is None:
        theoretical = pl.DataFrame(schema={"x": pl.Float64, "y": pl.Float64})
    return {"fastqc_curves": pl.concat([empty, *curves]), "fastqc_theoretical": theoretical,
            "fastqc_overrepresented": plots.overrepresented}  # fmt: skip


def inputs(registry: reg.Registry, options: "Options", measured: checkpoint.Measured | None,
           plots: fastqc.Plots | None) -> dict[str, pl.DataFrame]:  # fmt: skip
    """Every table above this run can write: the options always, the rest where measured."""
    found = {"figure_options": options_row(registry, options)}
    found |= pca_tables(measured) if measured is not None else {}
    return found | (plot_tables(plots) if plots is not None else {})


def plots_from(curves: pl.DataFrame, theoretical: pl.DataFrame,
               overrepresented: pl.DataFrame) -> fastqc.Plots:  # fmt: skip
    """The plots ``plot_tables`` wrote, read back as ``fastqc.read`` returned them."""
    key = curves.columns[1]
    found = {m: curves.filter(pl.col("metric") == m).drop("metric")
             for m in curves["metric"].unique(maintain_order=True)}  # fmt: skip
    return fastqc.Plots(key, found, None if theoretical.is_empty() else theoretical,
                        overrepresented)  # fmt: skip
