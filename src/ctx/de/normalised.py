"""M5's normalised expression: each contrast's size factors, and the one normalised value.

M5 writes the size factors of every contrast's fit beside its results (E.1b). A reader that
draws a gene's per-sample values, M5's named-gene boxes and M9's box and F1, reads that table
and calls :func:`log2_normalised`; it never fits size factors of its own.
"""

from collections.abc import Mapping, Sequence

import numpy as np
import polars as pl

from ctx.de.deseq import Fit

WHAT = "size_factors"
COLUMNS = {
    "contrast": ("none", "Contrast name from the dataset file."),
    "sample": ("none", "Sample in the contrast's fit."),
    "donor": ("none", "Individual."),
    "group": ("none", "reference or tested, as the fit coded the compared levels."),
    "size_factor": (
        "dimensionless",
        "DESeq2 median-of-ratios size factor over the contrast's fitted genes; a count over "
        "it is the normalised count.",
    ),
}


def size_factors(fits: Mapping[str, Fit]) -> pl.DataFrame:
    """Every contrast's samples with their size factors, in each fit's sample order.

    Args:
        fits: The fitted contrasts by name.

    Returns:
        One row per contrast and sample, the columns of ``COLUMNS``.
    """
    return pl.concat(
        [
            fit.samples.select("sample", "donor", "group")
            .join(fit.size_factors, on="sample", how="left", maintain_order="left")
            .select(pl.lit(name).alias("contrast"), "sample", "donor", "group", "size_factor")
            for name, fit in fits.items()
        ]
    )


def log2_normalised(
    counts: pl.DataFrame, factors: pl.DataFrame, gene_ids: Sequence[str]
) -> pl.DataFrame:
    """log2(count / size factor + 1) of named genes in the samples of one contrast.

    Args:
        counts: ``gene_id`` then one integer column per sample, as gather wrote it.
        factors: One contrast's ``sample`` and ``size_factor`` rows of the size-factor table.
        gene_ids: The genes; one absent from ``counts`` has no rows.

    Returns:
        ``gene_id``, ``sample`` and ``value``: sample by sample in ``factors``' order, each
        sample's genes in ``counts``' order.
    """
    samples = factors["sample"].to_list()
    long = (
        counts.filter(pl.col("gene_id").is_in(list(gene_ids)))
        .select("gene_id", *samples)
        .unpivot(index="gene_id", variable_name="sample", value_name="count")
        .join(factors.select("sample", "size_factor"), on="sample", maintain_order="left")
    )
    level = long["count"].to_numpy().astype(np.float64) / long["size_factor"].to_numpy()
    return long.select("gene_id", "sample").with_columns(value=pl.Series(np.log2(level + 1)))
