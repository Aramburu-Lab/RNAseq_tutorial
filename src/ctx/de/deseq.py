"""M5: one contrast through pydeseq2, and the blind-VST PCA.

The only module that imports pandas (D-007): pydeseq2 takes and returns pandas objects, so
frames cross into pandas here and come back as polars.

``fit()`` runs in a fixed, tested order: subset the samples; keep one sample per donor,
cell type and condition; drop unpaired donors; recode the compared levels; pre-filter on the
subset; check the group sizes and the design; fit; copy the raw results; shrink separately.
"""

import re
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
import polars as pl

from ctx.dataset import Contrast, Split
from ctx.de.filters import MIN_COUNT, expressed
from ctx.de.stats import ALPHA, classify
from ctx.gather.curation import REPLICATE_KEYS, rank_replicates

TOP_GENES = 500
SPLIT_LEVELS = ("a", "b")
# The coefficient pydeseq2 names for the tested level once "group" is recoded.
COEFFICIENT = "group[T.tested]"
RESULT_COLUMNS = {
    "gene_id": ("none", "Gene identifier."),
    "baseMean": ("normalised reads", "Mean of size-factor-normalised counts over the samples."),
    "log2FoldChange": ("log2 fold change", "Raw (unshrunk) LFC, tested over reference."),
    "lfcSE": ("log2 fold change", "Standard error of the raw LFC."),
    "stat": ("dimensionless", "Wald statistic, raw LFC over its SE; the GSEA ranking."),
    "pvalue": ("dimensionless", "Wald test p-value; missing for a Cook's outlier."),
    "padj": ("dimensionless", "Benjamini-Hochberg adjusted p; missing if not tested."),
    "log2FoldChange_shrunk": ("log2 fold change", "LFC after apeGLM-style shrinkage."),
    "lfcSE_shrunk": ("log2 fold change", "Posterior SE of the shrunk LFC."),
    "na_reason": ("none", "tested, cooks_outlier (no p-value) or independent_filter (no padj)."),
    "call": ("none", f"up, down or ns at padj < {ALPHA}, signed by the raw LFC."),
}


class ContrastError(ValueError):
    """A contrast that cannot be fitted as declared."""


@dataclass(frozen=True)
class Pca:
    """Principal components of the blind VST.

    Attributes:
        scores: ``sample``, ``pc1`` and ``pc2``.
        share: The share of the selected genes' variance each component explains.
        genes: The genes used: the ``TOP_GENES`` most variable after the VST.
        coords: Every component's scores, samples by components; Euclidean distances here
            equal those between the samples' centred VST values over ``genes``.
    """

    scores: pl.DataFrame
    share: tuple[float, float]
    genes: tuple[str, ...]
    coords: np.ndarray


@dataclass(frozen=True)
class Fit:
    """A fitted contrast."""

    contrast: Contrast
    results: pl.DataFrame
    samples: pl.DataFrame
    dropped: tuple[str, ...]
    dispersions: pl.DataFrame
    size_factors: pl.DataFrame


def _split(frame: pl.DataFrame, split: Split) -> pl.DataFrame:
    """Assign individuals alternately to a and b, shuffled within each stratum by the seed."""
    rng = np.random.default_rng(split.seed)
    donors = frame.select("donor", split.strata).unique("donor").sort(split.strata, "donor")
    labels: dict[str, str] = {}
    for _, group in donors.group_by(split.strata, maintain_order=True):
        for donor in rng.permutation(group["donor"].to_list()):
            labels[donor] = SPLIT_LEVELS[len(labels) % 2]
    return frame.with_columns(split=pl.col("donor").replace_strict(labels))


def select_samples(
    sheet: pl.DataFrame, contrast: Contrast, order: Mapping[str, tuple[str, ...]]
) -> tuple[pl.DataFrame, tuple[str, ...]]:
    """The samples a contrast compares, with ``group``, and the runs it dropped.

    Order: subset by the filters and the two levels; rank replicates on the subset and keep
    the first; split, if declared; drop donors without both levels, if paired; recode.

    Args:
        sheet: The harmonised sheet, curated.
        contrast: The contrast.
        order: The dataset's level orders, for ranking replicates.

    Returns:
        The used samples with a ``group`` column (reference or tested), and the dropped runs.
    """
    subset = sheet
    for column, allowed in contrast.filters.items():
        subset = subset.filter(pl.col(column).is_in(list(allowed)))
    variable, levels = contrast.variable, (contrast.reference, contrast.tested)
    if contrast.split is None:
        subset = subset.filter(pl.col(variable).is_in(list(levels)))
    first = rank_replicates(subset, REPLICATE_KEYS, dict(order)) == 1
    dropped = subset.filter(~first)["run_accession"].to_list()
    used = subset.filter(first)
    if contrast.split is not None:
        used, variable = _split(used, contrast.split), "split"
    if contrast.paired:
        levels_seen = used.group_by("donor").agg(pl.col(variable).n_unique().alias("n"))
        paired = levels_seen.filter(pl.col("n") == len(levels))["donor"].to_list()
        dropped += used.filter(~pl.col("donor").is_in(paired))["run_accession"].to_list()
        used = used.filter(pl.col("donor").is_in(paired))
    tested = SPLIT_LEVELS[1] if contrast.split else contrast.tested
    group = (
        pl.when(pl.col(variable) == tested).then(pl.lit("tested")).otherwise(pl.lit("reference"))
    )
    return used.with_columns(group=group), tuple(sorted(dropped))


def covariates(contrast: Contrast, sheet: pl.DataFrame) -> list[str]:
    """The sheet columns the design names besides ``group``; each must exist."""
    terms = [t for t in dict.fromkeys(re.findall(r"[A-Za-z_]\w*", contrast.design)) if t != "group"]
    if missing := [t for t in terms if t not in sheet.columns]:
        raise ContrastError(f"{contrast.name}: design names {missing}, not sheet columns")
    return terms


def _check_design(contrast: Contrast, metadata: Any) -> None:
    """Refuse a rank-deficient design; pydeseq2 only warns and then fits nonsense."""
    from formulaic import model_matrix  # noqa: PLC0415 - pydeseq2's own formula engine

    matrix = model_matrix(contrast.design, metadata)
    if np.linalg.matrix_rank(np.asarray(matrix, dtype=float)) < matrix.shape[1]:
        raise ContrastError(
            f"{contrast.name}: design {contrast.design} is rank deficient on these samples "
            f"({matrix.shape[1]} columns: {list(matrix.columns)})"
        )


def _shutdown_workers() -> None:
    """Stop loky's reusable workers, so their memory is reaped into RUSAGE_CHILDREN."""
    from joblib.externals.loky import get_reusable_executor  # noqa: PLC0415 - heavy

    get_reusable_executor().shutdown(wait=True)


def fit(  # noqa: PLR0913 - every argument is a separate input of the fit
    counts: pl.DataFrame,
    sheet: pl.DataFrame,
    contrast: Contrast,
    order: Mapping[str, tuple[str, ...]],
    *,
    threads: int,
    min_count: int = MIN_COUNT,
) -> Fit:
    """Fit one contrast with pydeseq2.

    Args:
        counts: ``gene_id`` then one integer column per sample.
        sheet: The harmonised sheet, curated.
        contrast: The contrast.
        order: The dataset's level orders.
        threads: CPUs for every pydeseq2 object: both the dataset and each stats object,
            which otherwise starts its own inference on every CPU.
        min_count: The pre-filter's reads at the median library size.

    Returns:
        The fit: one row per tested gene, with raw and shrunk LFC and why a value is missing.

    Raises:
        ContrastError: If the group sizes differ from the declared ones or the design is
            rank deficient.
    """
    import pandas as pd  # noqa: PLC0415 - heavy; the pydeseq2 boundary
    from pydeseq2.dds import DeseqDataSet  # noqa: PLC0415 - heavy
    from pydeseq2.default_inference import DefaultInference  # noqa: PLC0415 - heavy
    from pydeseq2.ds import DeseqStats  # noqa: PLC0415 - heavy

    used, dropped = select_samples(sheet, contrast, order)
    n = (used.filter(pl.col("group") == "reference").height,
         used.filter(pl.col("group") == "tested").height)  # fmt: skip
    if n != contrast.expected_n:
        raise ContrastError(f"{contrast.name}: n = {n}, declared {contrast.expected_n}")
    names = used["sample"].to_list()
    matrix = counts.select(names).to_numpy()
    keep = expressed(matrix, list(n), min_count)
    genes = counts["gene_id"].filter(pl.Series(keep)).to_list()
    metadata = pd.DataFrame(
        {c: used[c].to_list() for c in ["group", *covariates(contrast, used)]}, index=names
    ).astype("category")
    _check_design(contrast, metadata)
    inference = DefaultInference(n_cpus=threads)
    frame = pd.DataFrame(matrix[keep].T, index=names, columns=genes)
    dds = DeseqDataSet(
        counts=frame, metadata=metadata, design=contrast.design, refit_cooks=True,
        inference=inference, quiet=True,
    )  # fmt: skip
    dds.deseq2()
    wald = ["group", "tested", "reference"]
    raw = DeseqStats(dds, contrast=wald, alpha=ALPHA, inference=inference, quiet=True)
    raw.summary()
    results = raw.results_df.copy()
    # Shrinkage overwrites results_df and needs the Wald SE, so it gets its own object.
    shrink = DeseqStats(dds, contrast=wald, alpha=ALPHA, inference=inference, quiet=True)
    shrink.summary()
    shrink.lfc_shrink(coeff=COEFFICIENT)
    _shutdown_workers()
    out = pl.from_pandas(results.reset_index(names="gene_id")).with_columns(
        log2FoldChange_shrunk=pl.Series(shrink.results_df["log2FoldChange"].to_numpy()),
        lfcSE_shrunk=pl.Series(shrink.results_df["lfcSE"].to_numpy()),
    )
    out = out.with_columns(pl.col(pl.Float64).fill_nan(None)).with_columns(
        na_reason=pl.when(pl.col("pvalue").is_null()).then(pl.lit("cooks_outlier"))
        .when(pl.col("padj").is_null()).then(pl.lit("independent_filter"))
        .otherwise(pl.lit("tested")),
        call=classify(pl.col("padj"), pl.col("log2FoldChange")),
    )  # fmt: skip
    var = dds.var
    dispersions = pl.DataFrame({
        "gene_id": genes,
        "baseMean": results["baseMean"].to_numpy(),
        "genewise": var["genewise_dispersions"].to_numpy(),
        "fitted": var["fitted_dispersions"].to_numpy(),
        "final": var["dispersions"].to_numpy(),
    })  # fmt: skip
    factors = pl.DataFrame({"sample": names, "size_factor": dds.obs["size_factors"].to_numpy()})
    return Fit(contrast, out.select(list(RESULT_COLUMNS)), used, dropped, dispersions, factors)


def pca(counts: pl.DataFrame, samples: list[str], *, threads: int) -> Pca:
    """Principal components of the blind VST, on the most variable genes after the transform.

    Genes are ranked by their variance after the VST, not before: raw counts rank genes by
    depth, and a deep gene with no biology would lead every component.

    Args:
        counts: ``gene_id`` then one integer column per sample.
        samples: The samples to include.
        threads: CPUs for pydeseq2.

    Returns:
        The scores, the variance shares of PC1 and PC2, and the genes used.
    """
    import pandas as pd  # noqa: PLC0415 - heavy; the pydeseq2 boundary
    from pydeseq2.dds import DeseqDataSet  # noqa: PLC0415 - heavy
    from pydeseq2.default_inference import DefaultInference  # noqa: PLC0415 - heavy

    matrix = counts.select(samples).to_numpy()
    nonzero = matrix.sum(axis=1) > 0
    genes = counts["gene_id"].filter(pl.Series(nonzero)).to_list()
    frame = pd.DataFrame(matrix[nonzero].T, index=samples, columns=genes)
    metadata = pd.DataFrame({"all": ["all"] * len(samples)}, index=samples)
    dds = DeseqDataSet(
        counts=frame, metadata=metadata, design="~1",
        inference=DefaultInference(n_cpus=threads), quiet=True,
    )  # fmt: skip
    dds.vst(use_design=False)
    _shutdown_workers()
    values = np.asarray(dds.layers["vst_counts"])
    top = np.argsort(-values.var(axis=0), kind="stable")[:TOP_GENES]
    centred = values[:, top] - values[:, top].mean(axis=0)
    u, s, _ = np.linalg.svd(centred, full_matrices=False)
    share = s**2 / (s**2).sum()
    scores = pl.DataFrame({"sample": samples, "pc1": u[:, 0] * s[0], "pc2": u[:, 1] * s[1]})
    genes_used = tuple(genes[i] for i in top)
    return Pca(scores, (float(share[0]), float(share[1])), genes_used, u * s)
