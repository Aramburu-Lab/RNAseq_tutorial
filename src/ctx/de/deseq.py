"""M5: one contrast through pydeseq2, and the blind-VST PCA.

The only module that imports pandas (D-007): pydeseq2 takes and returns pandas objects, so
frames cross into pandas here and come back as polars.

``fit()`` runs in a fixed, tested order: subset the samples; keep one sample per donor,
cell type and condition; drop unpaired donors; recode the compared levels; pre-filter on the
subset; check the group sizes and the design; fit; copy the raw results; shrink separately.
"""

import re
import warnings
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
import polars as pl

from ctx.dataset import Contrast, Split
from ctx.de.filters import MIN_COUNT, expressed
from ctx.de.stats import ADJUSTMENTS, ALPHA, LFC_MIN, adjust, classify
from ctx.gather.curation import LIBRARY, REPLICATE_KEYS, library_replicates, rank_replicates

TOP_GENES = 500
# Genes per block of the initial dispersions (lean_inference): 33 MB a temporary at 1,016
# samples.
GENE_BLOCK = 4096
SPLIT_LEVELS = ("a", "b")
# The coefficient pydeseq2 names for the tested level once "group" is recoded.
COEFFICIENT = "group[T.tested]"
# The shrunk estimate is added unless ``ctx de --no-lfc-shrink`` turns it off.
LFC_SHRINK = True
RESULT_COLUMNS = {
    "gene_id": ("none", "Gene identifier."),
    "baseMean": ("normalised reads", "Mean of size-factor-normalised counts over the samples."),
    "log2FoldChange": ("log2 fold change", "Raw (unshrunk) LFC, tested over reference."),
    "lfcSE": ("log2 fold change", "Standard error of the raw LFC."),
    "stat": ("dimensionless", "Wald statistic, raw LFC over its SE; the GSEA ranking."),
    "pvalue": ("dimensionless", "Wald test p-value; missing for a Cook's outlier."),
    "padj": ("dimensionless", "Benjamini-Hochberg adjusted p; missing if not tested."),
    **{
        column: ("dimensionless", f"{name} adjusted p over the genes with a padj, else missing.")
        for column, name in ADJUSTMENTS.items()
        if column != "padj"
    },
    "log2FoldChange_shrunk": (
        "log2 fold change",
        "LFC after apeGLM-style shrinkage; empty with --no-lfc-shrink.",
    ),
    "lfcSE_shrunk": (
        "log2 fold change",
        "Posterior SE of the shrunk LFC; empty with --no-lfc-shrink.",
    ),
    "na_reason": ("none", "tested, cooks_outlier (no p-value) or independent_filter (no padj)."),
    "call": (
        "none",
        f"padj < {ALPHA} and raw LFC >= {LFC_MIN} is up, <= -{LFC_MIN} down; else ns.",
    ),
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
        trend: The dispersion trend the VST used: ``parametric``, or ``mean`` where
            pydeseq2's parametric fit did not converge and it fell back.
    """

    scores: pl.DataFrame
    share: tuple[float, float]
    genes: tuple[str, ...]
    coords: np.ndarray
    trend: str = "parametric"


@dataclass(frozen=True)
class Fit:
    """A fitted contrast."""

    contrast: Contrast
    results: pl.DataFrame
    samples: pl.DataFrame
    dropped: tuple[str, ...]
    dispersions: pl.DataFrame
    size_factors: pl.DataFrame


def _adjusted(results: pl.DataFrame) -> pl.DataFrame:
    """Add every other adjustment, over exactly the genes pydeseq2's BH adjusted.

    The same hypotheses as padj, after Cook's and independent filtering, so the columns
    differ by method alone; a gene with no padj has none in any column.
    """
    tested = results["padj"].is_not_null().to_numpy()
    p = results["pvalue"].to_numpy()[tested]
    columns = []
    for method in (m for m in ADJUSTMENTS if m != "padj"):
        values = np.full(results.height, np.nan)
        values[tested] = adjust(p, method)
        columns.append(pl.Series(method, values).fill_nan(None))
    return results.with_columns(columns)


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
    first = _first_libraries(subset, order)
    dropped = subset.filter(~first)["run_accession"].to_list()
    # A kept library's other runs are its lanes, part of the one sample, never dropped (D-063);
    # runs of one sample that differ in a key are a library per value, as curation ranks them.
    lanes = ["sample", *(k for k in (*REPLICATE_KEYS, *order) if k in subset.columns)]
    used = subset.filter(first).unique(lanes, keep="first", maintain_order=True)
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


def _first_libraries(subset: pl.DataFrame, order: Mapping[str, tuple[str, ...]]) -> pl.Series:
    """Whether each run belongs to the first library of its donor, cell type and condition.

    Curation's rule, ranked again over the contrast's subset: a sheet's sample is one library,
    so every run of it, a lane, takes the library's rank (``library_replicates``, D-063). A GDC
    sheet has no cell type: its one file per patient and sample type is already chosen, and
    its runs are ranked as they are.
    """
    if not set(REPLICATE_KEYS) <= set(subset.columns):
        keys = tuple(k for k in REPLICATE_KEYS if k in subset.columns)
        return rank_replicates(subset, keys, dict(order)) == 1
    libraries = subset.with_columns(pl.col("sample").alias(LIBRARY))
    return library_replicates(libraries, dict(order)) == 1


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


def samples_by_genes(counts: pl.DataFrame, samples: list[str], rows: np.ndarray) -> Any:
    """The counts pydeseq2 fits, samples by genes, as the one int64 copy it keeps (S.5).

    pydeseq2 casts its counts with ``astype(int)``, which under pandas' copy-on-write passes an
    int64 frame through uncopied, so this array becomes its ``X``. Built from polars in one
    step, it replaces three copies of the matrix: polars' int32 export, its row-filtered copy
    and pydeseq2's int64 cast, 953 MB against 471 MB for 1,016 samples by 58,038 genes. Its
    layout is the one those copies gave, Fortran order, so every reduction rounds as before.

    Args:
        counts: ``gene_id`` then one integer column per sample.
        samples: The samples, in the order pydeseq2 gets them.
        rows: Which genes, a boolean per row of ``counts``.

    Returns:
        A pandas frame, samples by genes, int64 in Fortran order.
    """
    import pandas as pd  # noqa: PLC0415 - heavy; the pydeseq2 boundary

    kept = counts.filter(pl.Series(rows))
    values = kept.select(samples).to_numpy(order="c").astype(np.int64).T
    return pd.DataFrame(values, index=samples, columns=kept["gene_id"].to_list(), copy=False)


def _shrunk_gene(shrink_index: int, **gene: Any) -> tuple[np.ndarray, float, bool]:
    """One gene's apeGLM fit, keeping of its inverse Hessian the one entry DeseqStats reads."""
    from pydeseq2.utils import nbinomGLM  # noqa: PLC0415 - heavy; runs in a worker

    beta, inverse_hessian, converged = nbinomGLM(shrink_index=shrink_index, **gene)
    return beta, float(inverse_hessian[shrink_index, shrink_index]), converged


def lean_inference(threads: int) -> Any:
    """pydeseq2's default inference, the same numbers from less memory (S.5).

    Two of its steps hold far more than they keep (pydeseq2 0.5.4). ``DeseqStats.lfc_shrink``
    collects each gene's full inverse Hessian, p by p for a design of p columns, twice over,
    to read one entry, the shrunk coefficient's variance: for 113 donors paired, p = 114,
    104 kB per gene, 4.9 GB traced over some 22,000 genes. Here each worker returns that entry,
    and a view repeating it at every index of a p-by-p matrix stands in for the matrices. The
    rough and moments dispersions take several samples-by-genes temporaries at once, 2.36 GB
    above entry for 1,016 samples by 58,038 genes; each gene's estimate reads only its own
    column, so they run ``GENE_BLOCK`` genes at a time, equal to the bit on that matrix.

    Args:
        threads: CPUs for the inference's workers.

    Returns:
        A ``DefaultInference``.
    """
    # joblib through pydeseq2's module, so the threads test's spy on its Parallel sees this loop.
    from pydeseq2 import default_inference as pydeseq2_inference  # noqa: PLC0415 - heavy
    from pydeseq2 import utils  # noqa: PLC0415 - heavy

    def blocks(estimate: Any, normed_counts: np.ndarray, other: Any) -> np.ndarray:
        """One per-gene estimate, ``GENE_BLOCK`` genes at a time."""
        return np.concatenate([
            estimate(normed_counts[:, i : i + GENE_BLOCK], other)
            for i in range(0, normed_counts.shape[1], GENE_BLOCK)
        ])  # fmt: skip

    class LeanInference(pydeseq2_inference.DefaultInference):  # type: ignore[misc]
        """DefaultInference with the shrinkage and initial dispersions above."""

        def fit_rough_dispersions(self, normed_counts: np.ndarray, design: Any) -> np.ndarray:
            """pydeseq2's rough dispersions, in blocks of genes."""
            return blocks(utils.fit_rough_dispersions, normed_counts, design)

        def fit_moments_dispersions(self, normed_counts: np.ndarray, factors: Any) -> np.ndarray:
            """pydeseq2's moments dispersions, in blocks of genes."""
            return blocks(utils.fit_moments_dispersions, normed_counts, factors)

        def lfc_shrink_nbinom_glm(  # noqa: PLR0913, PLR0917 - pydeseq2's signature
            self,
            design_matrix: np.ndarray,
            counts: np.ndarray,
            size: np.ndarray,
            offset: np.ndarray,
            prior_no_shrink_scale: float,
            prior_scale: float,
            optimizer: str,
            shrink_index: int,
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            """pydeseq2's loop over genes, each returning one variance, not a matrix."""
            with pydeseq2_inference.parallel_backend(self._backend, inner_max_num_threads=1):
                fits = pydeseq2_inference.Parallel(
                    n_jobs=self.n_cpus, verbose=self._joblib_verbosity,
                    batch_size=self._batch_size,
                )(
                    pydeseq2_inference.delayed(_shrunk_gene)(
                        shrink_index, design_matrix=design_matrix, counts=counts[:, i],
                        size=size[i], offset=offset, prior_no_shrink_scale=prior_no_shrink_scale,
                        prior_scale=prior_scale, optimizer=optimizer,
                    )
                    for i in range(counts.shape[1])
                )  # fmt: skip
            lfcs, variances, converged = (np.array(m) for m in zip(*fits, strict=True))
            p = design_matrix.shape[1]
            stand_in = np.lib.stride_tricks.as_strided(
                variances, shape=(len(variances), p, p), strides=(variances.strides[0], 0, 0),
                writeable=False,
            )  # fmt: skip
            return lfcs, stand_in, converged

    return LeanInference(n_cpus=threads)


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
    shrink: bool = LFC_SHRINK,
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
        shrink: Add apeGLM-shrunk fold changes beside the raw ones. Off skips that fit and
            leaves the two shrunk columns empty; nothing tested or called depends on them.

    Returns:
        The fit: one row per tested gene, with raw and shrunk LFC and why a value is missing.

    Raises:
        ContrastError: If the group sizes differ from the declared ones or the design is
            rank deficient.
    """
    import pandas as pd  # noqa: PLC0415 - heavy; the pydeseq2 boundary
    from pydeseq2.dds import DeseqDataSet  # noqa: PLC0415 - heavy
    from pydeseq2.ds import DeseqStats  # noqa: PLC0415 - heavy

    used, dropped = select_samples(sheet, contrast, order)
    n = (used.filter(pl.col("group") == "reference").height,
         used.filter(pl.col("group") == "tested").height)  # fmt: skip
    if n != contrast.expected_n:
        raise ContrastError(f"{contrast.name}: n = {n}, declared {contrast.expected_n}")
    names = used["sample"].to_list()
    keep = expressed(counts.select(names).to_numpy(), list(n), min_count)
    genes = counts["gene_id"].filter(pl.Series(keep)).to_list()
    metadata = pd.DataFrame(
        {c: used[c].to_list() for c in ["group", *covariates(contrast, used)]}, index=names
    ).astype("category")
    _check_design(contrast, metadata)
    inference = lean_inference(threads)
    dds = DeseqDataSet(
        counts=samples_by_genes(counts, names, keep), metadata=metadata,
        design=contrast.design, refit_cooks=True, inference=inference, quiet=True,
    )  # fmt: skip
    dds.deseq2()
    wald = ["group", "tested", "reference"]
    raw = DeseqStats(dds, contrast=wald, alpha=ALPHA, inference=inference, quiet=True)
    raw.summary()
    results = raw.results_df.copy()
    empty = pl.lit(None, pl.Float64)
    shrunk: dict[str, pl.Expr | pl.Series] = {"log2FoldChange_shrunk": empty, "lfcSE_shrunk": empty}
    if shrink:
        # Shrinkage overwrites results_df and needs the Wald SE, so it gets its own object.
        stats = DeseqStats(dds, contrast=wald, alpha=ALPHA, inference=inference, quiet=True)
        stats.summary()
        stats.lfc_shrink(coeff=COEFFICIENT)
        shrunk = {
            "log2FoldChange_shrunk": pl.Series(stats.results_df["log2FoldChange"].to_numpy()),
            "lfcSE_shrunk": pl.Series(stats.results_df["lfcSE"].to_numpy()),
        }
    _shutdown_workers()
    out = pl.from_pandas(results.reset_index(names="gene_id")).with_columns(**shrunk)
    out = out.with_columns(pl.col(pl.Float64).fill_nan(None)).with_columns(
        na_reason=pl.when(pl.col("pvalue").is_null()).then(pl.lit("cooks_outlier"))
        .when(pl.col("padj").is_null()).then(pl.lit("independent_filter"))
        .otherwise(pl.lit("tested")),
        call=classify(pl.col("padj"), pl.col("log2FoldChange")),
    )  # fmt: skip
    out = _adjusted(out)
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

    The VST runs on every gene with a count, decided (P2S.7), not inherited from the fits'
    expression filter: the view is unsupervised and design-blind, and a filter set per
    contrast has no contrast here; the recorded PERMANOVA metrics stand on this set (D-046).
    Where pydeseq2's parametric trend does not converge (Step-1's ss_ndn, 17 samples) it falls
    back to the mean of the dispersions, and ``trend`` says so.

    Args:
        counts: ``gene_id`` then one integer column per sample.
        samples: The samples to include.
        threads: CPUs for pydeseq2.

    Returns:
        The scores, the variance shares of PC1 and PC2, and the genes used.
    """
    import pandas as pd  # noqa: PLC0415 - heavy; the pydeseq2 boundary
    from pydeseq2.dds import DeseqDataSet  # noqa: PLC0415 - heavy

    nonzero = counts.select(pl.any_horizontal(pl.col(samples) > 0)).to_series().to_numpy()
    genes = counts["gene_id"].filter(pl.Series(nonzero)).to_list()
    metadata = pd.DataFrame({"all": ["all"] * len(samples)}, index=samples)
    dds = DeseqDataSet(
        counts=samples_by_genes(counts, samples, nonzero), metadata=metadata, design="~1",
        inference=lean_inference(threads), quiet=True,
    )  # fmt: skip
    with warnings.catch_warnings():
        # The fallback is recorded in ``trend`` and printed in the caption, not warned.
        warnings.filterwarnings("ignore", message="The dispersion trend curve fitting did not")
        # dds.vst(use_design=False) in its own steps, so the fitted means, 471 MB at 1,016
        # samples by 58,038 genes and read by nothing after the fit, go before the transform.
        dds.vst_fit_type = dds.fit_type
        dds.vst_fit(use_design=False)
    _shutdown_workers()
    del dds.layers["_vst_mu_hat"]
    values = np.asarray(dds.vst_transform())
    trend = str(dds.vst_fit_type)
    # The counts and normalised counts go before the variance ranking's samples-by-genes
    # temporary is made, not after.
    del dds
    top = np.argsort(-values.var(axis=0), kind="stable")[:TOP_GENES]
    centred = values[:, top] - values[:, top].mean(axis=0)
    u, s, _ = np.linalg.svd(centred, full_matrices=False)
    share = s**2 / (s**2).sum()
    scores = pl.DataFrame({"sample": samples, "pc1": u[:, 0] * s[0], "pc2": u[:, 1] * s[1]})
    genes_used = tuple(genes[i] for i in top)
    return Pca(scores, (float(share[0]), float(share[1])), genes_used, u * s, trend)
