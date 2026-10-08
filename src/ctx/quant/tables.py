"""M4's numbers, each computed once from what Salmon and tximport wrote; figures draw these.

Salmon's fates of a sample's fragments are disjoint and add up to what it processed (Salmon
documentation, ``meta_info.json``): mapped, mapped better to a decoy, aligned below the score
threshold (``num_fragments_filtered_vm``), and no hit at all, the remainder. Detection, the
sample correlation and the units figure read gene tables; the isoform figure reads transcript
TPM. Every random draw is seeded with ``SEED``.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import polars as pl
from numpy.typing import NDArray
from scipy.stats import rankdata

SEED = 42
# A gene is detected at this many reads; a transcript is expressed at this TPM.
MIN_READS = 10
MIN_TPM = 1.0
PERMUTATIONS = 9_999
# Depths each detection curve is drawn at: shares of the sample's own reads, then the common.
DEPTH_SHARES = tuple(float(x) for x in np.geomspace(0.02, 1.0, 15))
FATES = ("mapped", "decoy", "below_score", "no_hit")
# Salmon's MultiQC columns behind each fate; no_hit is what is left of num_processed.
FATE_COLUMNS = {"mapped": "num_mapped", "decoy": "num_decoy_fragments",
                "below_score": "num_fragments_filtered_vm"}  # fmt: skip
# Annotated isoforms per gene, binned for the isoform bars: (low, high), high inclusive.
ISOFORM_BINS = ((2, 2), (3, 4), (5, 8), (9, 16), (17, 1_000_000))
GENE_TYPES = ("protein_coding", "lncRNA")

Matrix = NDArray[np.float64]


def stat(series: pl.Series, how: str) -> float:
    """A column's min, max or median as a float; Polars types the result loosely."""
    value: Any = getattr(series, how)()
    return float(value)


def fates(salmon: pl.DataFrame) -> pl.DataFrame:
    """Per sample: processed fragments and each fate's count; no_hit is the remainder."""
    return salmon.select(
        "sample", processed=pl.col("num_processed"),
        **{fate: pl.col(column) for fate, column in FATE_COLUMNS.items()},
    ).with_columns(no_hit=pl.col("processed") - pl.sum_horizontal(list(FATE_COLUMNS)))  # fmt: skip


def fates_long(table: pl.DataFrame, groups: pl.DataFrame) -> pl.DataFrame:
    """One row per sample and fate: reads, and percent of the fragments processed."""
    return (table.unpivot(list(FATES), index=["sample", "processed"], variable_name="fate",
                          value_name="reads")
            .with_columns(percent=100 * pl.col("reads") / pl.col("processed"))
            .join(groups, on="sample", maintain_order="left"))  # fmt: skip


def library(salmon: pl.DataFrame) -> pl.DataFrame:
    """Per sample: the library type Salmon detected, its strand bias and its run's settings."""
    return salmon.select(
        "sample",
        library_type=pl.col("library_types").str.replace_all(r"[\[\]' ]", ""),
        strand_bias=pl.col("strand_bias"),
        fragment_mean=pl.col("frag_length_mean").cast(pl.Float64),
        fragment_sd=pl.col("frag_length_sd").cast(pl.Float64),
        fragment_source=pl.col("frag_length_source"),
        gc_bias=pl.col("gc_bias_correct"), seq_bias=pl.col("seq_bias_correct"),
        optimiser=pl.col("opt_type"), salmon_version=pl.col("salmon_version"),
        targets=pl.col("num_valid_targets").cast(pl.Int64),
        decoys=pl.col("num_decoy_targets").cast(pl.Int64),
    )  # fmt: skip


def sums(table: pl.DataFrame, scaled: Mapping[str, float]) -> pl.DataFrame:
    """Per sample: Salmon's mapped fragments beside the length-scaled matrix's column sum."""
    return table.select("sample", "mapped").with_columns(
        length_scaled_sum=pl.col("sample").replace_strict(scaled, return_dtype=pl.Float64),
    ).with_columns(relative_difference=(pl.col("length_scaled_sum") - pl.col("mapped"))
                   / pl.col("mapped"))  # fmt: skip


def spearman(matrix: Matrix) -> tuple[Matrix, int]:
    """Spearman's rho between samples over the genes detected on average, and how many."""
    detected = matrix[matrix.mean(axis=1) >= MIN_READS]
    ranks = rankdata(detected, axis=0)
    rho: Matrix = np.corrcoef(ranks, rowvar=False)
    return rho, int(detected.shape[0])


@dataclass(frozen=True)
class Permutation:
    """Mean rho within label groups minus between, and its permutation p.

    Attributes:
        difference: The observed statistic.
        p: (1 + permutations at least as large) / (1 + permutations); NaN when the strata
            leave no label free to move.
        within: Median rho of pairs that share a label.
        between: Median rho of pairs that do not.
        permutations: Label permutations drawn.
        strata: Strata with two labels or more, where labels could move.
    """

    difference: float
    p: float
    within: float
    between: float
    permutations: int
    strata: int


def within_minus_between(rho: Matrix, labels: Sequence[str], strata: Sequence[str],
                         permutations: int = PERMUTATIONS) -> Permutation | None:  # fmt: skip
    """Mean rho within labels minus between, labels permuted only among one stratum's samples.

    Samples of one donor are not exchangeable with another donor's: a free permutation mixes
    the donor effect into the null and calls it a cell-type effect. Permuting within donors
    keeps each donor's samples together, as ``ctx.viz.annotate.permanova`` does.

    Returns:
        The test, or None with one label: no pair lies between labels.
    """
    levels, codes = np.unique(np.asarray(labels), return_inverse=True)
    if levels.size < 2:  # noqa: PLR2004 - within and between need two labels
        return None
    upper = np.triu_indices(codes.size, 1)
    values = rho[upper]

    def statistic(shuffled: NDArray[np.int64]) -> float:
        """The statistic for one labelling."""
        same = shuffled[upper[0]] == shuffled[upper[1]]
        return float(values[same].mean() - values[~same].mean())

    blocks = [np.flatnonzero(np.asarray(strata) == s) for s in dict.fromkeys(strata)]
    movable = [b for b in blocks if np.unique(codes[b]).size > 1]
    observed = statistic(codes)
    same = codes[upper[0]] == codes[upper[1]]
    within, between = float(np.median(values[same])), float(np.median(values[~same]))
    if not movable:
        return Permutation(observed, float("nan"), within, between, 0, 0)
    rng, at_least = np.random.default_rng(SEED), 0
    for _ in range(permutations):
        shuffled = codes.copy()
        for block in movable:
            shuffled[block] = codes[rng.permutation(block)]
        at_least += statistic(shuffled) >= observed - 1e-12 * abs(observed)
    p = (1 + at_least) / (1 + permutations)
    return Permutation(observed, p, within, between, permutations, len(movable))


def thinned(counts: Matrix, samples: Sequence[str]) -> tuple[pl.DataFrame, float]:
    """Genes detected as each sample's reads are thinned binomially, and the common depth.

    Each sample is thinned to shares of its own reads and to the depth every sample reaches,
    the smallest library: a gene's reads are drawn Binomial(reads, depth / total), counts
    rounded half to even first, as tximport's estimates are fractional.

    Returns:
        sample, depth (reads kept), genes (detected at ``MIN_READS``), common (whether the
        point is the common depth); and that depth.
    """
    whole = np.rint(counts).astype(np.int64)
    totals = whole.sum(axis=0)
    common = float(totals.min())
    rng, rows = np.random.default_rng(SEED), []
    for j, sample in enumerate(samples):
        depths = np.union1d(np.round(np.asarray(DEPTH_SHARES) * totals[j]), [common])
        for depth in depths:
            kept = (
                rng.binomial(whole[:, j], depth / totals[j]) if depth < totals[j] else whole[:, j]
            )
            rows.append((sample, float(depth), int((kept >= MIN_READS).sum()), depth == common))
    frame = pl.DataFrame(rows, schema={"sample": pl.String, "depth": pl.Float64,
                                       "genes": pl.Int64, "common": pl.Boolean},
                         orient="row")  # fmt: skip
    return frame, common


def size_factors(counts: Matrix) -> Matrix:
    """DESeq2's median-of-ratios size factors, over genes with a count in every sample."""
    positive = counts[(counts > 0).all(axis=1)]
    log_mean = np.log(positive).mean(axis=1, keepdims=True)
    factors: Matrix = np.exp(np.median(np.log(positive) - log_mean, axis=0))
    return factors


def fold_change(counts: Matrix, tested: NDArray[np.bool_]) -> Matrix:
    """log2 of the tested group's mean normalised count over the reference's, each plus 1."""
    normal = counts / size_factors(counts)
    change: Matrix = np.log2(normal[:, tested].mean(axis=1) + 1) - np.log2(
        normal[:, ~tested].mean(axis=1) + 1)  # fmt: skip
    return change


def isoforms(expressed: pl.DataFrame, tx2gene: pl.DataFrame, types: pl.DataFrame) -> pl.DataFrame:
    """Per gene: transcripts in the index, those expressed, their share, and the gene's type.

    Args:
        expressed: tx and expressed (bool) per transcript.
        tx2gene: transcript_id, gene_id.
        types: gene_id, gene_type.

    Returns:
        gene_id, gene_type, annotated, expressed, share; genes with two transcripts or more.
    """
    joined = tx2gene.select("transcript_id", "gene_id").join(
        expressed.rename({"tx": "transcript_id"}), on="transcript_id", how="left")  # fmt: skip
    per_gene = joined.group_by("gene_id").agg(
        annotated=pl.len(), expressed=pl.col("expressed").fill_null(value=False).sum()
    )
    return (per_gene.filter(pl.col("annotated") >= ISOFORM_BINS[0][0])
            .join(types.select("gene_id", "gene_type"), on="gene_id", how="left")
            .with_columns(share=pl.col("expressed") / pl.col("annotated"))
            .sort("gene_id"))  # fmt: skip


def isoform_bin(annotated: int) -> int:
    """The index of the bin of ``ISOFORM_BINS`` holding a gene's annotated isoforms."""
    return next(i for i, (low, high) in enumerate(ISOFORM_BINS) if low <= annotated <= high)
