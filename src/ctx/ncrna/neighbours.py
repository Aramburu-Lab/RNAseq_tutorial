"""Each tested lncRNA and its nearest tested protein-coding gene, and how the two co-vary.

A lncRNA has no GO terms, so its neighbours describe it. Pairs are classed by where the two
genes lie: overlapping on opposite strands (antisense), overlapping on one strand, apart by at
most ``NEAR_NT``, or farther. Spearman's rho is taken on the values and on the residuals after
each stratum's mean (cell type by condition) is removed, since a shared response to the
condition correlates genes that never touch; BH runs over the pairs. Each lncRNA also meets a
protein-coding gene on another chromosome, drawn with a seed: the background a cis pair must
beat. ``ctxPairs`` (``templates/_pairs.js.j2``) is the page's twin of ``spearman`` and
``residuals``, tested under node against them (``tests/test_ncrna.py``).
"""

import numpy as np
import polars as pl
from scipy.stats import rankdata
from scipy.stats import t as student

from ctx.de.stats import adjust

ANTISENSE, SAME_STRAND, NEAR, FAR, BACKGROUND = "antisense", "same_strand", "near", "far", "trans"
KINDS = (ANTISENSE, SAME_STRAND, NEAR, FAR)
NEAR_NT = 10_000
SEED = 8
LOCUS = ("gene_id", "chrom", "start", "end", "strand")


def nearest(lncrnas: pl.DataFrame, coding: pl.DataFrame) -> pl.DataFrame:
    """The nearest protein-coding gene of each lncRNA on its chromosome, and the pair's kind.

    The gap is the nucleotides between the two genes, 0 for genes that touch or overlap;
    overlap decides the kind, not the gap. A tie goes to the gene that starts first. A lncRNA
    on a chromosome with no protein-coding gene has no pair.

    Args:
        lncrnas: ``LOCUS`` columns, 1-based inclusive coordinates.
        coding: The same, of the protein-coding genes.

    Returns:
        lncrna, neighbour, gap and kind, one row per lncRNA with a pair, in input order.
    """
    by_chrom = {chrom: rows.sort("start", "gene_id")
                for (chrom,), rows in coding.group_by("chrom")}  # fmt: skip
    found = []
    for gene, chrom, start, end, strand in lncrnas.select(LOCUS).rows():
        rows = by_chrom.get(chrom)
        if rows is None:
            continue
        between = np.maximum(rows["start"].to_numpy() - end, start - rows["end"].to_numpy()) - 1
        i = int(np.argmin(between))
        overlap = between[i] < 0
        kind = ((ANTISENSE if rows["strand"][i] != strand else SAME_STRAND) if overlap
                else NEAR if between[i] <= NEAR_NT else FAR)  # fmt: skip
        found.append((gene, rows["gene_id"][i], max(int(between[i]), 0), kind))
    schema = {"lncrna": pl.String, "neighbour": pl.String, "gap": pl.Int64, "kind": pl.String}
    return pl.DataFrame(found, schema=schema, orient="row")


def background(pairs: pl.DataFrame, coding: pl.DataFrame, lncrnas: pl.DataFrame) -> pl.DataFrame:
    """Each paired lncRNA with a protein-coding gene drawn, seeded, from another chromosome.

    A lncRNA with no protein-coding gene off its chromosome has no background pair.

    Returns:
        lncrna, neighbour, gap (null: no distance across chromosomes) and kind ``BACKGROUND``.
    """
    rng = np.random.default_rng(SEED)
    chrom = dict(lncrnas.select("gene_id", "chrom").rows())
    ids, chroms = coding["gene_id"].to_numpy(), coding["chrom"].to_numpy()
    drawn = []
    for gene in pairs["lncrna"]:
        elsewhere = ids[chroms != chrom[gene]]
        if elsewhere.size:
            drawn.append((gene, str(rng.choice(elsewhere))))
    return pl.DataFrame(drawn, schema={"lncrna": pl.String, "neighbour": pl.String},
                        orient="row").with_columns(
        gap=pl.lit(None, dtype=pl.Int64), kind=pl.lit(BACKGROUND)
    )  # fmt: skip


def residuals(values: np.ndarray, strata: np.ndarray) -> np.ndarray:
    """Each row minus its mean within each stratum of the columns.

    Args:
        values: Rows by samples.
        strata: One label per sample.
    """
    out = np.array(values, dtype=float)
    for level in np.unique(strata):
        cols = strata == level
        out[:, cols] -= out[:, cols].mean(axis=1, keepdims=True)
    return out


def spearman(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Spearman's rho of each row of ``a`` with the same row of ``b``: ties take mean ranks.

    A row with one value throughout has no rank order, and NaN.
    """
    ra, rb = rankdata(a, axis=1), rankdata(b, axis=1)
    ra -= ra.mean(axis=1, keepdims=True)
    rb -= rb.mean(axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.asarray(
            (ra * rb).sum(axis=1) / np.sqrt((ra**2).sum(axis=1) * (rb**2).sum(axis=1)),
            dtype=float,
        )


def p_value(rho: np.ndarray, n: int) -> np.ndarray:
    """Two-sided p of each rho over ``n`` samples, by the t approximation scipy's spearmanr uses."""
    with np.errstate(divide="ignore", invalid="ignore"):
        t = rho * np.sqrt((n - 2) / ((rho + 1.0) * (1.0 - rho)))
    return np.asarray(2 * student.sf(np.abs(t), n - 2), dtype=float)


def correlate(a: np.ndarray, b: np.ndarray, strata: np.ndarray) -> dict[str, np.ndarray]:
    """rho, p and BH q of each row pair, on values and on residuals within ``strata``.

    Rows whose rho has no value (a constant row) get no p and no q; BH runs over the rest.

    Returns:
        ``rho``, ``p``, ``q``, ``rho_within``, ``p_within``, ``q_within``, one value per row.
    """
    out: dict[str, np.ndarray] = {}
    n = a.shape[1]
    for suffix, (u, v) in (("", (a, b)), ("_within", (residuals(a, strata), residuals(b, strata)))):
        rho = spearman(u, v)
        p = p_value(rho, n)
        q = np.full(rho.size, np.nan)
        held = np.isfinite(p)
        q[held] = adjust(p[held], "padj")
        out |= {f"rho{suffix}": rho, f"p{suffix}": p, f"q{suffix}": q}
    return out
