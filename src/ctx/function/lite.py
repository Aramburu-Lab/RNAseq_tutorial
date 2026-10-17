"""M6 enrichment: over-representation and preranked GSEA, with numpy, scipy and the stdlib only.

Only those, so the same code runs in Pyodide for the JupyterLite notebooks (Phase 5). The
definitions follow GSEA (Subramanian et al. 2005) as gseapy 1.3.1 ``algorithm.py`` writes them,
and gseapy is the test oracle (D-008): ES is the running sum's extreme of larger magnitude;
NES divides by the mean null ES of the same sign; FDR comes from the pooled normalised nulls.
One deliberate difference: the nominal p counts the observed score among the permutations,
so it is never 0 (``(1 + more extreme) / (1 + same-sign nulls)``).
"""

import math
import warnings
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.stats import false_discovery_control, hypergeom

PERMUTATIONS = 1000
SEED = 42
MIN_SIZE = 15
MAX_SIZE = 500


def read_gmt(path: Path) -> dict[str, tuple[str, ...]]:
    """Gene sets from a GMT file: name, description, then genes, tab-separated.

    Args:
        path: The GMT file.

    Returns:
        Set name to its genes, each gene once, in file order.
    """
    sets = {}
    for line in path.read_text().splitlines():
        name, _, *genes = line.split("\t")
        sets[name] = tuple(dict.fromkeys(g for g in genes if g))
    return sets


def sized(
    sets: Mapping[str, Iterable[str]], universe: set[str], min_size: int, max_size: int
) -> dict[str, tuple[str, ...]]:
    """Each set cut to the universe first, then kept if its size is within the limits.

    The order matters: a set of 600 with 400 genes in the universe is tested at 400.
    """
    cut = {name: tuple(g for g in genes if g in universe) for name, genes in sets.items()}
    return {name: genes for name, genes in cut.items() if min_size <= len(genes) <= max_size}


@dataclass(frozen=True)
class Overrepresented:
    """One set's over-representation test."""

    name: str
    size: int
    hits: int
    expected: float
    p: float
    fdr: float
    genes: tuple[str, ...]


def ora(
    hits: Iterable[str],
    background: Iterable[str],
    sets: Mapping[str, Iterable[str]],
    min_size: int = MIN_SIZE,
    max_size: int = MAX_SIZE,
) -> list[Overrepresented]:
    """Hypergeometric over-representation of the hits in each set, against a background.

    The background must be the genes that could have been hits: the tested genes. Hits and
    sets are cut to it; P(X >= hits) is ``hypergeom.sf(hits - 1, N, set size, n hits)``;
    BH runs over every set that passed the size limits.

    Args:
        hits: The called genes.
        background: The tested genes.
        sets: Set name to genes.
        min_size: Smallest set, after the cut to the background.
        max_size: Largest set, after the cut.

    Returns:
        One result per set within the size limits, in the order of ``sets``.
    """
    universe = set(background)
    drawn = set(hits) & universe
    kept = sized(sets, universe, min_size, max_size)
    rows = []
    for name, genes in kept.items():
        overlap = tuple(g for g in genes if g in drawn)
        p = float(hypergeom.sf(len(overlap) - 1, len(universe), len(genes), len(drawn)))
        expected = len(drawn) * len(genes) / len(universe)
        rows.append((name, len(genes), len(overlap), expected, p, overlap))
    fdrs = false_discovery_control([r[4] for r in rows], method="bh") if rows else []
    return [
        Overrepresented(name, size, k, expected, p, float(fdr), overlap)
        for (name, size, k, expected, p, overlap), fdr in zip(rows, fdrs, strict=True)
    ]


@dataclass(frozen=True)
class Enriched:
    """One set's preranked GSEA result, with what the running-score plot draws."""

    name: str
    size: int
    es: float
    nes: float
    p: float
    fdr: float
    positions: tuple[int, ...]
    running: tuple[float, ...]


def _extremes(positions: np.ndarray, weights: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Max and min of the running sum for rows of sorted hit positions, from the hits alone.

    Between hits the running sum only falls, so its maximum sits on a hit and its minimum
    just before one (or at the end, where it returns to 0).

    Args:
        positions: Permutations by hits, sorted ascending, 0-based ranks.
        weights: |score|^weight at those positions, same shape.
        n: Genes in the ranked list.

    Returns:
        (max, min) per row.
    """
    k = positions.shape[1]
    hit_sum = np.cumsum(weights, axis=1) / weights.sum(axis=1, keepdims=True)
    missed = positions - np.arange(k)  # misses strictly before each hit
    at_hit = hit_sum - (missed + 0) / (n - k)
    before_hit = hit_sum - weights / weights.sum(axis=1, keepdims=True) - missed / (n - k)
    maxima = np.maximum(at_hit.max(axis=1), 0.0)
    minima = np.minimum(before_hit.min(axis=1), 0.0)
    return maxima, minima


def _pick(maxima: np.ndarray, minima: np.ndarray) -> np.ndarray:
    """The extreme of larger magnitude; the minimum on a tie, as gseapy does."""
    return np.where(np.abs(maxima) > np.abs(minima), maxima, minima)


def _running(positions: np.ndarray, weights: np.ndarray, n: int) -> tuple[float, ...]:
    """The observed running sum at every rank, for the plot."""
    tags = np.zeros(n)
    tags[positions] = weights / weights.sum()
    misses = np.full(n, 1.0 / (n - len(positions)))
    misses[positions] = 0.0
    return tuple(np.cumsum(tags - misses).tolist())


def collapse(
    genes: Sequence[str], symbols: Sequence[str], scores: Sequence[float], means: Sequence[float]
) -> tuple[list[str], np.ndarray]:
    """One score per symbol: the gene with the highest mean expression, ties to gene id.

    Args:
        genes: Gene identifiers.
        symbols: Their symbols; empty ones are dropped.
        scores: The ranking statistic per gene.
        means: Mean expression per gene, e.g. DESeq2's baseMean.

    Returns:
        Symbols and their scores, in no particular order.
    """
    best: dict[str, tuple[float, str, float]] = {}
    for gene, symbol, score, mean in zip(genes, symbols, scores, means, strict=True):
        if not symbol or math.isnan(score):
            continue
        held = best.get(symbol)
        if held is None or mean > held[0] or (mean == held[0] and gene < held[1]):
            best[symbol] = (mean, gene, score)
    names = sorted(best)
    return names, np.array([best[s][2] for s in names])


def background(genes: Sequence[str], padj: Sequence[float]) -> list[str]:
    """The genes ORA may count: those DESeq2 tested, so with an adjusted p-value."""
    return [g for g, q in zip(genes, padj, strict=True) if q is not None and not math.isnan(q)]


def ranking(
    names: Sequence[str], scores: np.ndarray, weight: float = 1.0
) -> tuple[list[str], np.ndarray]:
    """The list GSEA walks, and the weight each gene carries along it.

    Genes from the most positive score, a tie by name; each weighs |score|^weight. ``prerank``
    and M6's running score on the page (``ctx.function.running``) read this one ordering.

    Args:
        names: Genes, one score each.
        scores: The ranking statistic.
        weight: Exponent on |score|.

    Returns:
        The genes in rank order, and their weights in the same order.
    """
    order = sorted(range(len(names)), key=lambda i: (-scores[i], names[i]))
    return [names[i] for i in order], np.abs(np.asarray(scores, dtype=float)[order]) ** weight


def prerank(  # noqa: PLR0913 - the GSEA parameters, each with its CONVENTIONS default
    names: Sequence[str],
    scores: np.ndarray,
    sets: Mapping[str, Iterable[str]],
    *,
    permutations: int = PERMUTATIONS,
    seed: int = SEED,
    min_size: int = MIN_SIZE,
    max_size: int = MAX_SIZE,
    weight: float = 1.0,
) -> list[Enriched]:
    """Preranked GSEA with gene-set permutation.

    Args:
        names: Genes, one score each; no duplicates (``collapse`` first).
        scores: The ranking statistic, e.g. the Wald statistic.
        sets: Set name to genes.
        permutations: Random sets of the same size per set.
        seed: Seed of the permutation generator.
        min_size: Smallest set, after the cut to the ranked genes.
        max_size: Largest set, after the cut.
        weight: Exponent on |score| (1, the GSEA default).

    Returns:
        One result per set within the size limits, in the order of ``sets``.

    Raises:
        ValueError: If a gene is ranked twice.
    """
    if len(set(names)) != len(names):
        raise ValueError("a gene is ranked twice; collapse the ranking first")
    ranked, magnitude = ranking(names, scores, weight)
    rank_of = {g: i for i, g in enumerate(ranked)}
    kept = sized(sets, set(ranked), min_size, max_size)
    rng = np.random.default_rng(seed)
    n = len(ranked)
    observed, nulls, drawn = [], [], {}
    for name, genes in kept.items():
        positions = np.sort(np.array([rank_of[g] for g in genes]))
        es = _pick(*_extremes(positions[None, :], magnitude[positions][None, :], n))[0]
        random = np.sort(np.array([rng.choice(n, len(genes), replace=False)
                                   for _ in range(permutations)]), axis=1)  # fmt: skip
        observed.append(es)
        nulls.append(_pick(*_extremes(random, magnitude[random], n)))
        drawn[name] = positions
    return _significance(kept, drawn, np.array(observed), np.array(nulls), magnitude, n)


def _significance(  # noqa: PLR0913 - the pieces prerank computed, passed on as they are
    kept: Mapping[str, tuple[str, ...]],
    drawn: Mapping[str, np.ndarray],
    es: np.ndarray,
    null: np.ndarray,
    magnitude: np.ndarray,
    n: int,
) -> list[Enriched]:
    """Nominal p, NES and FDR from the observed scores and their permutation nulls."""
    if not len(es):
        return []
    positive = np.where(null >= 0, null, np.nan)
    negative = np.where(null < 0, null, np.nan)
    # A set with no null of one sign has no NES of that sign: NaN, not a warning.
    with np.errstate(invalid="ignore", divide="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mean_pos, mean_neg = np.nanmean(positive, axis=1), np.nanmean(negative, axis=1)
        nes = np.where(es >= 0, es / mean_pos, -es / mean_neg)
        nes_null = np.where(null >= 0, null / mean_pos[:, None], -null / mean_neg[:, None])
    p = np.where(
        es >= 0,
        (1 + (null >= es[:, None]).sum(axis=1)) / (1 + (null >= 0).sum(axis=1)),
        (1 + (null <= es[:, None]).sum(axis=1)) / (1 + (null < 0).sum(axis=1)),
    )
    pooled, fdr = nes_null.ravel(), []
    for value in nes:
        if value >= 0:
            share_null = (pooled >= value).sum() / max((pooled >= 0).sum(), 1)
            share_obs = (nes >= value).sum() / max((nes >= 0).sum(), 1)
        else:
            share_null = (pooled <= value).sum() / max((pooled < 0).sum(), 1)
            share_obs = (nes <= value).sum() / max((nes < 0).sum(), 1)
        fdr.append(float(min(share_null / share_obs, 1.0)))
    return [
        Enriched(name, len(genes), float(e), float(v), float(q), f,
                 tuple(drawn[name].tolist()), _running(drawn[name], magnitude[drawn[name]], n))
        for (name, genes), e, v, q, f in zip(kept.items(), es, nes, p, fdr, strict=True)
    ]  # fmt: skip
