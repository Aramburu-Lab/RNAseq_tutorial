"""Statistical annotation shared by the static and the interactive rendering of a figure.

One ``Caption`` object produces every string either renderer prints, so the PDF, the HTML
and the README cannot state different numbers. The words are ``content/annotate.toml``'s.
"""

import math
import textwrap
from dataclasses import dataclass
from typing import Any

import numpy as np
import polars as pl
from scipy.stats import false_discovery_control, mannwhitneyu

from ctx.report import text as content
from ctx.report.marks import plain
from ctx.viz.metrics import words, wrap_to
from ctx.viz.style import (
    AXIS_TITLE_SHARE,
    CAPTION_BOTTOM,
    CAPTION_INSET,
    LEGEND_MAX_COLUMNS,
    LINE_SPACING,
    SIZES,
)

PAGE = "annotate"
TEST_NAME = content.word(PAGE, "test_name")
PERMANOVA_NAME = content.word(PAGE, "permanova_name")


@dataclass(frozen=True)
class FeatureTest:
    """One feature's comparison of the two groups."""

    feature: str
    n_reference: int
    n_comparison: int
    u: float
    p: float
    p_adj: float
    r: float


@dataclass(frozen=True)
class Caption:
    """Test results for every plotted feature, and the text both renderers print."""

    groups: tuple[str, str]
    tests: tuple[FeatureTest, ...]

    def label(self, feature: str) -> str:
        """The annotation printed above one feature."""
        t = next(t for t in self.tests if t.feature == feature)
        return content.fill(PAGE, "label", p=_fmt(t.p_adj), r=t.r)

    def text(self, width: int = 0) -> str:
        """The caption block: test, correction, n per group, and the sign of r.

        Args:
            width: Characters per line; 0 leaves the caption on one line.

        Returns:
            The caption.
        """
        reference, comparison = self.groups
        ns = sorted({(t.n_reference, t.n_comparison) for t in self.tests})
        n_text = "; ".join(
            content.fill(
                PAGE,
                "sizes",
                n_reference=a,
                reference=reference,
                n_comparison=b,
                comparison=comparison,
            )
            for a, b in ns
        )
        body = content.fill(PAGE, "caption", test=TEST_NAME, features=len(self.tests),
                            sizes=n_text, comparison=comparison)  # fmt: skip
        return wrap(body, width)


def wrap(text: str, width: int = 0) -> str:
    """Wrap a caption to ``width`` characters; 0 leaves it on one line.

    Never splits at a hyphen: "BH-adjusted" broken across lines would be rejoined wrongly by a
    reader or a PDF text extractor.
    """
    return textwrap.fill(text, width, break_on_hyphens=False) if width else text


def fit_title(text: str, height: float) -> str:
    """Wrap a vertical axis title to the length its axis has: a share of the figure height.

    Measured as the files set it, so a line ends only where its next word would cross that
    length: an average glyph width broke M3's "Difference d (log2(count + 1))", 130.2 pt
    against 145.8 pt, after "d". A parenthesised span, such as that unit, is never broken:
    split at its spaces it left ``1))`` alone on a line. A span longer than a line takes its
    own. The page's titles break at the same places (``ctxExport.words``).

    Args:
        text: The title, with its unit.
        height: Figure height in inches.

    Returns:
        The title, broken into lines that fit.
    """
    room = height * 72 * AXIS_TITLE_SHARE
    return "\n".join(wrap_to(words(text), room, SIZES["axis_label"]))


def place_caption(
    fig: Any,
    text: str,
    width: float,
    height: float,
    *,
    ax: Any = None,
    handles: list[Any] | None = None,
) -> None:
    """Print a wrapped caption under a figure, then lay out the axes and legend above it.

    With ``handles``, the legend goes above the axes in the most columns that keep it within
    the axes' width, measured after layout: a wider legend makes the layout squeeze the axes.

    Args:
        fig: The matplotlib Figure.
        text: The caption, unwrapped.
        width: Figure width in inches; sets the wrap.
        height: Figure height in inches.
        ax: The axes that own the legend; required with ``handles``.
        handles: Legend handles, or None when the caller draws its own legend.
    """
    from matplotlib.backends.backend_agg import FigureCanvasAgg  # noqa: PLC0415 - heavy

    # Wrapped as printed: a glossary mark is its display text off the page (D-048).
    # Measured as the files set it, so a line ends only where its next word would cross the
    # inset: an average glyph width left about a word of space on every line (Rytis,
    # 2026-10-01). A page's own PDF of a changed view wraps by the same measure (D-033).
    usable = width * 72 * (1 - 2 * CAPTION_INSET)
    caption = "\n".join(wrap_to(plain(text), usable, SIZES["caption"]))
    fig.text(CAPTION_INSET, CAPTION_BOTTOM, caption, ha="left", va="bottom",
             fontsize=SIZES["caption"])  # fmt: skip
    caption_inches = (caption.count("\n") + 1) * SIZES["caption"] * LINE_SPACING / 72
    rect = (0, caption_inches / height, 1, 1)
    if not handles:
        fig.tight_layout(rect=rect)
        return
    renderer = FigureCanvasAgg(fig).get_renderer()  # type: ignore[no-untyped-call]
    for ncols in range(min(len(handles), LEGEND_MAX_COLUMNS), 0, -1):
        legend = ax.legend(
            handles=handles, loc="lower left", bbox_to_anchor=(0, 1), ncols=ncols,
            borderaxespad=0,
        )  # fmt: skip
        fig.tight_layout(rect=rect)
        if legend.get_window_extent(renderer).x1 <= ax.get_window_extent(renderer).x1 + 1:
            break


def _fmt(p: float) -> str:
    """An exact p-value to two significant digits, never a threshold or stars."""
    return f"{p:.2g}"


def compare_groups(
    frame: pl.DataFrame, *, feature: str, group: str, value: str, groups: tuple[str, str]
) -> Caption:
    """Compare two groups per feature: drop missing, count, test, effect size, then BH.

    BH runs over the features in ``frame`` only, the ones the figure plots (TODO.md,
    Figure and statistics conventions).

    Args:
        frame: Long table with one row per sample and feature.
        feature: Column naming the feature (e.g. gene).
        group: Column holding the group label.
        value: Column with the measured value.
        groups: (reference, comparison); r > 0 means the comparison group tends higher.

    Returns:
        The caption with one test per feature, in order of first appearance.

    Raises:
        ValueError: If a feature has fewer than one non-missing value in either group.
    """
    reference, comparison = groups
    clean = frame.filter(pl.col(value).is_not_null() & pl.col(value).is_not_nan())
    raw = []
    for name in frame[feature].unique(maintain_order=True).to_list():
        rows = clean.filter(pl.col(feature) == name)
        x = rows.filter(pl.col(group) == comparison)[value].to_numpy()
        y = rows.filter(pl.col(group) == reference)[value].to_numpy()
        if len(x) == 0 or len(y) == 0:
            raise ValueError(f"{name}: a group has no values after dropping missing ones")
        u, p = mannwhitneyu(x, y, alternative="two-sided")
        raw.append((name, len(y), len(x), float(u), float(p), 2 * float(u) / (len(x) * len(y)) - 1))
    adjusted = false_discovery_control([p for *_, p, _ in raw], method="bh")
    tests = tuple(
        FeatureTest(name, n_ref, n_cmp, u, p, float(p_adj), r)
        for (name, n_ref, n_cmp, u, p, r), p_adj in zip(raw, adjusted, strict=True)
    )
    if any(math.isnan(t.p) for t in tests):
        raise ValueError("a test returned NaN; check for constant groups")
    return Caption(groups, tests)


@dataclass(frozen=True)
class Permanova:
    """One PERMANOVA: pseudo-F, the share of variance the labels explain, and its p-value."""

    f: float
    r2: float
    p: float
    n: int
    groups: int
    permutations: int
    stratified: bool

    def text(self, factor: str) -> str:
        """The test as a caption sentence for one factor."""
        within = content.word(PAGE, "permanova.within") if self.stratified else ""
        values = {"name": PERMANOVA_NAME, "factor": factor, "r2": self.r2, "f": self.f,
                  "n": self.n, "groups": self.groups}  # fmt: skip
        if math.isnan(self.p):
            return content.fill(PAGE, "permanova.untested", **values)
        return content.fill(PAGE, "permanova.tested", p=_fmt(self.p), within=within,
                            permutations=self.permutations, **values)  # fmt: skip


def _pseudo_f(coords: np.ndarray, codes: np.ndarray, n_groups: int) -> tuple[float, float]:
    """Pseudo-F and R squared of one labelling, from within-group sums of squares."""
    total = float(((coords - coords.mean(axis=0)) ** 2).sum())
    counts = np.bincount(codes, minlength=n_groups)
    sums = np.zeros((n_groups, coords.shape[1]))
    np.add.at(sums, codes, coords)
    # Within = sum of squares minus each group's squared sum over its size.
    within = float((coords**2).sum() - ((sums**2).sum(axis=1) / counts).sum())
    between = total - within
    n = coords.shape[0]
    return (between / (n_groups - 1)) / (within / (n - n_groups)), between / total


def permanova(
    coords: np.ndarray,
    labels: list[str],
    strata: list[str] | None = None,
    permutations: int = 999,
    seed: int = 42,
) -> Permanova:
    """PERMANOVA on Euclidean distances, labels permuted freely or within strata.

    With Euclidean distances the pseudo-F is the trace ANOVA F of the coordinates. Samples
    that share an individual are not exchangeable across individuals, so pass the individual
    as ``strata``: labels then move only between samples of one individual.

    Args:
        coords: Samples by dimensions, e.g. principal component scores.
        labels: One group label per sample.
        strata: One stratum per sample, or None to permute freely.
        permutations: Label permutations for the null distribution.
        seed: Seed of the permutation generator.

    Returns:
        The test; p = (1 + permutations with F at least the observed) / (1 + permutations),
        or NaN when strata leave no label free to move: every permutation would be the
        observed labelling, and p = 1 would read as a result rather than as no test.

    Raises:
        ValueError: With fewer than two groups, or as many groups as samples.
    """
    levels, codes = np.unique(np.asarray(labels), return_inverse=True)
    n = len(codes)
    if not 1 < len(levels) < n:
        raise ValueError(f"{len(levels)} groups in {n} samples; PERMANOVA needs 2 to n-1")
    observed, r2 = _pseudo_f(coords, codes, len(levels))
    rng = np.random.default_rng(seed)
    blocks = (
        [np.flatnonzero(np.asarray(strata) == s) for s in np.unique(np.asarray(strata))]
        if strata is not None
        else [np.arange(n)]
    )
    if not any(len(set(codes[block])) > 1 for block in blocks):
        return Permanova(observed, r2, math.nan, n, len(levels), permutations, True)
    at_least = 0
    for _ in range(permutations):
        shuffled = codes.copy()
        for block in blocks:
            shuffled[block] = codes[rng.permutation(block)]
        at_least += _pseudo_f(coords, shuffled, len(levels))[0] >= observed - 1e-12 * abs(observed)
    p = (1 + at_least) / (1 + permutations)
    return Permanova(observed, r2, p, n, len(levels), permutations, strata is not None)
