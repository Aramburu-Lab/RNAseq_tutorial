"""Statistical annotation shared by the static and the interactive rendering of a figure.

One ``Caption`` object produces every string either renderer prints, so the PDF, the HTML
and the README cannot state different numbers.
"""

import math
import textwrap
from dataclasses import dataclass

import polars as pl
from scipy.stats import false_discovery_control, mannwhitneyu

TEST_NAME = "Mann-Whitney U, two-sided"


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
        return f"p={_fmt(t.p_adj)}\nr={t.r:.2f}"

    def text(self, width: int = 0) -> str:
        """The caption block: test, correction, n per group, and the sign of r.

        Args:
            width: Characters per line; 0 leaves the caption on one line.

        Returns:
            The caption.
        """
        reference, comparison = self.groups
        ns = sorted({(t.n_reference, t.n_comparison) for t in self.tests})
        n_text = "; ".join(f"n = {a} {reference}, {b} {comparison}" for a, b in ns)
        body = (
            f"{TEST_NAME}; p is BH-adjusted over {len(self.tests)} plotted features; "
            f"{n_text}; rank-biserial r > 0 when {comparison} is higher."
        )
        # Never split "BH-adjusted": a reader, or a PDF text extractor, would rejoin it wrongly.
        return textwrap.fill(body, width, break_on_hyphens=False) if width else body


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
