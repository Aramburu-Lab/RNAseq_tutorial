"""M5's lesson grid, agreement of two contrasts' calls and the boxes' DESeq2 p (WO6), on planted
rows: each rule against a hand count. Each test names the mutations it catches.
"""

import polars as pl

from ctx.de.figures import (
    AGREEMENT,
    BUILD,
    Setting,
    agreement_level,
    counted,
    settings,
)
from ctx.de.stats import ALPHA, LFC_MIN


class _Fit:
    """What the grid reads of a fit: its results."""

    def __init__(self, results: pl.DataFrame) -> None:
        self.results = results


def _results() -> pl.DataFrame:
    # One gene per case: tested and called up; tested at fold 0.6; tested at padj 0.03 under
    # BH but 0.2 under Holm; untested with the smallest raw p of all; tested and down.
    return pl.DataFrame({
        "gene_id": ["a", "b", "c", "d", "e"],
        "log2FoldChange": [2.0, 0.6, 1.2, 3.0, -1.5],
        "pvalue": [1e-6, 1e-5, 0.004, 1e-9, 0.001],
        "padj": [1e-4, 1e-3, 0.03, None, 0.02],
        "padj_by": [1e-3, 1e-2, 0.06, None, 0.04],
        "padj_holm": [1e-4, 1e-3, 0.2, None, 0.02],
        "padj_bonferroni": [1e-4, 1e-3, 0.3, None, 0.02],
    })  # fmt: skip


def test_the_grid_counts_tested_genes_by_the_one_rule():
    """Mutations caught: raw p counted over an untested gene; the fold change taken as strict;
    the cutoff taken as at most; the basis column ignored."""
    fit = _Fit(_results())
    assert counted(fit, BUILD) == (2, 1)
    assert counted(fit, Setting("pvalue", 0.05, 1.0)) == (2, 1)
    assert counted(fit, Setting("padj", 0.05, 0.6)) == (3, 1)
    assert counted(fit, Setting("padj_holm", 0.05, 1.0)) == (1, 1)
    assert counted(fit, Setting("padj", 0.02, 1.0)) == (1, 0)
    assert counted(fit, Setting("padj", 0.05, 2.0)) == (1, 0)


def test_the_build_leads_and_every_setting_is_named_once():
    """Mutations caught: the build not the first option; a setting dropped or repeated."""
    found = settings()
    assert found[0][0] == BUILD == Setting("padj", ALPHA, LFC_MIN)
    assert len({s for s, _ in found}) == len({p for _, p in found}) == 100


def test_agreement_names_each_kind_and_keeps_the_axes_apart():
    """Mutations caught: x-only and y-only swapped; opposite calls counted as the same way; a
    call in one contrast only counted as neither."""
    frame = pl.DataFrame({"x": ["up", "up", "ns", "down", "ns", "down"],
                          "y": ["up", "ns", "down", "up", "ns", "down"]})  # fmt: skip
    found = frame.select(agreement_level(pl.col("x"), pl.col("y")))["literal"].to_list()
    assert found == ["same", "x_only", "y_only", "opposite", "neither", "same"]
    swapped = frame.select(agreement_level(pl.col("y"), pl.col("x")))["literal"].to_list()
    assert swapped == ["same", "y_only", "x_only", "opposite", "neither", "same"]
    assert set(found) == set(AGREEMENT)
