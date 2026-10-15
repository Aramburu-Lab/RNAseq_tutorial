"""M9's exon panel (H.3): the reference over expressing tissues, steps without a pseudocount.

On the 2026-10-14 lung cards a median over all 68 GTEx tissues gave NKX2-1's E1, E2 and E4 a
reference of 0, so the gene's lung row drew empty, and exons with no reads printed "nan" as
their CV. Each test names the mutation it catches; each was planted and seen to fail.
"""

import math

import polars as pl

from ctx.genecard import annotation, figures, gtex
from ctx.genecard.cards import Card

# Tissues A, B and D express the gene (median TPM at or above 1); C does not.
TPM = {"A": 40.0, "B": 2.0, "C": 0.2, "D": 1.0}
# Exon shares by tissue. E1 has a reference of 0.02 over A, B, D, but 0.03 with C's 0.9 in;
# E2 has no reads where the gene is expressed; E3 is read in A only, so its reference is 0.
SHARES = {"E1": {"A": 0.04, "B": 0.0, "C": 0.9, "D": 0.02},
          "E2": {"A": 0.0, "B": 0.0, "C": 0.1, "D": 0.0},
          "E3": {"A": 0.05, "B": 0.0, "C": 0.0, "D": 0.0}}  # fmt: skip


def _spread(tpm: dict[str, float]) -> pl.DataFrame:
    rows = [(t, "all", 10, m, m, m) for t, m in tpm.items()]
    return pl.DataFrame(rows, schema=["tissue", "sex", "n", "median", "q1", "q3"], orient="row")


def _shares() -> pl.DataFrame:
    rows = [(e, t, 10, v) for e, by in SHARES.items() for t, v in by.items()]
    return pl.DataFrame(rows, schema=["exon", "tissue", "n", "share"], orient="row")


def _card(tpm: dict[str, float]) -> Card:
    return Card(1, "ENSG00000000001.1", "GENE", "named", tissues=_spread(tpm), exons=_shares())


def test_detected_tissues_reach_one_tpm():
    # Mutation: the floor taken as strictly above 1, which drops D at exactly 1 TPM.
    assert gtex.detected(_spread(TPM)) == ["A", "B", "D"]


def test_reference_is_taken_after_the_tissues_are_restricted():
    # Mutation: the median over every tissue, then the restriction (E1's reference 0.03, not
    # 0.02); the order is the rule, so this asserts it, not only the floor.
    use = annotation.exon_use(_shares(), gtex.detected(_spread(TPM)))
    reference = dict(use.frame.select("exon", "reference").unique().iter_rows())
    assert reference == {"E1": 0.02, "E2": 0.0, "E3": 0.0}
    assert use.tissues == ("A", "B", "D") and use.left_out == ("C",)
    assert "C" not in use.frame["tissue"].to_list()


def test_steps_take_the_outer_levels_without_a_pseudocount():
    # Mutations: a share over a reference of 0 left blank (NKX2-1's lung row) or given a
    # pseudocount; no share against a reference above 0 left blank; 0 over 0 given a step.
    use = annotation.exon_use(_shares(), ["A", "B", "D"])
    cells = {(e, t): s for (e, t), s in zip(use.frame.select("exon", "tissue").iter_rows(),
                                             use.steps, strict=True)}  # fmt: skip
    assert cells[("E1", "A")] == 1.0 and cells[("E1", "D")] == 0.0
    assert cells[("E1", "B")] == -math.inf
    assert cells[("E3", "A")] == math.inf
    assert math.isnan(cells[("E3", "B")]) and math.isnan(cells[("E2", "A")])
    ratios = dict(zip(use.frame.select("exon", "tissue").iter_rows(), use.frame["log2_ratio"],
                      strict=True))  # fmt: skip
    assert ratios[("E3", "A")] is None and ratios[("E1", "B")] is None  # never an infinity


def test_undefined_cv_is_none_and_printed_as_nothing():
    # Mutations: a CV of NaN for an exon without reads, printed "nan" beside its row; a CV
    # over one tissue; the qPCR pair taking an exon with no CV.
    use = annotation.exon_use(_shares(), ["A", "B", "D"])
    assert use.cv["E2"] is None and use.cv["E1"] is not None
    assert annotation.exon_cv(_shares().filter(pl.col("tissue") == "A")) == dict.fromkeys(SHARES)
    assert annotation.qpcr_pair({"E1": 0.3, "E2": None, "E3": 0.2}) is None
    drawn = figures.exons(_card(TPM))
    assert drawn is not None
    heat, numbers, _ = drawn
    assert "E2" not in (heat.row_stats or {})
    stats = heat.interactive().layout.yaxis2.ticktext
    assert "nan" not in " ".join(stats) and stats[1] == ""
    table = heat.table()
    assert table.filter(pl.col("exon") == "E2")["row_stat"].null_count() == 3  # noqa: PLR2004
    levels = dict(zip(table.select("exon", "tissue").iter_rows(), table["level"], strict=True))
    assert levels[("E3", "A")] == "≥ 2" and levels[("E1", "B")] == "< -2"
    assert levels[("E2", "A")] is None
    assert (numbers["tissues"], numbers["left"], numbers["unreferenced"]) == (3, 1, 2)


def test_fewer_than_two_expressing_tissues_draw_no_panel():
    # Mutation: a one-tissue panel whose every ratio is 0 by construction.
    assert figures.exons(_card({"A": 5.0, "B": 0.5, "C": 0.2, "D": 0.1})) is None
