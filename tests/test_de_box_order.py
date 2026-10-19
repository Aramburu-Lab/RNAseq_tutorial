"""M5's named-gene boxes stand in the dataset file's order on every run (TODO.md, M8 open item)."""

from types import SimpleNamespace

import polars as pl

from ctx.de.page import Words
from ctx.de.run import _box

SAMPLES = ("s1", "s2", "s3", "s4")


def fit() -> SimpleNamespace:
    """The parts of a fitted contrast that ``_box`` reads: two samples per level."""
    return SimpleNamespace(
        contrast=SimpleNamespace(reference="ref", tested="test", paired=False),
        results=pl.DataFrame({"gene_id": ["g1", "g2", "g3", "g4"], "padj": [0.1, 0.2, 0.3, 0.4]}),
        samples=pl.DataFrame(
            {
                "sample": list(SAMPLES),
                "group": ["reference", "reference", "tested", "tested"],
                "donor": ["d1", "d2", "d1", "d2"],
            }
        ),  # fmt: skip
        size_factors=pl.DataFrame({"sample": list(SAMPLES), "size_factor": [1.0] * 4}),
    )


def test_boxes_follow_the_named_order_and_the_first_gene_of_a_symbol():
    """The table holds C, A, B and a second A; the file names A, B, C and a missing Z.

    Mutation caught: the boxes in the gene table's order (C first), the order an unordered
    ``unique`` gave on some runs; and a symbol's second gene taken in place of its first.
    """
    genes = pl.DataFrame({"gene_id": ["g3", "g1", "g2", "g4"], "gene_name": ["C", "A", "B", "A"]})
    counts = pl.DataFrame({"gene_id": ["g1", "g2", "g3", "g4"],
                           **{s: [10.0 * (i + 1), 5.0 + 2 * i, 1.0 + i, 99.0]
                              for i, s in enumerate(SAMPLES)}})  # fmt: skip
    result = fit()
    data, missing = _box(result, counts, genes, ("A", "B", "C", "Z"), Words("t", "test", "ref"))
    assert data is not None and missing == ["Z"]
    assert data.frame["feature"].unique(maintain_order=True).to_list() == ["A", "B", "C"]
    assert [t.feature for t in data.caption.tests] == ["A", "B", "C"]
    # g1's counts reach 40, log2(41) = 5.36; g4's 99 would give log2(100) = 6.64.
    assert max(data.frame.filter(pl.col("feature") == "A")["value"].to_list()) < 6
    assert data.pairs == ""


def test_a_paired_contrasts_boxes_join_each_donor_and_say_so():
    """E.1e: a paired contrast's named-gene boxes draw a line per donor, and the reading guide
    says what the lines join while the test above them ignores the pairing.

    Mutations caught: the pairs never set for a paired contrast; the reading line left out."""
    from ctx.de.page import box_text  # noqa: PLC0415
    from ctx.viz.box import paired_lines  # noqa: PLC0415

    genes = pl.DataFrame({"gene_id": ["g1"], "gene_name": ["A"]})
    counts = pl.DataFrame({"gene_id": ["g1"], **{s: [float(i + 1)] for i, s in enumerate(SAMPLES)}})
    result = fit()
    result.contrast.paired = True
    data, _ = _box(result, counts, genes, ("A",), Words("t", "test", "ref"))
    assert data is not None and data.pairs == "donor"
    assert [u for u, _, _ in paired_lines(data)["A"]] == ["d1", "d2"]
    said = Words("t", "test", "ref")
    joined = box_text(said, 1, 0, paired=True).reading
    assert len(joined) == len(box_text(said, 1, 0).reading) + 1
    assert "line joins each individual's two samples" in joined[-1]
