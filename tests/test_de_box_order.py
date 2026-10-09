"""M5's named-gene boxes stand in the dataset file's order on every run (TODO.md, M8 open item)."""

from types import SimpleNamespace

import polars as pl

from ctx.de.run import _box

SAMPLES = ("s1", "s2", "s3", "s4")


def fit() -> SimpleNamespace:
    """The parts of a fitted contrast that ``_box`` reads: two samples per level."""
    return SimpleNamespace(
        contrast=SimpleNamespace(reference="ref", tested="test"),
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
    data, missing = _box(fit(), counts, genes, ("A", "B", "C", "Z"))
    assert data is not None and missing == ["Z"]
    assert data.frame["feature"].unique(maintain_order=True).to_list() == ["A", "B", "C"]
    assert [t.feature for t in data.caption.tests] == ["A", "B", "C"]
    # g1's counts reach 40, log2(41) = 5.36; g4's 99 would give log2(100) = 6.64.
    assert max(data.frame.filter(pl.col("feature") == "A")["value"].to_list()) < 6
