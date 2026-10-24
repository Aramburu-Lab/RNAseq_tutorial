"""M3's scatter views (WO8) offer a sample axis only up to ``SAMPLES_MAX`` samples, each view
carrying its genes again. The fixture page has fewer, so the cap is held here on a made-up
comparison. Mutation caught: the sample axis offered past its cap.
"""

from types import SimpleNamespace

import numpy as np
import polars as pl
from polars.testing import assert_frame_equal

from ctx import concordance
from ctx.figures.recorder import Recorder
from ctx.process import blocks, views
from ctx.process import tables as tab


def _part(samples: int) -> pl.DataFrame:
    rows = [(f"s{i}", f"g{j}", j + 1, j + 2) for i in range(samples) for j in range(4)]
    frame = pl.DataFrame(rows, schema=["sample", "gene_id", "a", "b"], orient="row")
    return frame.with_columns(type_group=pl.lit("protein_coding"))


def test_the_sample_axis_stops_at_its_cap():
    for n, units in ((views.SAMPLES_MAX, views.SAMPLES_MAX + 1), (views.SAMPLES_MAX + 1, 1)):
        found, axes, at = views.subsets(_part(n))
        assert len(axes[1][1]) == units
        assert {s.sample for s in found} >= {concordance.ALL}
        assert len(found) == len(at) == len(axes[0][1]) * units


def _typed_comparison(samples: int = 3, per_type: int = 12) -> tuple[pl.DataFrame, pl.DataFrame]:
    """A made-up comparison of three gene types, long form, and the genes' names and types."""
    rng = np.random.default_rng(7)
    types = [*tab.GENE_TYPES, "snRNA"]
    genes = pl.DataFrame({"gene_id": [f"g{t}{j}" for t in range(3) for j in range(per_type)],
                          "gene_name": [f"G{t}{j}" for t in range(3) for j in range(per_type)],
                          "gene_type": [types[t] for t in range(3)
                                        for _ in range(per_type)]})  # fmt: skip
    rows = [("salmon", f"s{i}", g, float(a), float(max(a + d, 0)))
            for i in range(samples) for g in genes["gene_id"]
            for a, d in [(rng.integers(1, 500), rng.integers(-20, 20))]]  # fmt: skip
    long = pl.DataFrame(rows, schema=["comparison", "sample", "gene_id", "a", "b"], orient="row")
    return long, genes


def test_each_gene_type_is_measured_once_and_its_rows_are_per_unit_own(monkeypatch):
    """Every view of F3 and F4 reads one ``per_unit`` per gene type, equal to its own call.

    Mutations caught: ``Agreed.of`` computing again on every call, so each type is measured once
    per view (the 2,579 s of the Breast-tissue-BRCA profile); and the every-gene rows computed
    again instead of taken from the concordance table.
    """
    long, genes = _typed_comparison()
    agreements = tab.agreements(long, genes)
    calls = []
    real = concordance.per_unit

    def counted(frame, *args, **kwargs):
        calls.append(frozenset(frame["gene_id"]))
        return real(frame, *args, **kwargs)

    monkeypatch.setattr(concordance, "per_unit", counted)
    found = SimpleNamespace(agreements=agreements, outcomes={f"s{i}": {} for i in range(3)},
                            comparison=long)  # fmt: skip
    blocks.agreement_figures(Recorder("03_process", "made-up"), found, SimpleNamespace(prefix=""))
    typed = views.typed(long, agreements.sums)
    by_type = [frozenset(views.genes_of(typed, t)["gene_id"]) for t in views.TYPES[1:]]
    assert sorted(calls, key=sorted) == sorted(by_type, key=sorted)
    agreed = views.Agreed(typed, {None: views.from_table(agreements.concordance, "salmon")})
    for gene_type in views.TYPES:
        before = real(views.genes_of(typed, gene_type).select("gene_id", "sample", "a", "b"),
                      "a", "b", unit="sample", gene="gene_id")  # fmt: skip
        assert_frame_equal(agreed.of(gene_type), before)
