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
from ctx.process import blocks, figures, views
from ctx.process import tables as tab
from ctx.viz.contour import DENSITY_MIN_POINTS
from ctx.viz.style import THEMES


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


def _sums(n_coding: int = 400, n_other: int = 12) -> pl.DataFrame:
    """Gene sums of two types: many protein-coding genes, a few others, b near a."""
    rng = np.random.default_rng(3)
    a = rng.integers(0, 5000, n_coding + n_other).astype(float)
    b = np.maximum(a * rng.lognormal(0, 0.3, a.size), 0).round()
    kinds = ["protein_coding"] * n_coding + ["other"] * n_other
    return pl.DataFrame({"gene_id": [f"g{i}" for i in range(a.size)], "a": a, "b": b,
                         "gene_name": [f"G{i}" for i in range(a.size)],
                         "type_group": kinds})  # fmt: skip


def test_the_page_draws_densities_and_the_genes_beyond_the_limits():
    """F3 on the page: each large type's regions and the genes outside as dots; files all genes.

    Mutations caught: the page drawing every gene again (outside ignored); a type under
    ``DENSITY_MIN_POINTS`` drawn as a density; the figure files losing the genes inside; the
    two renderers drawing different isolines; the caption not counting what the page draws.
    """
    sums = _sums()
    d = (pl.col("b") + 1).log(2) - (pl.col("a") + 1).log(2)
    kept = sums.filter(pl.col("a") + pl.col("b") > 0).select(d=d)["d"]
    limits = (float(kept.mean()) - 0.3, float(kept.mean()) + 0.3)
    data = figures.counter_scatter(sums, figures.comparison_words("salmon"), "A caption.",
                                   limits)  # fmt: skip
    frame = data.frame
    outside = frame.filter(pl.col("outside"))
    assert 0 < outside.height < frame.height
    assert outside.height == int(((kept < limits[0]) | (kept > limits[1])).sum())
    assert list(data.densities) == ["protein_coding"]
    assert 0 < frame.filter(pl.col("type_group") == "other").height < DENSITY_MIN_POINTS
    drawn = frame.filter(pl.col("outside") | (pl.col("type_group") == "other"))
    assert data.on_page()["gene_id"].to_list() == drawn["gene_id"].to_list()
    page = data.interactive()
    points = [t for t in page.data if t.mode == "markers"]
    assert sum(len(t.x) for t in points) == drawn.height
    lines = [t for t in page.data if t.mode == "lines" and t.legendgroup in ("inner", "outer")]
    assert {t.legendgroup for t in lines} == {"inner", "outer"}
    ax = data.static(THEMES["white"]).axes[0]
    assert sum(len(c.get_offsets()) for c in ax.collections) == frame.height
    drawn_lines = [ln.get_xydata() for ln in ax.lines][-len(data.regions) :]
    for region, line in zip(data.regions, drawn_lines, strict=True):
        np.testing.assert_array_equal(line[:, 0], region.x)
    on_page = [x for t in lines for x in t.x if x is not None]
    assert on_page == [x for r in data.regions for x in np.round(r.x, 4).tolist()]
    caption = data.caption_text()
    small = drawn.height - outside.height
    assert f"the {outside.height:,} genes beyond" in caption
    assert f"and {small:,} of smaller types" in caption and f"all {frame.height:,}" in caption
    layers = data.export_layers()
    assert layers["thinned"] == {} and layers["caption"] == caption


def test_without_limits_the_page_thins_as_before():
    """No limits (too few genes): no densities, every gene that differs kept, no region line."""
    data = figures.counter_scatter(_sums(), figures.comparison_words("salmon"), "A caption.")
    assert data.densities == {} and data.regions == []
    assert data.frame["outside"].null_count() == data.frame.height
    assert data.caption_text() == "A caption."
