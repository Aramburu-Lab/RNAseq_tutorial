"""M9's own samples (E.1c2): the box and F1 read M5's size factors and call log2_normalised.

The made-up root carries planted size factors, never median-of-ratios, so a reader that fitted
its own would miss them. Each test names the mutation it catches; each was planted and seen
to fail.
"""

import math

import numpy as np
import polars as pl
import pytest
from genecard_root import (
    ACTB,
    FACTORS,
    NPTX2,
    NPTX2_COUNTS,
    OURS,
    SAMPLES,
    SLUG,
    UPA_ID,
    own_files,
)
from scipy.stats import false_discovery_control, mannwhitneyu, wilcoxon
from test_genecard import _cards, _register, _settings, root  # noqa: F401 - the fixture

from ctx.dataset import load
from ctx.de.run import newest
from ctx.genecard import boxes, samples
from ctx.genecard import run as module
from ctx.genecard.cards import Card
from ctx.report.builder import unpacked
from ctx.viz.save import layout_problems, thin_tick_labels
from ctx.viz.style import THEMES

# Six contrasts in Step-1's words: the longest labels F1 draws side by side.
LONG = ("Neutrophils, G-CSF donors against steady state",
        "Monocytes, G-CSF donors against steady state",
        "Neutrophils, early after HSCT against steady state",
        "Monocytes, early after HSCT against steady state",
        "Neutrophils, G-CSF against steady state, paired by donor",
        "Neutrophils, G-CSF against steady state, the same donors unpaired")  # fmt: skip
CARD = Card(1, "ENSG00000000001.1", "GENE", "named")


def _rows(ref: list[float], tst: list[float], donors: list[str] | None = None) -> pl.DataFrame:
    """One gene's values in one contrast: the reference group, then the tested group."""
    groups = ["reference"] * len(ref) + ["tested"] * len(tst)
    donors = donors or [f"d{i}" for i in range(len(groups))]
    return pl.DataFrame({"donor": donors, "group": groups, "value": [*ref, *tst]})


def _values(contrasts: dict[str, tuple[list[float], list[float]]], paired=(),
            batch: dict[str, list[str]] | None = None) -> samples.Values:  # fmt: skip
    """A Values frame of one gene in each contrast, unpaired unless named in ``paired``."""
    frames = []
    for name, (ref, tst) in contrasts.items():
        n = len(ref) + len(tst)
        donors = [f"d{i}" for i in range(len(ref))] * 2 if name in paired else None
        frame = _rows(ref, tst, donors).with_columns(
            gene_id=pl.lit(CARD.gene_id), contrast=pl.lit(name),
            sample=pl.Series([f"{name}_{i}" for i in range(n)]),
            batch=pl.Series((batch or {}).get(name, [""] * n)),
        )  # fmt: skip
        frames.append(frame.select("gene_id", "contrast", "sample", "donor", "group", "batch",
                                   "value"))  # fmt: skip
    return samples.Values(pl.concat(frames), {c: c in paired for c in contrasts})


def test_values_are_m5s_size_factors_over_the_counts(root):  # noqa: F811
    # Mutations: a size factor fitted here; one contrast's factors used for another's samples;
    # the batch not joined from the sheet.
    folder, spec = root
    own_files(folder, ("a_case", "b_case"))
    got = samples.read(folder, load(spec), [OURS[NPTX2], OURS[ACTB]])
    assert isinstance(got, samples.Values)
    for gene, contrast, sample, donor, group, batch, value in got.frame.rows():
        count = NPTX2_COUNTS[sample] if gene == OURS[NPTX2] else 1000
        want_group, factor = FACTORS[contrast][sample]
        assert value == pytest.approx(math.log2(count / factor + 1), abs=1e-12)
        assert (donor, batch, group) == (*SAMPLES[sample], want_group)
    s1 = got.frame.filter((pl.col("sample") == "s1") & (pl.col("gene_id") == OURS[NPTX2]))
    assert sorted(s1["value"].to_list()) == pytest.approx([math.log2(11), math.log2(21)])
    assert got.frame.height == 2 * (len(FACTORS["a_case"]) + len(FACTORS["b_case"]))


def test_no_size_factor_table_names_why_and_every_card_says_so(root):  # noqa: F811
    # Mutations: a missing table crashing the run; the cards silent about it; a box drawn;
    # the reason naming the table's folder, a path that ran past a phone's page.
    folder, spec = root
    _register(folder, spec)
    why = samples.read(folder, load(spec), [OURS[NPTX2]])
    assert why == "no size_factors table for this dataset; run ctx de again"
    written = module.run(_settings(folder), _cards(), "ctx genecard", SLUG,
                         module.Options(spec, genes=("NPTX2",), offline=True))  # fmt: skip
    page = unpacked(written[0].read_text())
    assert f"M5 normalised values unavailable for NPTX2: {why}." in page
    assert not any("_box_" in p.name for p in written)


def test_unpaired_is_mann_whitney_of_tested_against_reference():
    # Mutations: U taken of the reference; r's sign flipped; r over n1 + n2.
    ref, tst = [1.0, 2.0, 3.0], [2.5, 5.0, 6.0, 7.0]
    got = samples.compare(_rows(ref, tst), paired=False)
    u, p = mannwhitneyu(tst, ref, alternative="two-sided")
    assert got is not None
    assert (got.statistic, got.p) == (pytest.approx(u), pytest.approx(p))
    assert got.r == pytest.approx(2 * u / (3 * 4) - 1) and got.r > 0
    assert (got.n_reference, got.n_tested) == (3, 4)


def test_paired_ranks_each_individuals_difference_matched_by_donor():
    # Mutations: pairs matched by row order; the matched-pairs r replaced by the unpaired one.
    ref, tst = [1.0, 5.0, 2.0, 8.0], [2.0, 4.0, 6.0, 9.5]
    donors = ["a", "b", "c", "d", "d", "c", "b", "a"]  # the tested rows in reverse
    got = samples.compare(_rows(ref, tst, donors), paired=True)
    by = dict(zip(donors[4:], tst, strict=True))
    pairs = [(by[d], r) for d, r in zip(donors[:4], ref, strict=True)]
    want = wilcoxon([t for t, _ in pairs], [r for _, r in pairs])
    assert got is not None
    assert (got.statistic, got.p) == (pytest.approx(want.statistic), pytest.approx(want.pvalue))
    d = np.array([t - r for t, r in pairs])
    assert got.r == pytest.approx(samples.rank_biserial(d))
    assert (got.n_reference, got.n_tested) == (4, 4)


def test_rank_biserial_by_hand():
    # Mutation: zeros kept in the ranks; the negative ranks added.
    assert samples.rank_biserial(np.array([1.0, -2.0, 3.0, 0.0])) == pytest.approx(1 / 3)


def test_too_few_values_or_nothing_to_rank_is_no_test():
    # Mutations: MIN_PER_GROUP lowered to 1; a tied group tested to NaN; zero differences tested.
    assert samples.MIN_PER_GROUP == 2  # noqa: PLR2004 - the rule as registered
    assert samples.compare(_rows([1.0], [2.0, 3.0]), paired=False) is None
    assert samples.compare(_rows([4.0, 4.0], [4.0, 4.0]), paired=False) is None
    assert samples.compare(_rows([1.0, 2.0], [1.0, 2.0], ["a", "b", "a", "b"]), paired=True) is None
    assert samples.compare(_rows([1.0, 2.0], [3.0], ["a", "b", "a"]), paired=True) is None


def test_the_box_family_is_bh_over_the_listed_genes_with_a_test():
    # Mutations: BH over the one gene drawn; a gene no test can rank kept in the family.
    rows = []
    for gene, (ref, tst) in {"g1": ([1.0, 2.0, 3.0], [7.0, 8.0, 9.0]),
                             "g2": ([1.0, 2.0, 3.0], [2.0, 3.0, 4.0]),
                             "g3": ([5.0, 5.0, 5.0], [5.0, 5.0, 5.0])}.items():  # fmt: skip
        rows.append(_rows(ref, tst).with_columns(gene_id=pl.lit(gene), contrast=pl.lit("c"),
                                                 sample=pl.lit("s"), batch=pl.lit("")))  # fmt: skip
    values = samples.Values(pl.concat(rows), {"c": False})
    family = samples.family(values, "c", ["g1", "g2", "g3"])
    assert list(family) == ["g1", "g2"]
    raw = [family[g].p for g in ("g1", "g2")]
    assert [family[g].p_adj for g in ("g1", "g2")] == pytest.approx(
        list(false_discovery_control(raw, method="bh"))
    )
    assert family["g1"].p_adj != pytest.approx(family["g1"].p)


def test_pages_split_contrasts_evenly_at_most_six():
    # Mutations: the cap raised; the last page left with one contrast; a contrast lost.
    assert boxes.PER_FIGURE == 6  # noqa: PLR2004 - six labels fit 170 mm, measured
    sizes = {
        n: [len(p) for p in boxes.pages([str(i) for i in range(n)])] for n in (2, 6, 7, 11, 13)
    }
    assert sizes == {2: [2], 6: [6], 7: [4, 3], 11: [6, 5], 13: [5, 5, 3]}
    assert sum(boxes.pages([str(i) for i in range(11)]), []) == [str(i) for i in range(11)]


def test_f1_adjusts_over_the_views_family_not_the_page():
    # Mutations: BH over the contrasts of one page; a contrast with too few values counted as
    # tested; the n over each box left out.
    contrasts = {f"c{i}": ([1.0, 2.0, 3.0, 4.0], [3.0 + i, 4.0 + i, 5.0 + i, 6.0 + i])
                 for i in range(7)}  # fmt: skip
    contrasts["short"] = ([1.0], [2.0, 3.0])
    values = _values(contrasts)
    labels = {c: f"Contrast {c}" for c in contrasts}
    views = boxes.views(CARD, values, labels, [])
    assert [len(v[1].caption.tests) for v in views] == [4, 3]
    family, left = samples.across(values, CARD.gene_id, list(labels))
    assert left == ["short"]
    for _, data, numbers in views:
        assert numbers["family"] == 7 and numbers["left"] == 1  # noqa: PLR2004
        for t in data.caption.tests:
            name = next(c for c in labels if boxes.label(labels[c]) == t.feature)
            assert t.p_adj == pytest.approx(family[name].p_adj)
            assert data.caption.label(t.feature).startswith("n=4/4\np=")


def test_a_batch_view_is_offered_only_where_it_tests_two_contrasts():
    # Mutations: a batch view with one contrast offered; the batch filter not applied; more
    # batch levels than MAX_BATCHES offered.
    two = {"x": ([1.0, 2.0, 3.0, 4.0], [5.0, 6.0, 7.0, 8.0]),
           "y": ([1.0, 2.0, 3.0, 4.0], [2.0, 3.0, 9.0, 9.5])}  # fmt: skip
    batch = {"x": ["B1", "B1", "B2", "B2", "B1", "B1", "B2", "B2"],
             "y": ["B1", "B1", "B1", "B1", "B1", "B1", "B2", "B2"]}  # fmt: skip
    values = _values(two, batch=batch)
    assert samples.batches(values, boxes.MAX_BATCHES) == ["B1", "B2"]
    views = boxes.views(CARD, values, {"x": "X", "y": "Y"}, ["B1", "B2"])
    assert [n["view"] for _, _, n in views] == ["All samples", "B1 only"]
    b1 = views[1][1].caption.tests
    assert [(t.n_reference, t.n_comparison) for t in b1] == [(2, 2), (4, 2)]
    many = _values(two, batch={c: [f"B{i}" for i in range(8)] for c in two})
    assert samples.batches(many, boxes.MAX_BATCHES) == []


def test_six_long_contrast_labels_draw_whole_in_the_files():
    # Mutations: PER_FIGURE raised to seven or more; LABEL_CHARS widened so a label hides its
    # neighbour's: thin_tick_labels then hides every other label, which layout_problems
    # cannot see.
    contrasts = {f"c{i}": ([1.0, 2.0, 3.0, 4.0], [3.0, 4.0, 5.0 + i, 6.0]) for i in range(6)}
    labels = dict(zip(contrasts, LONG, strict=True))
    (_, data, _), *_ = boxes.views(CARD, _values(contrasts), labels, [])
    for theme in THEMES.values():
        fig = data.static(theme)
        thin_tick_labels(fig)
        shown = [t.get_text() for t in fig.axes[0].get_xticklabels() if t.get_visible()]
        assert shown == [boxes.label(x) for x in LONG]
        assert layout_problems(fig) == []


def test_the_page_plot_names_each_line_and_leaves_room_for_n():
    # Mutations: the category's line breaks left to plotly, which prints them as nothing; the
    # least width measured on the whole name; the n line given no room; the phone's layout kept
    # where it needs more room than upright.
    contrasts = {f"c{i}": ([1.0, 2.0, 3.0, 4.0], [3.0, 4.0, 5.0 + i, 6.0]) for i in range(6)}
    labels = dict(zip(contrasts, LONG, strict=True))
    (_, data, _), *_ = boxes.views(CARD, _values(contrasts), labels, [])
    fig = data.interactive()
    assert list(fig.layout.xaxis.ticktext) == [boxes.label(x).replace("\n", "<br>") for x in LONG]
    plain = boxes.box_interactive(data)
    assert sum(fig.layout.meta["least_width"][0]) < plain.layout.meta["narrow"]["below"]
    assert "narrow" not in fig.layout.meta
    low, high = fig.layout.yaxis.range
    low0, high0 = plain.layout.yaxis.range
    assert low == pytest.approx(low0) and high > high0


def test_the_page_draws_the_box_and_f1_and_ships_every_number(root):  # noqa: F811
    # Mutations: the box's BH over the drawn gene alone; F1 drawn for one contrast; the values
    # or the views' tests not written; F4's page opened alone without its tissue boxes.
    folder, spec = root
    _register(folder, spec)
    own_files(folder, ("a_case", "b_case"))
    results = newest(folder / "outputs" / "05_de", SLUG, "a_case_results")
    for suffix in (".parquet", ".csv"):
        source = results.with_suffix(suffix)
        source.with_name(source.name.replace("a_case", "b_case")).write_bytes(source.read_bytes())
    spec.write_text(spec.read_text() + '\n[[contrasts]]\nname = "b_case"\nvariable = "condition"'
                    '\nreference = "control"\ntested = "case"\nn = [4, 4]\npaired = true\n'
                    'label = "Case against control, paired"\n')  # fmt: skip
    written = module.run(_settings(folder), _cards(), "ctx genecard", SLUG,
                         module.Options(spec, genes=("NPTX2", "ACTB", "UPA"),
                                        offline=True))  # fmt: skip
    names = [p.name for p in written]
    page = unpacked(written[0].read_text())
    assert "NPTX2 in every contrast of this dataset" in page
    assert "No box for NPTX2 in" in page  # M5 tested ACTB and UPA in a_case, not NPTX2
    box = pl.read_csv(next(p for p in written if "_box_upa_" in p.name and p.suffix == ".csv"))
    values = samples.read(folder, load(spec), [OURS[ACTB], UPA_ID])
    assert isinstance(values, samples.Values)
    family = samples.family(values, "a_case", [OURS[ACTB], UPA_ID])
    assert box["p_adj"].unique().to_list() == [pytest.approx(family[UPA_ID].p_adj)]
    assert family[UPA_ID].p_adj == pytest.approx(2 * family[UPA_ID].p)  # BH over the two
    assert not any("_box_nptx2_" in n for n in names)
    tests = pl.read_csv(next(p for p in written if p.name.endswith("_contrast_tests_v1.csv")))
    assert set(tests["test"]) == {"Mann-Whitney U, two-sided",
                                  "Wilcoxon signed-rank over individuals, two-sided"}  # fmt: skip
    own = pl.read_csv(next(p for p in written if p.name.endswith("_own_values_v1.csv")))
    assert own.height == 3 * (len(FACTORS["a_case"]) + len(FACTORS["b_case"]))
    assert any(n.endswith(".pdf") and "_contrasts_nptx2_" in n for n in names)
    alone = next(p for p in written if "_exons_nptx2_" in p.name and p.suffix == ".html")
    assert 'form class="calls subset"' in alone.read_text()
