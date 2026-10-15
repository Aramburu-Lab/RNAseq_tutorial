"""M9's F2 (figure controls plan): the tissue panel's donor views, each from its own summaries.

D-056 keeps GTEx as summaries per tissue and sex, so the donor control is a select between
views Python drew (the D-051 view kind; the page's twin is a lookup). Each test names the
mutation it catches; each was planted and seen to fail.
"""

import polars as pl

from ctx.genecard import figures
from ctx.genecard.cards import Card
from ctx.report.builder import FigureText
from ctx.report.views import viewed
from ctx.viz.save import figure_block

# Tissue, sex, n, median: B has female donors only (a uterus), C male only (a prostate). Male
# donors rank B absent and A above D; female donors rank D above A.
ROWS = [("A", "all", 30, 10.0), ("A", "female", 10, 4.0), ("A", "male", 20, 12.0),
        ("B", "all", 8, 50.0), ("B", "female", 8, 50.0),
        ("C", "all", 6, 2.0), ("C", "male", 6, 2.0),
        ("D", "all", 12, 6.0), ("D", "female", 6, 9.0), ("D", "male", 6, 3.0)]  # fmt: skip


def _card() -> Card:
    frame = pl.DataFrame(
        [(t, s, n, m, m * 0.5, m * 1.5) for t, s, n, m in ROWS],
        schema=["tissue", "sex", "n", "median", "q1", "q3"],
        orient="row",
    )
    return Card(1, "ENSG00000000001.1", "GENE", "named", tissues=frame)


def test_each_view_draws_its_own_donors_in_the_all_donor_order():
    # Mutations: a sex view drawing the other sex's rows; the all-donor view drawing the
    # sexes; a view reordered by its own medians, so a bar's place changes between views.
    views = {label: (data, numbers) for label, data, numbers in figures.tissue_views(_card())}
    assert list(views) == ["By sex", "All donors", "Female donors", "Male donors"]
    for label, sexes in (("All donors", {"all"}), ("Female donors", {"female"}),
                         ("Male donors", {"male"}), ("By sex", {"female", "male"})):  # fmt: skip
        frame = views[label][0].frame
        assert set(frame["sex"]) == sexes
        tissues = list(dict.fromkeys(frame["tissue"]))
        order = [t for t in ("B", "A", "D", "C") if t in tissues]
        assert tissues == order, label


def test_each_view_states_its_own_numbers():
    # Mutations: a sex view's lead read from the all-donor medians; its samples counted over
    # both sexes; a sex with no sample in a tissue given a zero bar.
    views = {label: numbers for label, _, numbers in figures.tissue_views(_card())}
    female, male, built = views["Female donors"], views["Male donors"], views["By sex"]
    assert (female["top"], female["bottom"], female["tissues"], female["samples"]) == (
        "B",
        "A",
        3,
        "24",
    )
    assert (male["top"], male["bottom"], male["tissues"], male["samples"]) == ("A", "C", 3, "32")
    assert (built["top"], built["bottom"], built["samples"]) == ("B", "C", "56")
    assert "female donors" in figures.tissues(_card(), ("female",))[0].caption


def test_the_views_become_one_select_with_the_built_figure_first():
    # Mutation: a view list whose first entry is not the figure as built, which viewed refuses.
    views = [(label, data, FigureText(title=label)) for label, data, _ in
             figures.tissue_views(_card())]  # fmt: skip
    block = viewed(figure_block(views[0][1], "x", views[0][2]), views, "Donors")
    controls = block.controls or {}
    assert controls["kind"] == "view"
    assert [o["label"] for o in controls["options"]] == [v[0] for v in views]
    assert controls["options"][2]["caption"].count("female donors") == 1


def test_a_card_with_one_sex_offers_no_view_of_the_other():
    # Mutation: a male view drawn for a gene whose GTEx rows hold no male donor.
    card = _card()
    card.tissues = card.tissues.filter(pl.col("sex") != "male")  # type: ignore[union-attr]
    labels = [label for label, _, _ in figures.tissue_views(card)]
    assert labels == ["By sex", "All donors", "Female donors"]
