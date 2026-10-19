"""M5's lesson grid, agreement colours and the boxes' DESeq2 p (WO6, F8 to F12), in Chromium.

Every option of every view on the pipeline fixture's M5 page is set as a reader sets it and read
back from plotly, at 1,280 px light and 375 px dark, and checked against numbers computed here
from the results tables with the one call rule (``ctx.de.stats.classify``), never against the
page's own JSON alone. Each test names the mutations it catches.
"""

import re

import polars as pl
import pytest
from test_export_browser import CHROMIUM
from views_browser import exported, walk

from ctx.de.figures import ALPHAS, BASES, FOLDS
from ctx.de.stats import classify

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
# The rule a lesson's caption states: the basis, the cutoff and the fold change.
RULE = r"\((\S+)(?:-adjusted)? p below ([\d.]+), log2 FC at or beyond ±([\d.]+)\)"


@pytest.fixture(scope="module")
def m5(pipeline, tmp_path_factory):
    """(the M5 folder, its overview page, every view walked at both sizes)."""
    folder = pipeline[0] / "outputs" / "05_de"
    (page,) = folder.glob("*_de_desynth_overview_v1.html")
    return folder, page, walk(page, tmp_path_factory.mktemp("walk"))


def _table(folder, what: str) -> pl.DataFrame:
    (path,) = folder.glob(f"*_de_desynth_{what}_v1.parquet")
    return pl.read_parquet(path)


def _tested(folder, contrast: str) -> pl.DataFrame:
    return _table(folder, f"{contrast}_results").filter(pl.col("padj").is_not_null())


def _views(seen: dict, label: str) -> list[dict]:
    return [v for v in seen["acted"]["views"] if v["label"] == label]


def _counts(option: dict) -> dict[str, dict[str, int]]:
    """Up and down per contrast label, as plotly drew a horizontal stacked bar."""
    return {t["name"]: dict(zip(t["y"], t["x"], strict=True)) for t in option["traces"]}


@pytest.mark.parametrize("size", ["wide_light", "phone_dark"])
def test_every_grid_setting_draws_what_classify_counts(m5, size):
    """Mutations caught: an option drawn under another's key; the build not first; raw p read
    over untested genes; the cutoff applied as at most rather than below; a lead naming the
    second contrast."""
    folder, _, seen = m5
    (grid,) = _views(seen[size], "Calls by")
    labels = dict(pl.read_csv(next(folder.glob("*_de_desynth_calls_batch_v1.csv")))
                  .select("label", "contrast").unique(maintain_order=True).rows())  # fmt: skip
    keys = [(b, a, f) for b in BASES for a in ALPHAS for f in FOLDS]
    # Rytis's plan: 5 bases, 4 cutoffs, 5 fold changes.
    assert len(grid["options"]) == len(keys) == 100
    build = grid["options"][0]
    assert build["line"] == "" and not build["moved"]
    for option in grid["options"][1:]:
        assert option["moved"] and option["line"].startswith("This view shows ")
        found = re.search(r"shows (\S+(?: p)?) below ([\d.]+), log2 FC at least ([\d.]+);",
                          option["line"])  # fmt: skip
        assert found, option["line"]
    for option in grid["options"]:
        rule = re.search(RULE, option["caption"])
        assert rule, option["caption"]
        word, alpha, fold = rule.group(1), float(rule.group(2)), float(rule.group(3))
        basis = {"BH": "padj", "BY": "padj_by", "Holm": "padj_holm",
                 "Bonferroni": "padj_bonferroni", "raw": "pvalue"}[word.split("-")[0]]  # fmt: skip
        drawn = _counts(option)
        totals = {}
        for label, contrast in labels.items():
            calls = (
                _tested(folder, contrast)
                .select(classify(pl.col(basis), pl.col("log2FoldChange"), alpha, fold))
                .to_series()
            )
            up, down = int((calls == "up").sum()), int((calls == "down").sum())
            assert drawn["up"][label] == up and drawn["down"][label] == down, (basis, alpha, fold)
            totals[label] = up + down
        most = max(totals.values())
        if most:
            first = next(lab for lab in labels if totals[lab] == most)
            assert option["lead"].startswith(f"{first} calls the most genes: {most:,}")
        else:
            assert option["lead"] == "No contrast calls a gene at this setting."
    assert grid["back"]["traces"] == build["traces"] and not grid["back"]["moved"]


def _agreement(folder, x: str, y: str) -> dict[str, int]:
    a = _tested(folder, x).select("gene_id", pl.col("call").alias("cx"),
                                  pl.col("log2FoldChange").alias("lx"))  # fmt: skip
    b = _tested(folder, y).select("gene_id", pl.col("call").alias("cy"),
                                  pl.col("log2FoldChange").alias("ly"))  # fmt: skip
    both = a.join(b, on="gene_id").drop_nulls()
    rows = both.select("cx", "cy").rows()
    kinds = {"same": 0, "x_only": 0, "y_only": 0, "opposite": 0, "neither": 0}
    for cx, cy in rows:
        if cx != "ns" and cy != "ns":
            kinds["same" if cx == cy else "opposite"] += 1
        elif cx != "ns":
            kinds["x_only"] += 1
        elif cy != "ns":
            kinds["y_only"] += 1
        else:
            kinds["neither"] += 1
    return kinds


@pytest.mark.parametrize("size", ["wide_light", "phone_dark"])
def test_the_agreement_view_counts_each_kind_of_call(m5, size):
    """Mutations caught: x-only and y-only swapped; opposite calls counted as the same; the
    counts left out of the legend; the build's colours kept."""
    folder, _, seen = m5
    views = _views(seen[size], "Colour by")
    pairs = [("a_naive", "a_adjusted"), ("a_case", "b_case")]
    agreeing = [v for v in views if any("same way in both" in (t["name"] or "")
                                        for t in v["options"][1]["traces"])]  # fmt: skip
    assert len(agreeing) == len(pairs)
    for view, (x, y) in zip(agreeing, pairs, strict=True):
        expected = _agreement(folder, x, y)
        names = [t["name"] for t in view["options"][1]["traces"]]
        for kind, phrase in (("same", "same way in both"),
                             ("opposite", "opposite ways"),
                             ("neither", "neither")):  # fmt: skip
            assert f"{phrase} ({expected[kind]:,})" in names, (kind, names)
        only = [n for n in names if "contrast only (" in n]
        assert [int(n.rsplit("(", 1)[1].rstrip(")").replace(",", "")) for n in only] == [
            expected["x_only"],
            expected["y_only"],
        ]
        assert view["options"][0]["traces"] == view["back"]["traces"]
        assert "rho = " in view["options"][1]["caption"]


@pytest.mark.parametrize("size", ["wide_light", "phone_dark"])
def test_the_boxes_carry_deseq2s_padj_over_every_tested_gene(m5, size):
    """Mutations caught: the Mann-Whitney p kept in the DESeq2 view; a gene's padj taken from
    another gene; the family's size left out of the caption."""
    folder, _, seen = m5
    (view,) = _views(seen[size], "p above each gene")
    results = _table(folder, "a_case_results")
    tested = int(results["padj"].is_not_null().sum())
    box = pl.read_csv(next(folder.glob("*_de_desynth_box_a_case_v1.csv")))
    genes = box.select("feature").unique(maintain_order=True)["feature"].to_list()
    padj = {f: p for f, p in box.select("feature", "deseq2_padj").unique().rows()}
    notes = view["options"][1]["notes"]
    assert len(notes) == len(genes)
    for gene, note in zip(genes, notes, strict=True):
        assert note.startswith(f"padj={padj[gene]:.2g}<br>lfc="), (gene, note)
    assert f"over all {tested:,} tested genes" in view["options"][1]["caption"]
    assert all(n.startswith("p=") for n in view["options"][0]["notes"])
    if size == "phone_dark":
        # A phone turns each gene's p into its tick label and hides the labels over the
        # boxes: the view's own phone layout, not the build's, after the view is drawn.
        deseq, mwu = view["options"][1], view["options"][0]
        assert deseq["ticks"] and all("padj=" in t for t in deseq["ticks"]), deseq["ticks"]
        assert mwu["ticks"] and all("p=" in t and "padj" not in t for t in mwu["ticks"])
        assert deseq["hidden"] == len(genes) == mwu["hidden"]
    else:
        # Upright, each tick names its gene alone: the p sits over its box.
        assert all(o["ticks"] == genes and o["hidden"] == 0 for o in view["options"])


@pytest.mark.parametrize("size", ["wide_light", "phone_dark"])
def test_each_view_draws_its_python_figure_and_the_page_fits(m5, size):
    """Mutations caught: a view drawn from the build's traces; a select wider than a phone; the
    dark theme not applied at load; a grid's selects loaded off the build's setting."""
    _, _, seen = m5
    page, errors = seen[size]["acted"], seen[size]["errors"]
    assert errors == []
    assert page["scroll"] <= page["box"], (page["scroll"], page["box"])
    assert page["theme"] == ("black" if size == "phone_dark" else "white")
    for view in page["views"]:
        assert not view["loaded"]["moved"] and view["loaded"]["line"] == ""
        assert view["loaded"]["traces"] == view["options"][0]["traces"]
        for drawn, built in zip(view["options"], view["built"], strict=True):
            assert [(t["name"], t["y"]) for t in drawn["traces"]] == [
                (t["name"], t["y"]) for t in built
            ]
        assert all(o["wide"] <= page["box"] for o in view["options"])


def test_a_grid_setting_exports_its_own_caption(m5, tmp_path):
    """Mutations caught: the build's caption printed under another setting's bars; the export
    switch left off for a grid."""
    _, page, seen = m5
    (grid,) = _views(seen["wide_light"], "Calls by")
    option = next(i for i, o in enumerate(grid["options"]) if "Holm-adjusted p below 0.01" in
                  o["caption"] and "±2)" in o["caption"])  # fmt: skip
    result = exported(page, tmp_path, f"#{grid['figure']}", option)
    assert "Holm-adjusted p below 0.01, log2 FC at or beyond ±2" in result["acted"]
    svg = (tmp_path / result["files"]["svg"]["guid"]).read_text()
    text = " ".join(re.findall(r">([^<>]+)</text>", svg))
    assert "Holm-adjusted p below 0.01" in text and "BH-adjusted" not in text
    assert (tmp_path / result["files"]["pdf"]["guid"]).read_bytes().startswith(b"%PDF")


def test_the_deseq2_boxes_export_their_own_labels_and_caption(m5, tmp_path):
    """Mutations caught: a box figure with no export, whose buttons kept the build's files under
    the DESeq2 view; the build's Mann-Whitney caption printed under DESeq2's labels."""
    _, page, seen = m5
    (view,) = _views(seen["wide_light"], "p above each gene")
    result = exported(page, tmp_path, f"#{view['figure']}", 1)
    assert result["acted"].startswith("DESeq2 Wald test; padj is BH-adjusted over all")
    svg = (tmp_path / result["files"]["svg"]["guid"]).read_text()
    text = " ".join(re.findall(r">([^<>]+)</text>", svg))
    assert "padj=" in svg and "p=0" not in svg.replace("padj=", "")
    assert "DESeq2 Wald test" in text and "Mann-Whitney U" not in text
