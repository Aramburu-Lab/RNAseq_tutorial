"""M8 ``ctx ncrna``: biotype groups, nearest neighbours, Spearman and its page twin, the page.

Each test names the mutation it catches. ``ctxPairs`` (templates/_pairs.js.j2) runs under node
on the same cases as ``ctx.ncrna.neighbours``.
"""

import inspect
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from scipy.stats import false_discovery_control, spearmanr

from ctx.ncrna import NcrnaError, biotypes, figures, neighbours

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src" / "ctx" / "report" / "templates" / "_pairs.js.j2"
CTX = Path(sys.executable).parent / "ctx"


def test_group_decides_rrna_and_small_rna_before_pseudogene():
    """Catches: the pseudogene rule first (rRNA_pseudogene, 5S copies, read as pseudogene), the
    mitochondrial rRNA pooled with the small RNAs, or an unlisted type left without a group."""
    types = ["rRNA_pseudogene", "rRNA", "Mt_rRNA", "miRNA", "snoRNA", "IG_V_pseudogene",
             "transcribed_unprocessed_pseudogene", "protein_coding", "lncRNA", "Mt_tRNA",
             "TEC"]  # fmt: skip
    got = pl.DataFrame({"t": types}).select(g=biotypes.group(pl.col("t")))["g"].to_list()
    assert got == ["small_RNA", "small_RNA", "Mt_rRNA", "small_RNA", "small_RNA", "pseudogene",
                   "pseudogene", "protein_coding", "lncRNA", "other", "other"]  # fmt: skip


def test_annotate_keeps_count_order_and_refuses_unknown_genes():
    """Catches: rows returned in the annotation's order (every share then lands on the wrong
    gene), or a gene missing from the annotation passed through untyped."""
    types = pl.DataFrame({"gene_id": ["a", "b", "c"], "gene_name": ["A", "B", "C"],
                          "gene_type": ["lncRNA", "protein_coding", "miRNA"]})  # fmt: skip
    got = biotypes.annotate(pl.Series(["c", "a", "b"]), types)
    assert got["gene_id"].to_list() == ["c", "a", "b"]
    assert got["biotype"].to_list() == ["small_RNA", "lncRNA", "protein_coding"]
    with pytest.raises(NcrnaError, match="1 genes are not in the annotation"):
        biotypes.annotate(pl.Series(["a", "z"]), types)


def test_shares_by_hand():
    """Catches: shares over the group's own total instead of the sample's, or genes by samples
    transposed."""
    matrix = np.array([[90, 10], [5, 30], [5, 60]])
    got = biotypes.shares(matrix, np.array(["x", "y", "y"]), ("x", "y"))
    assert got.tolist() == [[90.0, 10.0], [10.0, 90.0]]


def test_parse_biotypes():
    """Catches: an unknown group accepted, or the groups kept in typed order, not GROUPS order."""
    assert biotypes.parse(("lncRNA", "protein_coding")) == ("protein_coding", "lncRNA")
    with pytest.raises(NcrnaError, match="small RNA"):
        biotypes.parse(("small RNA",))


def loci(rows: list[tuple]) -> pl.DataFrame:
    return pl.DataFrame(rows, schema=list(neighbours.LOCUS), orient="row")


def test_nearest_kinds_at_their_boundaries():
    """Each lncRNA beside one protein-coding gene, at the edges of each kind.

    Catches: genes that touch counted as overlapping (gap 0 is not overlap), the gap off by one
    (10,000 nt between is near, 10,001 far), the strand comparison inverted, a lncRNA on a
    chromosome without protein-coding genes given a pair, or a tie not going to the gene that
    starts first.
    """
    coding = loci([("P1", "chr1", 1001, 2000, "+"), ("P2", "chr1", 50_001, 51_000, "-"),
                   ("P3", "chr1", 100_001, 101_000, "+"), ("P4", "chr1", 200_001, 201_000, "+"),
                   ("P5", "chr1", 300_001, 301_000, "+"),
                   ("P6", "chr1", 304_000, 305_000, "+")])  # fmt: skip
    lnc = loci([
        ("antisense", "chr1", 1500, 3000, "-"),
        ("same", "chr1", 50_500, 50_600, "-"),
        ("touching", "chr1", 101_001, 102_000, "-"),  # starts the base after P3 ends
        ("near", "chr1", 211_001, 212_000, "+"),  # 10,000 nt after P4
        ("before", "chr1", 190_000, 190_000, "+"),  # 10,000 nt before P4, 88,999 after P3
        ("lonely", "chr9", 1, 100, "+"),
        ("tie", "chr1", 302_001, 302_999, "+"),  # 1,000 nt after P5 and before P6
    ])  # fmt: skip
    got = {r[0]: r[1:] for r in neighbours.nearest(lnc, coding).rows()}
    assert got["antisense"] == ("P1", 0, "antisense")
    assert got["same"] == ("P2", 0, "same_strand")
    assert got["touching"] == ("P3", 0, "near")
    assert got["near"] == ("P4", 10_000, "near")
    assert got["before"] == ("P4", 10_000, "near")
    assert got["tie"] == ("P5", 1000, "near")
    assert "lonely" not in got
    farther = loci([("far2", "chr1", 211_002, 211_100, "+")])
    assert neighbours.nearest(farther, coding).row(0)[1:] == ("P4", 10_001, "far")


def test_background_draws_from_other_chromosomes_with_a_seed():
    """Catches: a background partner on the lncRNA's own chromosome, or an unseeded draw."""
    coding = loci([(f"P{i}", f"chr{i % 3}", 1, 10, "+") for i in range(30)])
    lnc = loci([(f"L{i}", f"chr{i % 3}", 100, 200, "+") for i in range(30)])
    pairs = pl.DataFrame({"lncrna": lnc["gene_id"]})
    first = neighbours.background(pairs, coding, lnc)
    chrom = dict(coding.select("gene_id", "chrom").rows()) | dict(
        lnc.select("gene_id", "chrom").rows()
    )
    assert all(chrom[a] != chrom[b] for a, b in first.select("lncrna", "neighbour").rows())
    assert first.equals(neighbours.background(pairs, coding, lnc))


TIES = np.array(
    [[1.0, 2.0, 2.0, 3.0, 5.0, 5.0, 5.0, 0.5], [3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0]]
)
OTHER = np.array(
    [[2.0, 1.0, 3.0, 3.0, 4.0, 6.0, 5.0, 1.0], [2.0, 7.0, 1.0, 8.0, 2.0, 8.0, 1.0, 8.0]]
)
STRATA = np.array(["a", "a", "b", "b", "a", "b", "a", "b"])


def test_spearman_and_p_equal_scipy_with_ties():
    """Catches: ties ranked by position rather than at their mean rank, or a p from the normal
    approximation instead of scipy's t."""
    rho = neighbours.spearman(TIES, OTHER)
    p = neighbours.p_value(rho, TIES.shape[1])
    for i in range(2):
        want = spearmanr(TIES[i], OTHER[i])
        assert rho[i] == pytest.approx(want.statistic, abs=1e-12)
        assert p[i] == pytest.approx(want.pvalue, rel=1e-9)


def test_residuals_take_each_stratum_mean():
    """Catches: one mean over every sample removed, or the strata's means swapped."""
    got = neighbours.residuals(np.array([[1.0, 3.0, 10.0, 20.0]]), np.array(["a", "a", "b", "b"]))
    assert got.tolist() == [[-1.0, 1.0, -5.0, 5.0]]


def test_correlate_residuals_remove_a_shared_condition_effect():
    """Two genes that move only with the condition: rho on values is high, on residuals not.

    Catches: the values used where the residuals belong (rho_within equals rho), or BH pooled
    across the values and the residuals.
    """
    rng = np.random.default_rng(1)
    condition = np.repeat([0.0, 3.0], 20)
    a = condition + rng.normal(0, 0.3, (50, 40))
    b = condition + rng.normal(0, 0.3, (50, 40))
    got = neighbours.correlate(a, b, np.repeat(["ref", "test"], 20))
    assert np.median(got["rho"]) > 0.5 > abs(np.median(got["rho_within"]))
    assert np.allclose(got["q"], false_discovery_control(got["p"]))


HARNESS = r"""
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
console.log(JSON.stringify(input.cases.map(([a, b, strata]) => {
  const ra = ctxPairs.residuals(a, strata), rb = ctxPairs.residuals(b, strata);
  return {rho: ctxPairs.spearman(a, b), within: ctxPairs.spearman(ra, rb), ra: ra,
          ranks: ctxPairs.rank(a)};
})));
"""


def test_page_twin_equals_python():
    """ctxPairs under node against ctx.ncrna.neighbours: rho on values and residuals, ranks.

    Catches, in the page script: ties ranked by position, the residuals taken over all
    samples, a constant side not given null (Python's NaN), or the ranks counted from 0.
    """
    node = shutil.which("node")
    assert node, "nodejs is pinned in the ctx-dev environment; run tests with pixi"
    script = SCRIPT.read_text().replace("{% raw %}", "").replace("{% endraw %}", "")
    constant = np.full(8, 2.0)
    cases = [(TIES[0], OTHER[0]), (TIES[1], OTHER[1]), (constant, OTHER[0])]
    payload = {"cases": [[a.tolist(), b.tolist(), STRATA.tolist()] for a, b in cases]}
    done = subprocess.run([node, "-e", script + HARNESS], input=json.dumps(payload),
                          capture_output=True, text=True, check=True)  # fmt: skip
    page = json.loads(done.stdout)
    for (a, b), js in zip(cases, page, strict=True):
        rho = neighbours.spearman(a[None], b[None])[0]
        within = neighbours.spearman(neighbours.residuals(a[None], STRATA),
                                     neighbours.residuals(b[None], STRATA))[0]  # fmt: skip
        for got, want in ((js["rho"], rho), (js["within"], within)):
            assert (got is None) if np.isnan(want) else got == pytest.approx(want, abs=1e-12)
        assert js["ra"] == pytest.approx(neighbours.residuals(a[None], STRATA)[0].tolist())
    assert page[0]["ranks"] == [2.0, 3.5, 3.5, 5.0, 7.0, 7.0, 7.0, 1.0]


def test_by_biotype_adjusts_within_each_biotype_over_tested_genes():
    """Catches: M5's pooled padj reused for the call, an untested gene (no padj) adjusted with
    the rest, or a biotype outside --biotypes kept."""
    results = pl.DataFrame({
        "gene_id": ["c1", "c2", "c3", "l1", "l2", "x1", "u1"], "gene_name": list("ABCDEFG"),
        "baseMean": [1.0] * 7, "log2FoldChange": [2.0] * 7,
        "pvalue": [0.001, 0.01, 0.02, 0.02, 0.04, 0.001, 0.5],
        "padj": [0.007, 0.035, 0.035, 0.035, 0.047, 0.007, None],
    })  # fmt: skip
    genes = pl.DataFrame({"gene_id": ["c1", "c2", "c3", "l1", "l2", "x1", "u1"],
                          "biotype": ["protein_coding"] * 3 + ["lncRNA"] * 2
                          + ["pseudogene", "lncRNA"]})  # fmt: skip
    got = figures.by_biotype(results, genes, ("protein_coding", "lncRNA"))
    assert got["gene_id"].to_list() == ["c1", "c2", "c3", "l1", "l2"]
    within = dict(got.select("gene_id", "padj_within").rows())
    assert within["l1"] == pytest.approx(0.04) and within["l2"] == pytest.approx(0.04)
    assert within["c3"] == pytest.approx(0.02)
    assert got.filter(pl.col("gene_id") == "l2")["call"].item() == "up"
    assert got.filter(pl.col("gene_id") == "l2")["call_pooled"].item() == "up"


def test_composition_mirna_share_is_the_hairpin_genes_alone():
    """Catches: the miRNA share taken over the small_RNA group (snRNA counted as miRNA)."""
    sheet = pl.DataFrame({"sample": ["s1"], "cell_type": ["A"], "condition": ["c"],
                          "library_selection": ["Oligo-dT"]})  # fmt: skip
    genes = pl.DataFrame({"gene_type": ["miRNA", "snRNA", "protein_coding"],
                          "biotype": ["small_RNA", "small_RNA", "protein_coding"]})  # fmt: skip
    _, wide = figures.composition(sheet, genes, np.array([[1], [9], [990]]), "caption")
    assert wide["mirna_percent"].item() == pytest.approx(0.1)
    assert wide["small_RNA_percent"].item() == pytest.approx(1.0)


def test_composition_filters_by_biotype_too():
    """F1's subset by biotype (M8 plan, F1's control): the page offers a biotype select holding
    every group, beside cell type and condition, and a GDC sheet's bars filter by it as well.
    Catches: the biotype filter left off, which leaves only the sheet's columns to filter by."""
    genes = pl.DataFrame({"gene_type": ["miRNA", "lncRNA", "protein_coding"],
                          "biotype": ["small_RNA", "lncRNA", "protein_coding"]})  # fmt: skip
    counts = np.array([[1, 2], [9, 8], [990, 990]])
    step1 = pl.DataFrame({"sample": ["s1", "s2"], "cell_type": ["A", "B"], "condition": ["c", "d"],
                          "library_selection": ["Oligo-dT", "Oligo-dT"]})  # fmt: skip
    gdc = pl.DataFrame({"sample": ["s1", "s2"], "condition": ["tumour", "normal"]})
    for sheet, facets in ((step1, ["cell_type", "condition"]), (gdc, ["condition"])):
        bars, _ = figures.composition(sheet, genes, counts, "caption")
        controls = bars.controls()
        assert controls is not None
        assert list(controls["filters"]) == [*facets, "biotype"]
        assert controls["filters"]["biotype"] == sorted(biotypes.GROUPS)


def test_shown_pairs_are_cis_pairs_by_q_then_strength():
    """Catches: background pairs offered in the select, a pair without q (a constant gene)
    offered, or a tie in q broken by signed rho, which puts the strongest negative pair last."""
    pairs = pl.DataFrame({
        "lncrna_id": ["bg", "nan", "neg", "pos", "weak"],
        "kind": ["trans", "near", "antisense", "near", "far"],
        "q_within": [1e-9, float("nan"), 1e-4, 1e-4, 0.3],
        "rho_within": [0.99, float("nan"), -0.9, 0.8, 0.1],
    })  # fmt: skip
    assert figures.shown_pairs(pairs, 4)["lncrna_id"].to_list() == ["neg", "pos", "weak"]


def test_cli_defaults_are_the_module_defaults():
    """Catches: the CLI's typed defaults drifting from biotypes.TESTED and M5's MIN_COUNT."""
    from ctx.cli import ncrna  # noqa: PLC0415
    from ctx.de.filters import MIN_COUNT  # noqa: PLC0415

    defaults = {p.name: p.default for p in inspect.signature(ncrna).parameters.values()}
    assert tuple(defaults["biotypes"].split(",")) == biotypes.TESTED
    assert defaults["min_count"] == MIN_COUNT


@pytest.fixture(scope="session")
def built(ncrna_built):
    """ctx ncrna on the made-up root, the session's run: (output folder, paths printed)."""
    root, _, printed = ncrna_built
    return root / "outputs" / "08_ncrna", printed


def table(out: Path, what: str) -> pl.DataFrame:
    return pl.read_parquet(next(out.glob(f"*_ncrna_{'ncsynth'}_{what}_v1.parquet")))


def test_page_end_to_end(built):
    """The page, its figures and tables; the filter at M5's cutoff passes M5's tested genes.

    Catches: a consumer reading another filter than M5's (filter.kept differs from m5.tested),
    rule N1 decided on a missing metric, the planted pair not first in the select, or a page
    over the 15 MB cap.
    """
    out, written = built
    page = Path(written[0])
    assert page.name.endswith("_ncrna_ncsynth_overview_v1.html")
    assert page.stat().st_size < 15 * 10**6
    metrics = dict(table(out, "metrics").select("metric", "value").rows())
    assert metrics["filter.kept"] == metrics["m5.tested"] > 0
    assert table(out, "sanity")["verdict"].to_list() == ["PASS"]
    pairs = table(out, "pairs")
    kinds = set(pairs["kind"])
    assert {"antisense", "same_strand", "near", "far", "trans"} <= kinds
    html = page.read_text()
    blocks = re.findall(r'<script type="application/json" id="ctx-controls-\d+">(.*?)</script>',
                        html, flags=re.S)  # fmt: skip
    (controls,) = [c for c in map(json.loads, blocks) if c["kind"] == "pairs"]
    assert controls["options"][0].startswith("G00201 and G00001")
    assert len(controls["options"]) == 5


def test_the_filter_figure_offers_each_preset_as_m5_counts_it(built):
    """M8 F2's control: the run's minimum count first, then the plan's presets, each kept
    share M5's own filter at that count; a higher count keeps no more genes.

    Mutations caught: the presets' kept shares all taken at the run's count; the run's count
    not the first view, so the page opens on another figure than its files."""
    out, written = built
    html = Path(written[0]).read_text()
    blocks = re.findall(r'<script type="application/json" id="ctx-controls-\d+">(.*?)</script>',
                        html, flags=re.S)  # fmt: skip
    (view,) = [c for c in map(json.loads, blocks) if c["kind"] == "view"]
    assert [o["label"] for o in view["options"]] == [
        "10 reads", "1 read", "2 reads", "5 reads", "20 reads", "50 reads"]  # fmt: skip
    presets = table(out, "filter_presets")
    for biotype in set(presets["biotype"]):
        shares = presets.filter(pl.col("biotype") == biotype).sort("min_count")["kept_pct"]
        assert shares.to_list() == sorted(shares.to_list(), reverse=True)
        assert shares.n_unique() > 1
    metrics = dict(table(out, "metrics").select("metric", "value").rows())
    at_ten = presets.filter(pl.col("min_count") == 10)  # noqa: PLR2004 - the run's own
    assert len(at_ten) and metrics["filter.kept"] > 0


def test_refuses_without_ncrna(tmp_path, built):
    """Catches: a dataset file without [ncrna] run with some default contrast."""
    out, _ = built
    root = out.parents[1]
    spec = tmp_path / "plain.toml"
    text = next(root.parent.glob("ncsynth.toml")).read_text()
    spec.write_text(text.split("[ncrna]")[0])
    done = subprocess.run([str(CTX), "--root", str(root), "ncrna", "--dataset", "ncsynth",
                           "--dataset-file", str(spec)], capture_output=True, text=True,
                          check=False)  # fmt: skip
    assert done.returncode == 1 and "needs [ncrna]" in done.stderr
