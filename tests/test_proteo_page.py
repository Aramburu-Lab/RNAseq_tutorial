"""``ctx proteo`` end to end on made-up roots (tests/proteo_root.py), both designs, and the
rules' registration order. Each test names the mutation it catches.
"""

import html as markup
import re
import subprocess
from pathlib import Path

import polars as pl
import pytest
from proteo_root import TUMOURS

from ctx.report import text
from ctx.report.builder import unpacked
from ctx.report.cap import MAX_BYTES

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def across_root(proteo_across):
    """The across-samples root, ctx proteo run."""
    return proteo_across[0]


@pytest.fixture(scope="session")
def cross_root(proteo_cross):
    """The cross-dataset root, ctx proteo run."""
    return proteo_cross[0]


def output(root: Path, what: str) -> Path:
    """The one proteo output of this kind."""
    (path,) = (root / "outputs" / "07_proteo").glob(f"*_proteo_*_{what}")
    return path


@pytest.mark.parametrize("which", ["across_root", "cross_root"])
def test_page_offline_and_every_table_described(which, request):
    """Catches a page over the cap or fetching a script, and a table without its README."""
    root = request.getfixturevalue(which)
    html = unpacked(output(root, "overview_v1.html").read_text())
    assert len(html.encode()) <= MAX_BYTES
    assert not re.search(r"<script[^>]+\bsrc\s*=", html, re.IGNORECASE)
    for table in (root / "outputs" / "07_proteo").glob("*.csv"):
        assert table.with_name(table.stem + "_README.md").is_file(), table.name


@pytest.mark.parametrize("which", ["across_root", "cross_root"])
def test_missingness_first_and_rules_pass(which, request):
    """Catches the class or fold figure drawn before the missingness figures (K3 fails), and
    the planted class shift or fold agreement lost (K2, P1 or P2 fails)."""
    root = request.getfixturevalue(which)
    verdicts = pl.read_parquet(output(root, "sanity_v1.parquet"))
    assert set(verdicts["rule_verdict"]) == {"PASS"}, verdicts.select("rule", "left").rows()
    html = unpacked(output(root, "overview_v1.html").read_text())
    # A figure's Reproduce panel names its data table again, after its own link.
    stems = list(dict.fromkeys(re.findall(r"_proteo_fake_\w+?_(\w+?)_v1\.csv", html)))
    assert stems[:2] == ["quantified", "coverage"], stems


def test_across_missing_never_zero_and_one_row_per_gene(across_root):
    """Catches a missing protein read as 0 (G10 would count every sample), the extra G1 row
    kept as a second gene, and samples joined by position instead of case."""
    genes = pl.read_parquet(output(across_root, "genes_v1.parquet"))
    assert genes["gene"].n_unique() == genes.height
    g10 = genes.filter(pl.col("gene") == "G10").row(0, named=True)
    assert g10["n_quantified"] == TUMOURS // 2 and not g10["tested"]
    g9 = genes.filter(pl.col("gene") == "G9").row(0, named=True)
    assert g9["rho"] > 0.95
    metrics = pl.read_parquet(output(across_root, "metrics_v1.parquet"))
    assert metrics["metric"].n_unique() == metrics.height


def test_across_f4_f5_and_the_sample_id_column(across_root):
    """F4 compares M5's calls with protein on the same tumours; F5 pairs each gene's rows.

    Mutations caught: the column named by GDC sample ID left out (29 of 30 paired); F4 not
    drawn when M5's results exist; G11's planted protein shift not called (both is 0); G10,
    untested at protein, counted as an RNA-only call; F5 a gene's protein row dropped."""
    html = unpacked(output(across_root, "overview_v1.html").read_text())
    assert f"{TUMOURS} of {TUMOURS} protein samples" in html and "1 of them are named" in html
    metrics = dict(pl.read_parquet(output(across_root, "metrics_v1.parquet"))
                   .select("metric", "value").rows())  # fmt: skip
    assert (metrics["f4.basal_vs_luma.both"], metrics["f4.basal_vs_luma.rna"]) == (1.0, 1.0)
    assert metrics["f4.basal_vs_luma.protein_up"] == 1.0
    f4 = pl.read_csv(output(across_root, "contrast_basal_vs_luma_v1.csv"), comment_prefix="#")
    assert "G10" not in f4["gene"].to_list() and f4.filter(pl.col("called"))["gene"][0] == "G11"
    paired = pl.read_csv(output(across_root, "paired_zscores_v1.csv"), comment_prefix="#")
    rows = paired["row"].unique().to_list()
    assert len(rows) % 2 == 0 and {r.rsplit(" ", 1)[1] for r in rows} == {"RNA", "protein"}
    assert paired["sample"].n_unique() == TUMOURS


def test_cross_untested_below_the_share(cross_root):
    """Catches the 70 % filter skipped: G30 has 5 of 8 values."""
    genes = pl.read_parquet(output(cross_root, "genes_v1.parquet"))
    g30 = genes.filter(pl.col("gene") == "G30").row(0, named=True)
    assert g30["n_samples"] == 5 and not g30["tested"]
    html = unpacked(output(cross_root, "overview_v1.html").read_text())
    assert "Cross-study comparison: RNA from fake_cells" in html


def git(*args: str) -> str:
    """git in this checkout."""
    return subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True,
                          check=True).stdout.strip()  # fmt: skip


@pytest.mark.parametrize("name", ["Non-cancerous-blood", "cptac_brca"])
def test_rules_registered_before_the_protein_reader(name):
    """The registration commit of [proteo.sanity] precedes the commit that added the reader.

    Catches a rule edited after protein values could be read without registering again."""
    if git("rev-parse", "--is-shallow-repository") == "true":
        pytest.skip("shallow clone: the history is absent")
    path = f"conf/datasets/{name}.toml"
    registered = git("log", "--format=%H", "-G", r"^registered_at = \"2026-10", "--", path).split()
    reader = git("log", "--reverse", "--diff-filter=A", "--format=%H", "--",
                 "src/ctx/proteo/tables.py").split()  # fmt: skip
    assert registered and reader
    assert registered[0] != reader[0]
    ancestor = subprocess.run(["git", "-C", str(ROOT), "merge-base", "--is-ancestor",
                               registered[0], reader[0]], check=False)  # fmt: skip
    assert ancestor.returncode == 0


# Written here, not read from ctx.proteo.page.VARIANTS, so a picker that falls back fails.
DESIGNS = {"across_samples": "variants.across_samples", "tumour_normal": "variants.tumour_normal",
           "cross_dataset": None}  # fmt: skip
SECTIONS = ("purpose", "input", "method", "interpretation", "pitfalls", "terms", "questions")


def _plain(marked: str) -> str:
    """Page text as a reader sees it: a glossary mark is its display text, spaces collapsed."""
    shown = re.sub(r"\[\[([^\]|]+?)!?(?:\|([^\]]+))?\]\]", lambda m: m[2] or m[1], marked)
    return " ".join(shown.replace("`", "").split())


def _sentences(design: str) -> set[str]:
    """Every sentence of a design's sections, as the page prints it."""
    made = text.skeleton("proteo", DESIGNS[design])
    parts = [made[k] for k in ("purpose", "input", "method", "interpretation")]
    parts += [*made["pitfalls"], *(m for _, m in made["terms"])]
    parts += [x for q in made["questions"] for x in (q.prompt, q.answer)]
    return {s for part in parts for s in re.split(r"(?<=[.?])\s+", _plain(part)) if s}


def _seen(root: Path) -> str:
    """The overview page's words, tags dropped."""
    html = unpacked(output(root, "overview_v1.html").read_text())
    words = " ".join(markup.unescape(re.sub(r"<[^>]+>", " ", html)).split())
    # A mark's closing tag sits before the punctuation after it; the reader sees none.
    return re.sub(r" ([,.;:?])", r"\1", words)


def test_every_design_writes_every_section():
    """L.1: no design falls back to another design's text: each variant replaces every section."""
    for key in filter(None, DESIGNS.values()):
        assert set(SECTIONS) <= set(text.raw("proteo", key)), key


@pytest.mark.parametrize("design, fixture", [("across_samples", "proteo_across"),
                                             ("cross_dataset", "proteo_cross"),
                                             ("tumour_normal", "proteo_made")])  # fmt: skip
def test_page_text_is_its_own_designs(design, fixture, request):
    """L.1: each design's page shows every sentence of its own sections and none of the
    sentences only another design's sections hold (granules, donors who differ, normals).

    Catches: the variant picker falling back to the page's own, cross-study, sections."""
    seen = _seen(request.getfixturevalue(fixture)[0])
    own = _sentences(design)
    others = set().union(*(_sentences(d) for d in DESIGNS if d != design)) - own
    assert others and not [s for s in others if s in seen]
    assert not [s for s in own if s not in seen]
