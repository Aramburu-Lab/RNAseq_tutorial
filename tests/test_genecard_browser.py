"""E.1a in headless Chromium: M9's F3 and F4 recompute in the page as Python does, and a list
from a link reaches the page through its folder's entry page.

The made-up root's NPTX2 card carries the real HPA v25 immune panel and a two-tissue exon
panel. Each test names the mutation it catches; each was planted and seen to fail. Skipped
where no Chromium is found; set ``CTX_CHROMIUM``.
"""

import polars as pl
import pytest
from genecard_root import CALLED, SLUG, Stub, make, responses
from test_export_browser import CHROMIUM, _browse
from test_genecard import _cards, _register, _settings

from ctx.genecard import annotation, figures, hpa
from ctx.genecard import run as module
from ctx.report import genes
from ctx.viz.heatmap import diverging_levels

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
NPTX2, ACTB, UPA, TWIN = 106236, 75624, 900011, 900001
# Untick the first box of each checkbox set (and, for F3, take tau on nTPM), then read what
# the page recomputed and drew. The theme toggle is pressed first when the argument is true.
MOVE = """(async () => {
  if (%s) document.getElementById("ctx-theme-toggle").click();
  await new Promise(r => setTimeout(r, 1500));
  const out = {};
  for (const form of document.querySelectorAll("form.subset")) {
    const fig = form.closest(".fig"), cfg = JSON.parse(document.getElementById(
      form.getAttribute("data-controls")).textContent);
    const box = form.querySelector("input[name=kept]");
    box.checked = false;
    if (form.elements.scale) form.elements.scale.value = "linear";
    form.dispatchEvent(new Event("change"));
    const d = fig.querySelector(".js-plotly-plot");
    out[cfg.twin] = {dropped: box.value, got: fig.ctxM9, options: cfg.options,
      lead: fig.querySelector(".fig-lead").textContent,
      summary: fig.querySelector(".subset-summary").textContent,
      visible: d.data.map(t => t.visible === undefined ? true : t.visible),
      categories: (d.layout.yaxis.categoryarray || []).slice(),
      z: d.data[0].z, ticks: d.layout.yaxis2 ? d.layout.yaxis2.ticktext : null,
      buttons: Array.from(fig.querySelectorAll("a[data-export]")).map(a => a.hasAttribute("href"))};
  }
  out.theme = document.documentElement.getAttribute("data-theme");
  out.scroll = document.documentElement.scrollWidth <= document.documentElement.clientWidth;
  return out;
})()"""
# What the list panel shows once the page has read its list.
LISTED = """(async () => {
  let panel = null;
  for (let i = 0; i < 300 && !(panel = document.querySelector(".genelist-panel")); i++) {
    await new Promise(r => setTimeout(r, 100));
  }
  if (!panel) return {url: location.href, list: null};
  for (let i = 0; i < 100 && !panel.ctxList; i++) await new Promise(r => setTimeout(r, 100));
  return {url: location.href, list: panel.ctxList || null,
          summary: panel.querySelector(".list-summary").textContent,
          rows: Array.from(panel.querySelectorAll("tbody tr")).map(tr =>
            Array.from(tr.children).map(td => td.textContent)),
          links: Array.from(panel.querySelectorAll("tbody a")).map(a => a.getAttribute("href")),
          command: panel.querySelector(".list-command code").textContent,
          hidden: panel.querySelector(".list-command").hidden};
})()"""


@pytest.fixture(scope="module")
def ran(tmp_path_factory):
    """The made-up root's page with NPTX2 and ACTB cards, driven at both widths and themes."""
    folder = tmp_path_factory.mktemp("m9") / "root"
    spec = make(folder)
    _register(folder, spec)
    options = module.Options(spec, genes=("NPTX2", "ACTB"))
    written = module.run(_settings(folder), _cards(), "ctx genecard", SLUG, options,
                         Stub(responses()))  # fmt: skip
    page = written[0]
    entry = page.parent / "index.html"
    listed = genes.Selection(f"{SLUG}:a_case", (ACTB, NPTX2, TWIN, UPA))
    tasks = [{"name": f"move_{w}_{dark}", "url": page.as_uri(), "width": w, "height": 900,
              "act": MOVE % ("true" if dark else "false")}
             for w in (1280, 375) for dark in (False, True)]  # fmt: skip
    tasks.append({"name": "link", "url": entry.as_uri() + genes.encode(listed),
                  "until": "_overview_v1.html", "act": LISTED})  # fmt: skip
    return folder, _browse(tasks, tmp_path_factory.mktemp("downloads"))


def _immune_rows() -> list[tuple[str, str, float]]:
    """NPTX2's HPA immune cells, the pooled total left out, as F3 draws them."""
    xml = hpa.parse_xml(responses()[hpa.url("v25", "ENSG00000106236", "xml")], "v25")
    return [(c, lineage, v) for c, lineage, v in xml.immune if c != hpa.POOLED]


def test_f3_retakes_tau_and_the_top_cell_over_the_lineages_kept(ran):
    """With the first lineage left out and tau on nTPM, the page's tau, top cell type and
    count equal ``figures.immune_subset``; that lineage's bars are hidden and its cell types
    leave the axis; the PDF and SVG buttons draw the view. Mutations caught: the scale
    ignored; the lineage's trace left visible; the axis keeping its cell types; the export
    left at the build's files."""
    _, results = ran
    rows = _immune_rows()
    lineages = list(dict.fromkeys(r[1] for r in rows))
    for name, r in results.items():
        if not name.startswith("move"):
            continue
        assert r["errors"] == [], name
        page = r["acted"]["immune"]
        kept = set(lineages) - {page["dropped"]}
        want = figures.immune_subset(rows, kept, log=False)
        assert page["got"]["cells"] == want["cells"] and page["got"]["top"] == want["top"]
        # NPTX2's immune cells are all under 1 nTPM, so tau is undefined on both sides; the
        # node twins (test_genecard_twins) compare defined ones.
        assert page["got"]["tau"] is None and want["tau"] is None
        assert f"tau {figures.shown(want['tau'])}" in page["lead"] and want["top"] in page["lead"]
        assert page["visible"] == [lineage in kept for lineage in lineages]
        assert page["categories"] == [c for c, lineage, _ in rows if lineage in kept]
        assert page["summary"].endswith("; tau on nTPM.")
        assert page["summary"].startswith(f"{want['cells']} of {len(rows)} cell types")
        assert page["buttons"] and not any(page["buttons"])
        assert r["acted"]["scroll"], name
    assert {r["acted"]["theme"] for r in results.values() if "acted" in r and r["acted"]
            and "theme" in r["acted"]} == {"white", "black"}  # fmt: skip


def test_f4_retakes_the_exon_panel_over_the_tissues_kept(ran):
    """With the first tissue left out, each kept cell's level and each exon's CV equal
    ``annotation.exon_use`` over the rest, the left-out column draws nothing, and the lead
    falls back to no pair over one tissue. Mutations caught: the reference over every tissue;
    the left-out column still drawn; the CV axis left as built."""
    folder, results = ran
    out = folder / "outputs" / "11_genecard"
    shares = pl.read_parquet(next(out.glob("*_exons_v1.parquet"))).filter(
        pl.col("symbol") == "NPTX2"
    )
    names = [s.level for s in diverging_levels([])[1]]
    for name, r in results.items():
        if not name.startswith("move"):
            continue
        page = r["acted"]["exons"]
        # The panel's own tissues, those at 1 TPM or more; the table holds every tissue.
        tissues = page["options"]
        assert set(tissues) < set(shares["tissue"])
        kept = [t for t in tissues if t != page["dropped"]]
        use = annotation.exon_use(shares.select("exon", "tissue", "share"), kept)
        levels, _ = diverging_levels(list(use.steps))
        want = {(e, t): lv for (e, t), lv in zip(use.frame.select("exon", "tissue").rows(),
                                                  levels, strict=True)}  # fmt: skip
        exons = list(dict.fromkeys(shares["exon"].to_list()))
        for i, e in enumerate(exons):
            for j, t in enumerate(tissues):
                cell = page["z"][i][j]
                if t == page["dropped"]:
                    assert cell is None
                else:
                    assert (names[cell] if cell is not None else None) == want[(e, t)]
        assert page["ticks"] == ["" if use.cv[e] is None else f"{use.cv[e]:.2f}" for e in exons]
        pair = annotation.qpcr_pair(use.cv)
        assert (page["got"]["pair"] is None) is (pair is None)
        assert f"{len(kept)} of {len(tissues)} tissues" in page["summary"]


def test_a_list_from_a_link_opens_the_page_through_its_entry_with_each_gene(ran):
    """The entry page forwards the fragment; the page counts the list against every tested
    gene: UPA with its evidence, ACTB tested and carded, NPTX2 untested but linked to its
    card, TWIN neither and not
    shown; the command builds UPA's card alone. Mutations caught: the fragment dropped by the
    entry page; a gene neither tested nor carded shown as a row; a carded gene M5 did not test
    left out; the command naming a gene that has a card."""
    _, results = ran
    got = results["link"]["acted"]
    assert results["link"]["errors"] == []
    assert "_overview_v1.html#v=1&" in got["url"]
    assert got["list"] == {"found": [ACTB, UPA], "missing": 2, "cards": 2}
    _, p, lfc, call = CALLED["ENSG00000900011.1"]
    assert (p, lfc, call) == (1e-9, 3.0, "up")  # as genecard_root plants them
    assert got["rows"] == [["UPA", "3.00", "1.0e-9", "up", "none here"],
                           ["ACTB", "0.400", "0.20", "ns", "Card"],
                           ["NPTX2", "", "", "not tested", "Card"]]  # fmt: skip
    assert got["links"] == ["#group-3", "#group-2"]
    assert got["command"].endswith(" --genes ENSG00000900011") and not got["hidden"]
    assert "4 genes from a link: 2 tested by M5, 2 not; 2 with a card" in got["summary"]
