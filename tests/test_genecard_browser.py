"""E.1a in headless Chromium: M9's F3 and F4 recompute in the page as Python does, and a list
from a link reaches the page through its folder's entry page.

The made-up root's NPTX2 card carries the real HPA v25 immune panel and a two-tissue exon
panel. Each test names the mutation it catches; each was planted and seen to fail. Skipped
where no Chromium is found; set ``CTX_CHROMIUM``.
"""

import json
import math
import re

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
  await new Promise(r => setTimeout(r, 500));
  return {url: location.href, list: panel.ctxList || null,
          top: panel.closest(".group").getBoundingClientRect().top, window: innerHeight,
          scrolled: scrollY,
          summary: panel.querySelector(".list-summary").textContent,
          rows: Array.from(panel.querySelectorAll("tbody tr")).map(tr =>
            Array.from(tr.children).map(td => td.textContent)),
          links: Array.from(panel.querySelectorAll("tbody a")).map(a => a.getAttribute("href")),
          command: panel.querySelector(".list-command code").textContent,
          hidden: panel.querySelector(".list-command").hidden,
          heat: panel.ctxHeat || null,
          keys: Array.from(panel.querySelector(".list-heat-plot").data || []).slice(1)
            .map(t => t.name),
          caption: panel.querySelector(".list-heat-caption").textContent,
          least: panel.querySelector(".list-heat-plot").style.minWidth,
          shown: !panel.querySelector(".list-heat").hidden};
})()"""
# F4 opened alone (E.1d): untick its first tissue, then read what the plot drew.
ALONE = """(async () => {
  await new Promise(r => setTimeout(r, 1500));
  const form = document.querySelector("form.subset");
  if (!form) return null;
  const fig = form.closest(".fig"), box = form.querySelector("input[name=kept]");
  box.checked = false;
  form.dispatchEvent(new Event("change"));
  const d = fig.querySelector(".js-plotly-plot"), at = Array.from(form.querySelectorAll(
    "input[name=kept]")).indexOf(box);
  return {dropped: box.value, got: fig.ctxM9 ? fig.ctxM9.cells.length : null,
          column: d.data[0].z.map(row => row[at]),
          summary: fig.querySelector(".subset-summary").textContent};
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
    alone = next(p for p in written if "_exons_nptx2_" in p.name and p.suffix == ".html")
    tasks += [{"name": f"alone_{w}", "url": alone.as_uri(), "width": w, "height": 900,
               "act": ALONE} for w in (1280, 375)]  # fmt: skip
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
    ``annotation.tissue_use`` over the rest, the left-out column draws nothing, and the lead
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
        use = annotation.tissue_use(shares.select("exon", "tissue", "share"), kept)
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


def test_a_list_from_a_link_brings_its_panel_into_view(ran):
    """The fragment is no anchor, so the page opened at its top and the list sat screens down
    (PK, 2026-10-20): the panel's group is in the top half of the window. Mutation caught:
    the reveal left out."""
    got = ran[1]["link"]["acted"]
    assert got["scrolled"] > 0 and 0 <= got["top"] < got["window"] / 2, got


def test_the_list_draws_its_evidence_as_the_evidence_view_does(ran):
    """E.1d: the list's genes in list order by every contrast: UPA's up call at 3.0 in the top
    step, ACTB tested and not called grey, NPTX2 untested blank; the key holds the steps drawn;
    the plot is no narrower than its turned labels need. Mutations caught: an untested cell
    drawn grey; the key listing steps not drawn; the rows in another order than the list's;
    the cap below the list; the least width left to plotly."""
    got = results_link = ran[1]["link"]["acted"]
    page = next((ran[0] / "outputs" / "11_genecard").glob("*_overview_v1.html")).read_text()
    heat = json.loads(re.search(r'id="ctx-panel-1">([^<]*)</script>', page)[1])["heat"]
    assert got["least"] == f"{math.ceil(heat['least_px'])}px"
    steps = [s.level for s in diverging_levels([])[1]]
    assert results_link["shown"]
    assert got["heat"]["genes"] == [ACTB, NPTX2, UPA]
    assert got["heat"]["z"] == [[len(steps)], [None], [steps.index("≥ 2")]]
    assert got["keys"] == ["≥ 2", "not called"]
    assert got["caption"].startswith("Columns are M5 contrasts, 1 here: DESeq2 Wald test")
    assert "The first 3 of 3 genes in list order" in got["caption"]


def test_f4_opened_alone_keeps_its_tissue_boxes(ran):
    """E.1d: F4's own page carries its checkbox set, and unticking a tissue empties its column
    there as on the module page, at 1,280 and 375 px. Mutation caught: the page opened alone
    written from the contract, which carries no controls."""
    _, results = ran
    for width in (1280, 375):
        r = results[f"alone_{width}"]
        assert r["errors"] == [] and r["acted"] is not None, width
        assert r["acted"]["got"] and all(c is None for c in r["acted"]["column"])
        assert r["acted"]["summary"].startswith("Medians and CVs over 1 of 2 tissues")
