"""M10.G in Chromium: a gene typed on Figure 1 redraws every ER figure as that gene's patients at
or above its median against below, each on its own endpoint, equal to ctx.clinical.gene.split
on the page's own data; the screen's Cox numbers for the gene beside it; the view's PDF and SVG;
the ER groups back on Reset; a list from M5's link. At 1,280 and 375 px, light and dark, with no
console error and no sideways scroll. Each test names the mutation it catches.
"""

import json

import numpy as np
import pytest
from test_clinical_page import clinical_root, controls, outputs  # noqa: F401
from test_export_browser import CHROMIUM, _browse
from test_gene_survival import page_patients

from ctx.clinical import gene, panel
from ctx.dataset import ENDPOINTS
from ctx.report import genes, text
from ctx.report.builder import unpacked

# Every ER figure's plot after a gene is typed on the first, then the first after Reset.
ACT = """(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(3000);
  const forms = Array.from(document.querySelectorAll('form.gene-km'));
  const read = form => {
    const fig = form.closest('.fig'), div = fig.querySelector('.fig-plot .js-plotly-plot');
    const own = document.getElementById(form.getAttribute('data-controls'));
    return {endpoint: JSON.parse(own.textContent).endpoint,
      traces: div.data.map(t => ({x: Array.from(t.x || []), y: Array.from(t.y || []), name: t.name,
                                  visible: t.visible, text: t.text ? Array.from(t.text) : null})),
      ticks: Array.from(div.layout.xaxis2.tickvals || []),
      summary: fig.querySelector('.gene-summary').textContent,
      notes: fig.querySelector('.gene-notes').textContent,
      caption: fig.querySelector('figcaption').textContent,
      title: fig.querySelector('.fig-title').textContent,
      svg: fig.querySelector('a[data-export=svg]').textContent};
  };
  const out = {built: read(forms[0])};
  forms[0].elements.gene.value = 'NOTAGENE';
  forms[0].querySelector('[type=submit]').click();
  await wait(500); out.missing = read(forms[0]);
  forms[0].elements.gene.value = 'SYMBOL';
  forms[0].querySelector('[type=submit]').click();
  await wait(1500);
  out.drawn = forms.map(read);
  out.wide = document.documentElement.scrollWidth, out.window = innerWidth;
  out.theme = document.documentElement.getAttribute('data-theme');
  return out;
})()"""
RESET = """(async () => {
  const form = document.querySelector('form.gene-km'), fig = form.closest('.fig');
  form.querySelector('[data-reset]').click();
  await new Promise(r => setTimeout(r, 1000));
  const div = fig.querySelector('.fig-plot .js-plotly-plot');
  return {names: div.data.map(t => t.name), caption: fig.querySelector('figcaption').textContent,
          title: fig.querySelector('.fig-title').textContent,
          svg: fig.querySelector('a[data-export=svg]').textContent};
})()"""
LINKED = """(async () => {
  await new Promise(r => setTimeout(r, 3500));
  const form = document.querySelector('form.gene-km'), fig = form.closest('.fig');
  const menu = form.elements.listed;
  return {gene: form.elements.gene.value, menu: Array.from(menu.options).map(o => o.textContent),
          shown: !menu.closest('label').hidden, notes: fig.querySelector('.gene-notes').textContent,
          names: fig.querySelector('.fig-plot .js-plotly-plot').data.map(t => t.name)};
})()"""
FIGURE = ".fig:has(form.gene-km)"
SIZES = {"wide_light": (1280, 900, "light"), "wide_dark": (1280, 900, "dark"),
         "phone_light": (375, 812, "light"), "phone_dark": (375, 812, "dark")}  # fmt: skip


def chosen(cfg: dict) -> int:
    """The first panel gene whose split on every endpoint has two sides."""
    values = panel.decode(cfg["block"], cfg["genes"], cfg["patients"])
    patients = page_patients(cfg)
    for row in range(cfg["genes"]):
        if all(gene.split(patients, e, values[row]).logrank for e in ENDPOINTS):
            return row
    raise AssertionError("no gene splits in two on every endpoint")


@pytest.fixture(scope="module")
def ran(clinical_root, tmp_path_factory):  # noqa: F811 - the fixture
    page = outputs(clinical_root, "overview_v1.html")
    cfg = controls(unpacked(page.read_text()), "survival")
    row = chosen(cfg)
    act = ACT.replace("SYMBOL", cfg["symbols"][row])
    tasks = [{"name": name, "url": page.as_uri(), "width": w, "height": h, "scheme": scheme,
              "act": act, "figure": FIGURE, "reset": RESET}
             for name, (w, h, scheme) in SIZES.items()]  # fmt: skip
    listed = sorted(cfg["numbers"][row : row + 3])
    link = genes.encode(genes.Selection("fake_brca:a_case", tuple(listed)))
    tasks.append({"name": "linked", "url": page.as_uri() + link, "act": LINKED})
    downloads = tmp_path_factory.mktemp("downloads")
    return cfg, row, _browse(tasks, downloads), downloads, listed


def expected(cfg: dict, row: int, endpoint: str) -> gene.GeneSplit:
    """Python's split of the gene on the page's own data."""
    values = panel.decode(cfg["block"], cfg["genes"], cfg["patients"])
    return gene.split(page_patients(cfg), endpoint, values[row], gene.group_names())


@pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
@pytest.mark.parametrize("name", list(SIZES))
def test_a_typed_gene_redraws_every_er_figure_as_python_splits_it(ran, name):
    """Curves, censor ticks, at-risk rows and log-rank p equal Python's on all four ER
    figures. Mutations caught: the twin's split at > rather than >=; its log-rank variance
    without the ties factor; an at-risk count from the time before a tick; the box drawing the
    first figure alone; another gene's plane read."""
    cfg, row, results, _, _ = ran
    seen = results[name]
    assert seen["errors"] == []
    a = seen["acted"]
    assert [d["endpoint"] for d in a["drawn"]] == list(ENDPOINTS)
    below, above = gene.group_names()
    for drawn in a["drawn"]:
        found = expected(cfg, row, drawn["endpoint"])
        t = drawn["traces"]
        assert [x["name"] for x in t] == [below] * 3 + [above] * 3 + [below, above]
        for slot, side in enumerate((below, above)):
            steps = found.figure.steps(side)
            np.testing.assert_allclose(t[3 * slot + 1]["x"], steps["x"], rtol=1e-12)
            np.testing.assert_allclose(t[3 * slot + 1]["y"], steps["y"], rtol=1e-12)
            rows = found.figure.frame.filter(found.figure.frame["group"] == side)
            censored = rows.filter(rows["censored"] > 0)
            np.testing.assert_allclose(t[3 * slot + 2]["x"], censored["time"].to_numpy(),
                                       rtol=1e-12)  # fmt: skip
            assert t[6 + slot]["text"] == [f"{n:,}" for n in found.figure.at_risk(side)]
        assert drawn["ticks"] == list(found.figure.risk_times)
        chi, _, p = found.logrank
        assert f"chi-squared {chi:.2f}, p = {format(p, '.2g')}" in drawn["summary"]
        assert f"chi-squared {chi:.2f}, p = {format(p, '.2g')}" in drawn["caption"]
        assert cfg["symbols"][row] in drawn["title"]


@pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
@pytest.mark.parametrize("name", list(SIZES))
def test_the_screens_numbers_stand_beside_the_curves(ran, name):
    """The screen's HR per SD, p and q of the gene, alone and adjusted, as its table prints
    them, and the endpoint they were fitted on where the figure draws another. Mutations
    caught: another row's cells; the p values left out of screen_rows."""
    cfg, row, results, _, _ = ran
    cells, ps = cfg["screen"][row][2], cfg["screen"][row][4]
    for drawn in results[name]["acted"]["drawn"]:
        line = drawn["summary"]
        assert (f"HR {cells[1]} alone, p = {ps[0]}, q = {cells[2]}; HR {cells[3]} ({cells[4]}) "
                f"adjusted, p = {ps[1]}, q = {cells[5]}") in line  # fmt: skip
        other = drawn["endpoint"] != cfg["endpoint"]
        assert (f"The screen fits {cfg['endpoint']} and this figure draws" in line) == other


@pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
@pytest.mark.parametrize("name", list(SIZES))
def test_no_sideways_scroll_and_an_unknown_gene_is_named(ran, name):
    """Mutations caught: the box wider than a phone; an unknown gene drawing nothing and
    saying nothing."""
    cfg, _, results, _, _ = ran
    a = results[name]["acted"]
    assert a["wide"] <= a["window"], (a["wide"], a["window"])
    assert a["theme"] == ("black" if name.endswith("dark") else "white")
    assert a["missing"]["notes"] == text.fill("clinical", "gene.missing", gene="NOTAGENE",
                                              genes=cfg["genes"])  # fmt: skip
    assert [t["name"] for t in a["missing"]["traces"]] == [t["name"] for t in a["built"]["traces"]]


@pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
@pytest.mark.parametrize("name", list(SIZES))
def test_the_pdf_and_svg_draw_the_gene_and_reset_restores_the_er_groups(ran, name):
    """The buttons save the view as set: the gene's groups and caption in the SVG, a PDF.
    Mutations caught: the view's caption left at the ER figure's; Reset leaving the gene."""
    cfg, row, results, downloads, _ = ran
    seen = results[name]
    files = seen["files"]
    svg = (downloads / files["svg"]["guid"]).read_text()
    below, above = gene.group_names()
    assert below in svg and above in svg
    assert f"split at the median of {cfg['symbols'][row]}" in svg
    assert (downloads / files["pdf"]["guid"]).read_bytes().startswith(b"%PDF")
    assert files["svg"]["name"].endswith("_view.svg")
    a, reset = seen["acted"], seen["reset"]
    assert a["drawn"][0]["svg"] != a["built"]["svg"]
    assert reset["names"] == [t["name"] for t in a["built"]["traces"]]
    assert reset["caption"] == a["built"]["caption"] and reset["title"] == a["built"]["title"]
    assert reset["svg"] == a["built"]["svg"]


@pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
def test_a_list_from_a_link_fills_the_menu_and_draws_its_first_gene(ran):
    """D-049's fragment: the list's genes in the box's menu, the first drawn. Mutation caught:
    the box ignoring the fragment."""
    cfg, _, results, _, listed = ran
    seen = results["linked"]
    assert seen["errors"] == []
    a = seen["acted"]
    rows = [cfg["numbers"].index(n) for n in listed]
    assert a["menu"] == [cfg["symbols"][r] for r in rows] and a["shown"]
    assert a["gene"] == cfg["symbols"][rows[0]]
    assert a["names"][:1] == [gene.group_names()[0]]
    assert json.dumps(len(rows)) in a["notes"]
