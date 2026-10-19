"""M6's running score for any set, the background view and the live ORA's contrast and set
count (WO7), in Chromium on the pipeline fixture's M6 page, at 1,280 px light and 375 px dark.

Each figure is read back from plotly and checked against numbers computed here: the running
sum by ``lite``, the ORA rows by the tables ``ctx function`` wrote, every p within D-049's
1.52e-10. Each test names the mutations it catches.
"""

import math
import re

import numpy as np
import polars as pl
import pytest
from test_export_browser import CHROMIUM, _browse
from views_browser import walk

from ctx.function import lite, live, running
from ctx.function.page import set_words

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
SIZES = (("wide_light", 1280, "light"), ("phone_dark", 375, "dark"))
RUNNING = """(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(1500);
  const form = document.querySelector('form.running'), fig = form.closest('.fig');
  const d = fig.querySelector('.js-plotly-plot');
  const cfg = JSON.parse(document.getElementById(form.getAttribute('data-controls')).textContent);
  const words = s => (fig.querySelector(s) || {}).textContent || '';
  const read = () => ({line: Array.from(d.data.find(t => t.name === cfg.line).y),
    ranks: Array.from(d.data.find(t => t.name === cfg.line).x),
    dots: [Array.from(d.data.find(t => t.name === cfg.member).x),
           Array.from(d.data.find(t => t.name === cfg.member).y)],
    title: words('.fig-title'), lead: words('.fig-lead'), caption: words('figcaption'),
    line_words: words('.running-summary'), es: fig.ctxRunning ? fig.ctxRunning.es : null,
    moved: !!fig.querySelector('a[data-export][role=button]'),
    export: fig.ctxRunning ? fig.ctxRunning.caption : cfg.export.caption,
    label: cfg.sets[Number(form.elements.set.value)].label});
  const out = {loaded: read(), sets: []};
  for (let i = 0; i < cfg.sets.length; i++) {
    form.elements.set.value = String(i); form.dispatchEvent(new Event('change'));
    await wait(80); out.sets.push(read());
  }
  form.elements.set.value = '0'; form.dispatchEvent(new Event('change'));
  await wait(80); out.back = read();
  out.box = document.documentElement.clientWidth; out.scroll = document.documentElement.scrollWidth;
  return out;
})()"""
ORA = """(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(1500);
  const form = document.querySelector('form.ora'), fig = form.closest('.fig');
  const d = fig.querySelector('.js-plotly-plot');
  const read = () => ({tsv: fig.querySelector('.ora-tsv').textContent,
    bars: d.data.filter(t => t.type === 'bar').map(t => [t.name, Array.from(t.y)]),
    summary: fig.querySelector('.ora-summary').textContent,
    background: form.elements.background.value,
    moved: !!fig.querySelector('a[data-export][role=button]')});
  const out = {contrasts: []};
  for (const option of Array.from(form.elements.contrast.options)) {
    form.elements.contrast.value = option.value;
    form.elements.contrast.dispatchEvent(new Event('change', {bubbles: true}));
    await wait(150); out.contrasts.push([option.value, read()]);
  }
  form.elements.contrast.value = form.elements.contrast.options[0].value;
  form.elements.contrast.dispatchEvent(new Event('change', {bubbles: true}));
  form.elements.sets.value = 'top';
  form.elements.sets.dispatchEvent(new Event('change', {bubbles: true}));
  await wait(150); out.top = read();
  form.querySelector('[data-reset]').click();
  await wait(150); out.reset = read();
  return out;
})()"""


CONTRASTS = """(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(1500);
  const form = document.querySelector('form.running'), fig = form.closest('.fig');
  const d = fig.querySelector('.js-plotly-plot');
  const cfg = JSON.parse(document.getElementById(form.getAttribute('data-controls')).textContent);
  const words = s => (fig.querySelector(s) || {}).textContent || '';
  const read = () => ({contrast: form.elements.contrast.selectedOptions[0].textContent,
    sets: Array.from(form.elements.set.options).map(o => o.textContent),
    set: form.elements.set.selectedOptions[0].textContent,
    line: Array.from(d.data.find(t => t.name === cfg.line).y),
    dots: Array.from(d.data.find(t => t.name === cfg.member).x),
    es: fig.ctxRunning ? fig.ctxRunning.es : null, caption: words('figcaption'),
    line_words: words('.running-summary'),
    moved: !!fig.querySelector('a[data-export][role=button]'),
    export: fig.ctxRunning ? fig.ctxRunning.caption : cfg.export.caption});
  const out = {loaded: read(), contrasts: []};
  /* Each contrast in turn, then the build's second set kept across a contrast that tests it. */
  for (const option of Array.from(form.elements.contrast.options)) {
    form.elements.contrast.value = option.value;
    form.elements.contrast.dispatchEvent(new Event('change', {bubbles: true}));
    await wait(80); out.contrasts.push(read());
  }
  form.elements.contrast.value = '0';
  form.elements.contrast.dispatchEvent(new Event('change', {bubbles: true}));
  await wait(80); out.back = read();
  out.box = document.documentElement.clientWidth; out.scroll = document.documentElement.scrollWidth;
  return out;
})()"""


@pytest.fixture(scope="module")
def m6(pipeline, tmp_path_factory):
    """(the M6 and M5 folders, the page, and what each act read at both sizes)."""
    root = pipeline[0]
    folder = root / "outputs" / "06_function"
    (page,) = folder.glob("*_function_desynth_overview_v1.html")
    tasks = [{"name": f"{what}_{name}", "url": page.as_uri(), "width": width, "height": 900,
              "scheme": scheme, "act": act}
             for what, act in (("running", RUNNING), ("ora", ORA), ("contrasts", CONTRASTS))
             for name, width, scheme in SIZES]  # fmt: skip
    seen = _browse(tasks, tmp_path_factory.mktemp("downloads"))
    views = walk(page, tmp_path_factory.mktemp("walk"), SIZES)
    return folder, root / "outputs" / "05_de", page, seen, views


def _newest(folder, what: str) -> pl.DataFrame:
    (path,) = folder.glob(f"*_desynth_{what}_v1.parquet")
    return pl.read_parquet(path)


def _sets(folder) -> dict[str, tuple[str, ...]]:
    """The fixture's gene sets, as ``ctx function`` read them."""
    (gmt,) = (folder.parents[1] / "resources").rglob("*.gmt")
    return lite.read_gmt(gmt)


@pytest.mark.parametrize("size", ["wide_light", "phone_dark"])
def test_every_set_draws_lites_running_sum_and_its_own_words(m6, size):
    """Mutations caught: the twin's misses counting the members; the dots left on the build's
    set; a set's caption not swapped; the export kept on the build's caption; the ES check
    never run."""
    folder, de, _, seen, _ = m6
    result = seen[f"running_{size}"]
    assert result["errors"] == []
    acted = result["acted"]
    assert acted["scroll"] <= acted["box"]
    gsea = _newest(folder, "a_case_gsea")
    weights = running.weights(_newest(de, "a_case_results"))
    names, scores = running.ranking(_newest(de, "a_case_results"))
    ranked, _ = lite.ranking(names, scores)
    rank_of = {g: i for i, g in enumerate(ranked)}
    members = {set_words(s): genes for s, genes in _sets(folder).items()}
    by_words = {set_words(r["gene_set"]): r for r in gsea.iter_rows(named=True)}
    assert len(acted["sets"]) == gsea.height
    for i, drawn in enumerate(acted["sets"]):
        row = by_words[drawn["label"]]
        positions = np.sort([rank_of[g] for g in members[drawn["label"]] if g in rank_of])
        assert drawn["dots"][0] == positions.tolist()
        expected = lite._running(positions, weights[positions], len(weights))  # noqa: SLF001
        assert np.allclose(drawn["line"], expected, rtol=0, atol=1e-12)
        assert drawn["ranks"] == list(range(len(weights)))
        assert np.allclose(drawn["dots"][1], np.asarray(expected)[positions], atol=1e-12)
        assert math.isclose(drawn["es"], row["es"], rel_tol=0, abs_tol=1e-9)
        assert f"ES = {row['es']:.2f}, NES = {row['nes']:.2f}" in drawn["caption"]
        assert drawn["title"].startswith(f"Running score: {drawn['label']} in ")
        assert "gives ES" not in drawn["line_words"]
        assert drawn["moved"] == (i != 0)
        assert drawn["export"].startswith(f"{drawn['label']} in ")
        assert drawn["caption"].startswith(f"{drawn['label']} in ")
    assert acted["back"]["line"] == acted["sets"][0]["line"] and not acted["back"]["moved"]
    assert acted["loaded"]["caption"] == acted["sets"][0]["caption"]


@pytest.mark.parametrize("size", ["wide_light", "phone_dark"])
def test_every_contrast_draws_its_own_ranking_and_lists_its_own_sets(m6, size):
    """The running score's contrast select (V2): each contrast's sets in ``running.order``,
    the set kept by ``running.keep``, the line ``lite._running`` on that contrast's ranking.

    Mutations caught: the build's ranking kept under another contrast; the build's set list
    kept; the shown set reset to the first where the next contrast tested it; the export left
    on the build's caption; the contrast select's return leaving another contrast's line."""
    from conftest import DE_SPEC  # noqa: PLC0415

    from ctx.dataset import load  # noqa: PLC0415
    from ctx.de.page import every  # noqa: PLC0415

    folder, de, _, seen, _ = m6
    result = seen[f"contrasts_{size}"]
    assert result["errors"] == []
    acted = result["acted"]
    assert acted["scroll"] <= acted["box"]
    spec = load(DE_SPEC)
    by_label = {w.label: n for n, w in every(spec.contrasts, spec.columns).items()}
    members = _sets(folder)
    shown = "UP_SET"
    assert by_label[acted["contrasts"][0]["contrast"]] == "a_case"
    assert len(acted["contrasts"]) > 1
    for i, drawn in enumerate(acted["contrasts"]):
        name = by_label[drawn["contrast"]]
        gsea = _newest(folder, f"{name}_gsea")
        tested = {r["gene_set"]: r for r in gsea.iter_rows(named=True)}
        order = running.order(tested, "UP_SET")
        assert drawn["sets"] == [set_words(s) for s in order]
        kept = order[running.keep(order, shown, "UP_SET")]
        assert drawn["set"] == set_words(kept)
        shown = kept
        results = _newest(de, f"{name}_results")
        weights = running.weights(results)
        ranked, _ = lite.ranking(*running.ranking(results))
        rank_of = {g: r for r, g in enumerate(ranked)}
        positions = np.sort([rank_of[g] for g in members[kept] if g in rank_of])
        assert drawn["dots"] == positions.tolist()
        expected = lite._running(positions, weights[positions], len(weights))  # noqa: SLF001
        assert np.allclose(drawn["line"], expected, rtol=0, atol=1e-12)
        assert math.isclose(drawn["es"], tested[kept]["es"], rel_tol=0, abs_tol=1e-9)
        assert drawn["moved"] == (i != 0) and drawn["export"].startswith(set_words(kept))
        assert drawn["caption"].startswith(f"{set_words(kept)} in {drawn['contrast']}")
    assert acted["back"]["line"] == acted["contrasts"][0]["line"] and not acted["back"]["moved"]


@pytest.mark.parametrize("size", ["wide_light", "phone_dark"])
def test_another_contrasts_calls_rerun_as_its_own_table(m6, size):
    """Mutations caught: the build's calls kept for another contrast; the background left on
    the build's tested genes; the set count ignored; reset leaving another contrast."""
    folder, _, _, seen, _ = m6
    result = seen[f"ora_{size}"]
    assert result["errors"] == []
    acted = result["acted"]
    for contrast, drawn in acted["contrasts"]:
        assert drawn["background"] == contrast
        table = _newest(folder, f"{contrast}_ora")
        rows = [line.split("\t") for line in drawn["tsv"].splitlines()[1:]]
        assert len(rows) == table.height, contrast
        expected = {(r["direction"], r["gene_set"]): r for r in table.iter_rows(named=True)}
        for direction, gene_set, size_, hits, _, p, fdr, _ in rows:
            want = expected[(direction, gene_set)]
            assert (int(size_), int(hits)) == (want["size"], want["hits"])
            assert math.isclose(float(p), want["p"], rel_tol=1.52e-10, abs_tol=0)
            assert math.isclose(float(fdr), want["fdr"], rel_tol=1.52e-10, abs_tol=0)
    assert acted["contrasts"][1][1]["moved"] and not acted["contrasts"][0][1]["moved"]
    top = acted["top"]
    # The plan's seven smallest FDR per family, whether they pass or not.
    assert top["moved"] and [len(y) for _, y in top["bars"]] == [7, 7]
    assert not acted["reset"]["moved"]
    assert acted["reset"]["tsv"] == acted["contrasts"][0][1]["tsv"]


@pytest.mark.parametrize("size", ["wide_light", "phone_dark"])
def test_the_background_views_draw_lites_p_for_each_setting(m6, size):
    """Mutations caught: the sets drawn capped at seven in the every-set view; the down calls
    tested as the up; a view's wider background drawn from the build's."""
    folder, de, _, _, views = m6
    acted = views[size]["acted"]
    assert views[size]["errors"] == []
    (view,) = [v for v in acted["views"] if v["label"] == "Background"]
    results = _newest(de, "a_case_results")
    sets = _sets(folder)
    tested = live.symbols(results)
    union = sorted({g for gs in sets.values() for g in gs} | set(tested))
    seen = set()
    for option in view["options"]:
        direction = "down" if "down calls" in (option["line"] or "") else "up"
        every_set = "every set tested" in (option["line"] or "")
        seen.add((direction, every_set))
        hits = live.symbols(results, direction)
        right = {r.name: r for r in lite.ora(hits, tested, sets)}
        wrong = {r.name: r for r in lite.ora(hits, union, sets)}
        order = sorted((s for s in right if s in wrong), key=lambda s: right[s].p)
        order = order if every_set else order[:7]
        drawn = {t["name"]: (list(t["y"]), list(t["x"])) for t in option["traces"]}
        mine, wide = drawn["tested genes"], drawn["gene-set union"]
        assert mine[0] == wide[0] == [set_words(s) for s in order]
        assert np.allclose(mine[1], [-math.log10(right[s].p) for s in order])
        assert np.allclose(wide[1], [-math.log10(wrong[s].p) for s in order])
        assert f" {direction} calls" in option["caption"]
    assert seen == {("up", False), ("up", True), ("down", False), ("down", True)}


def test_another_set_exports_its_own_curve_and_caption(m6, tmp_path):
    """Mutations caught: the export switch left off for another set; the build's caption
    printed under another set's curve."""
    _, _, page, _, _ = m6
    act = """(async () => {
      await new Promise(r => setTimeout(r, 1500));
      const form = document.querySelector('form.running');
      form.elements.set.value = '1'; form.dispatchEvent(new Event('change'));
      await new Promise(r => setTimeout(r, 200));
      return [form.closest('.fig').id, form.closest('.fig').ctxRunning.caption];
    })()"""
    probe = _browse([{"name": "probe", "url": page.as_uri(), "act": act}], tmp_path)["probe"]
    figure, caption = probe["acted"]
    (result,) = _browse([{"name": "export", "url": page.as_uri(), "act": act,
                          "figure": f"#{figure}"}], tmp_path).values()  # fmt: skip
    svg = (tmp_path / result["files"]["svg"]["guid"]).read_text()
    words = " ".join(re.findall(r">([^<>]+)</text>", svg))
    assert caption.split(":")[0] in words and "This view shows" in words
    assert (tmp_path / result["files"]["pdf"]["guid"]).read_bytes().startswith(b"%PDF")
