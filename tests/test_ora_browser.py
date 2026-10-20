"""M6's live over-representation in Chromium (GS.5), on the pipeline fixture's M6 page.

At the build's settings the figure stays as built; a direction, a background or a pasted list
reruns ctxOra in the page and redraws the passing sets. Skipped where no Chromium is found.
"""

import pytest
from test_export_browser import CHROMIUM, _browse

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
ACT = """(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(2000);
  const form = document.querySelector('form.ora'), fig = form.closest('.fig');
  const d = fig.querySelector('.js-plotly-plot');
  const cfg = JSON.parse(document.getElementById(form.getAttribute('data-controls')).textContent);
  const read = () => ({bars: d.data.filter(t => t.type === 'bar').map(t => [t.name, t.y.length]),
    hover: d.data.filter(t => t.type === 'bar').map(t => t.hovertemplate).join('|'),
    summary: fig.querySelector('.ora-summary').textContent,
    tsv: fig.querySelector('.ora-tsv').textContent.split('\\n').length,
    moved: !!fig.querySelector('a[data-export][role=button]')});
  const out = {start: read()};
  form.elements.direction.value = 'up';
  form.dispatchEvent(new Event('change'));
  await wait(400); out.up = read();
  form.elements.background.value = 'union';
  form.dispatchEvent(new Event('change'));
  await wait(400); out.union = read();
  form.elements.pasted.value = cfg.symbols.slice(0, 60).join(' ') + ' NOT_A_GENE';
  form.querySelector('[data-paste]').click();
  await wait(400); out.pasted = read();
  out.notes = fig.querySelector('.ora-notes').textContent;
  form.querySelector('[data-reset]').click();
  await wait(400); out.reset = read();
  return out;
})()"""


@pytest.fixture(scope="module")
def ran(pipeline, tmp_path_factory):
    root, _, _ = pipeline
    (page,) = (root / "outputs" / "06_function").glob("*_function_desynth_overview_v1.html")
    entry = root / "outputs" / "06_function" / "index.html"
    return entry, _browse([{"name": "ora", "url": page.as_uri(), "width": 1280, "height": 900,
                            "act": ACT}], tmp_path_factory.mktemp("downloads"))  # fmt: skip


def test_the_live_ora_reruns_on_each_setting(ran):
    """The build's bars stand until a setting moves; up only draws one family; the union
    background and a pasted list rerun; the table holds every kept set; the export follows.

    Mutations caught: the figure redrawn at load; the direction ignored; the paste counting
    a name with no symbol; the export switch left off."""
    entry, results = ran
    seen = results["ora"]
    assert seen["errors"] == []
    acted = seen["acted"]
    # Each direction's passing set and its grey fill to live.TOP (PK.3): four traces.
    assert not acted["start"]["moved"] and len(acted["start"]["bars"]) == 4  # noqa: PLR2004
    assert [n for n, _ in acted["up"]["bars"]] == ["up, FDR < 0.05 (BH)",
                                                   "up, not passing FDR < 0.05"]  # fmt: skip
    assert acted["up"]["moved"]
    live = "%{customdata[0]}<br>p %{customdata[4]"  # the live bars' hover; the build's its own
    assert live not in acted["start"]["hover"] and live in acted["up"]["hover"]
    assert "The featured contrast's calls" in acted["up"]["summary"]
    assert "gene-set union" in acted["union"]["summary"]
    assert acted["pasted"]["summary"].startswith("The pasted list: 60 genes")
    assert acted["notes"] == "60 pasted names are symbols on this page; names matching none: 1."
    assert acted["pasted"]["tsv"] > 1
    assert not acted["reset"]["moved"]
    assert entry.is_file()
