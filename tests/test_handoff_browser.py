"""M5's picker to M6 as a student opens the kit (PK): file://, a real mouse click on the onward
link or Open beside, the entry page, and what the reader sees on the page that opens.

Rytis's six genes reached M6 on 2026-10-20 and he saw none of them: the page opened at its
top, the live ORA panel 2,345 px down, and its notes counted the list without naming a gene.
Each test names the mutation it catches. Skipped where no Chromium is found; set
``CTX_CHROMIUM``.
"""

import re
from datetime import date

import pytest
from conftest import DE_SPEC, _ran, once
from synth import DE_GROUPS, design_counts, write_gather_outputs, write_gene_sets
from test_export_browser import CHROMIUM, _browse

from ctx.dataset import load
from ctx.report import genes

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
# Tick the picker's first `n` genes, or press a call button; `before` runs first.
PICK = """(async () => {{
  {before}
  await new Promise(r => setTimeout(r, 2500));
  const fig = document.querySelector('form.select').closest('.fig'), names = [];
  for (let i = 0; i < {n}; i++) {{
    // The table redraws on each tick, so each box is found again.
    const box = fig.querySelectorAll('table.pick tbody input[type=checkbox]')[i];
    names.push(box.getAttribute('aria-label')); box.click();
  }}
  {press}
  await new Promise(r => setTimeout(r, 600));
  return {{names, count: fig.querySelector('.select-count').textContent,
           href: fig.querySelector('a[data-onward]').getAttribute('href')}};
}})()"""
# A real click at the centre of the picker's first element `sel` names.
CLICK = """(() => {{
  const el = document.querySelector('form.select').closest('.fig').querySelector({sel!r});
  el.scrollIntoView({{block: 'center', behavior: 'instant'}});
  const r = el.getBoundingClientRect(), x = r.left + r.width / 2, y = r.top + r.height / 2;
  return [['mousePressed', x, y], ['mouseReleased', x, y]];
}})()"""
# What the reader of M6 sees once the page has loaded and settled.
SEEN = """(async () => {
  await new Promise(r => setTimeout(r, 2500));
  const fig = document.querySelector('form.ora').closest('.fig');
  const top = fig.getBoundingClientRect().top;
  return {url: location.href, hash: location.hash, opener: !!window.opener,
          notes: fig.querySelector('.ora-notes').textContent,
          summary: fig.querySelector('.ora-summary').textContent,
          top, window: innerHeight, scrolled: scrollY};
})()"""
# The window beside forwards through its entry page before `others` reads it.
WAIT = "new Promise(r => setTimeout(r, 3000)).then(() => null)"
# The address left unwritable, as a browser may leave a file:// page's.
REFUSED = "history.replaceState = () => { throw new DOMException('refused', 'SecurityError'); };"


def _task(name, m5, n=6, press="", before="", beside=False):
    """Pick on M5, then click the M6 link (or Open beside) with the mouse."""
    sel = "[data-beside]" if beside else "a[data-onward]"
    step = {"at": CLICK.format(sel=sel), "until": "/" if beside else "_function_desynth_overview",
            "then": WAIT if beside else SEEN}  # fmt: skip
    task = {"name": name, "url": m5.as_uri(), "width": 1280, "height": 900,
            "act": PICK.format(n=n, press=press, before=before), "mouse": [step]}  # fmt: skip
    return {**task, "others": SEEN} if beside else task


def _root(root):
    """de and function on the synthetic design, its genes under made-up Ensembl IDs: a made-up
    gene's g00001 is no Ensembl gene, and M5 then draws no picker."""
    write_gather_outputs(root, "desynth", design_counts(DE_GROUPS, batch_lfc=2.0), date.today(),
                         ensembl=True)  # fmt: skip
    write_gene_sets(root / load(DE_SPEC).resource("sets").dest)
    for module in ("de", "function"):
        _ran(root, module, "--dataset", "desynth", "--dataset-file", str(DE_SPEC))
    (m5,) = (root / "outputs" / "05_de").glob("*_de_desynth_overview_v1.html")
    return m5


@pytest.fixture(scope="module")
def ran(tmp_path_factory):
    """Each way of sending the list, in one browser."""
    m5 = once(tmp_path_factory, "handoff", _root)
    up = "fig.querySelector('[data-pick=up]').click();"
    tasks = [_task("link", m5), _task("beside", m5, beside=True), _task("many", m5, 0, up),
             _task("refused", m5, before=REFUSED)]  # fmt: skip
    return _browse(tasks, tmp_path_factory.mktemp("downloads"))


def _listed(notes: str) -> list[str]:
    """The names the notes give after 'The list: '."""
    return re.search(r"The list: (.*)\.$", notes)[1].split(", ")


def _in_view(seen: dict) -> bool:
    return 0 <= seen["top"] < seen["window"] / 2 and seen["scrolled"] > 0


def _where(seen: dict) -> tuple:
    return seen["top"], seen["window"], seen["scrolled"]


def test_a_clicked_link_shows_the_genes_on_m6(ran):
    """The click goes through the entry page to M6's overview with the six genes; M6 names
    each and brings its live ORA into the top half of the window.

    Mutations caught: the reveal left out (the panel 2,000 px down, as Rytis saw it); the
    names left out of the notes; the entry page dropping the fragment."""
    seen = ran["link"]
    assert seen["errors"] == []
    picked, got = seen["acted"], seen["moused"][0]
    assert picked["count"].startswith("6 genes selected")
    assert got["url"].split("#")[0].endswith("_function_desynth_overview_v1.html")
    assert genes.decode(got["hash"]).genes and "&n=6&" in got["hash"]
    assert got["notes"].startswith("6 genes of the list have a symbol on this page")
    assert sorted(_listed(got["notes"])) == sorted(picked["names"])
    assert got["summary"].startswith("The list of the link: 6 genes")
    assert _in_view(got), _where(got)


def test_open_beside_shows_the_genes_in_the_window_it_opens(ran):
    """Open beside, clicked: the named window opens M6 with the list, names the genes and shows
    the panel. Mutation caught: the reveal left out."""
    seen = ran["beside"]
    assert seen["errors"] == []
    (got,) = (v for url, v in seen["others"].items() if "_function_desynth_overview" in url)
    assert got["opener"]
    assert sorted(_listed(got["notes"])) == sorted(seen["acted"]["names"])
    assert _in_view(got), _where(got)


def test_a_long_list_names_its_first_genes_and_counts_the_rest(ran):
    """Every up call sent: the notes name ``genes.PAGE_ROWS`` genes, then the rest by count.
    Mutation caught: the cap ignored, a list of thousands of names in one line."""
    got = ran["many"]["moused"][0]
    n = int(re.match(r"(\d+) genes of the list", got["notes"])[1])
    assert n > genes.PAGE_ROWS, n
    shown = re.search(r"The list: (.*) and (\d+) more\.$", got["notes"])
    assert len(shown[1].split(", ")) == genes.PAGE_ROWS
    assert int(shown[2]) == n - genes.PAGE_ROWS


def test_a_refused_address_still_sends_the_list(ran):
    """A page whose address cannot be rewritten (a browser's file:// rule) still links the
    list and M6 receives it. Mutation caught: replaceState left unguarded, which threw
    before the onward link and the windows beside were updated."""
    seen = ran["refused"]
    assert seen["errors"] == []
    got = seen["moused"][0]
    assert "&n=6&" in seen["acted"]["href"] and "&n=6&" in got["hash"]
    assert sorted(_listed(got["notes"])) == sorted(seen["acted"]["names"])
