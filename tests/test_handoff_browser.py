"""M5's picker to M6 as a student opens the kit (PK): file://, a real mouse click on the onward
link or Open beside, the entry page, and what the reader sees on the page that opens.

Rytis's six genes reached M6 on 2026-10-20 and he saw none of them: the page opened at its
top, the live ORA panel 2,345 px down, and its notes counted the list without naming a gene.
Then the panel drew an empty chart, since no set passed FDR (PK.3): it now draws the sets of
smallest p in grey, checked here against ``live.bars`` at 1,280 and 375 px in both themes.
Each test names the mutation it catches. Skipped where no Chromium is found; set
``CTX_CHROMIUM``.
"""

import html
import math
import re
from datetime import date
from itertools import pairwise

import polars as pl
import pytest
from conftest import DE_SPEC, _ran, once
from synth import DE_GROUPS, design_counts, write_gather_outputs, write_gene_sets
from test_export_browser import CHROMIUM, _browse

from ctx.dataset import load
from ctx.function import lite, live
from ctx.function.page import set_words
from ctx.report import genes, text

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
# What the reader of M6 sees once the page has loaded and settled: the live ORA's bars too.
SEEN = """(async () => {
  await new Promise(r => setTimeout(r, 2500));
  const form = document.querySelector('form.ora'), fig = form.closest('.fig');
  const top = fig.getBoundingClientRect().top, d = fig.querySelector('.js-plotly-plot');
  const cfg = JSON.parse(document.getElementById(form.getAttribute('data-controls')).textContent);
  return {url: location.href, hash: location.hash, opener: !!window.opener,
          notes: fig.querySelector('.ora-notes').textContent,
          summary: fig.querySelector('.ora-summary').textContent,
          bars: d.data.filter(t => t.type === 'bar').map(t => ({name: t.name,
            colour: t.marker.color, y: Array.from(t.y), x: Array.from(t.x)})),
          categories: Array.from(d.layout.yaxis.categoryarray || []),
          shapes: (d.layout.shapes || []).map(s => ({x0: s.x0, x1: s.x1, xref: s.xref,
            dash: s.line && s.line.dash, name: s.name})),
          legend: d.layout.showlegend, grey: cfg.colours.failing,
          moved: !!fig.querySelector('a[data-export][role=button]'),
          box: document.documentElement.clientWidth, scroll: document.documentElement.scrollWidth,
          ticks: Array.from(d.querySelectorAll('.yaxislayer-above .ytick text')).map(t => ({
            text: t.textContent, lines: Math.max(1, t.querySelectorAll('tspan.line').length),
            left: t.getBoundingClientRect().left, right: t.getBoundingClientRect().right,
            top: t.getBoundingClientRect().top, bottom: t.getBoundingClientRect().bottom})),
          height: d._fullLayout.height, plot: d._fullLayout._size.h,
          width: d.getBoundingClientRect().width, top, window: innerHeight, scrolled: scrollY};
})()"""
# The window beside forwards through its entry page before `others` reads it.
WAIT = "new Promise(r => setTimeout(r, 3000)).then(() => null)"
# The address left unwritable, as a browser may leave a file:// page's.
REFUSED = "history.replaceState = () => { throw new DOMException('refused', 'SecurityError'); };"
# The live ORA's figure on M6, whose PDF and SVG buttons the driver clicks.
ORA_FIG = ".fig:has(form.ora)"
# Rytis's 6 genes at both widths and both themes (PK.3): name to (width, height, scheme).
SIZES = {
    "link": (1280, 900, "light"),
    "link_wide_dark": (1280, 900, "dark"),
    "link_phone_light": (375, 812, "light"),
    "link_phone_dark": (375, 812, "dark"),
}
# M6 opened as built, no list (PK.4): name to (width, height, scheme).
BUILT = {
    "built_wide": (1280, 900, "light"),
    "built_phone_light": (375, 812, "light"),
    "built_phone_dark": (375, 812, "dark"),
}
# Step-1's labels are MSigDB Hallmark names, up to 51 characters in words: the fixture's other
# sets take real Hallmark names, so its labels run as wide (PK.4). UP_SET and DOWN_SET stay,
# the names its dataset file tests.
HALLMARK = (
    "INTERFERON_GAMMA_RESPONSE", "EPITHELIAL_MESENCHYMAL_TRANSITION", "INFLAMMATORY_RESPONSE",
    "IL6_JAK_STAT3_SIGNALING", "IL2_STAT5_SIGNALING", "OXIDATIVE_PHOSPHORYLATION",
    "UNFOLDED_PROTEIN_RESPONSE", "REACTIVE_OXYGEN_SPECIES_PATHWAY", "WNT_BETA_CATENIN_SIGNALING",
    "TGF_BETA_SIGNALING", "PI3K_AKT_MTOR_SIGNALING", "MTORC1_SIGNALING",
    "ESTROGEN_RESPONSE_EARLY", "ESTROGEN_RESPONSE_LATE", "ANDROGEN_RESPONSE",
    "CHOLESTEROL_HOMEOSTASIS", "FATTY_ACID_METABOLISM", "XENOBIOTIC_METABOLISM",
    "KRAS_SIGNALING_UP", "KRAS_SIGNALING_DN", "TNFA_SIGNALING_VIA_NFKB", "UV_RESPONSE_DN",
    "HEDGEHOG_SIGNALING", "NOTCH_SIGNALING", "APICAL_JUNCTION", "APICAL_SURFACE",
    "ADIPOGENESIS", "ALLOGRAFT_REJECTION", "COAGULATION", "COMPLEMENT", "HYPOXIA",
    "G2M_CHECKPOINT",
)  # fmt: skip
GREY_SVG = re.compile(r"fill:\s*rgb\(140,\s*140,\s*140\)|fill=\"#8c8c8c\"", re.IGNORECASE)


def _task(name, m5, n=6, press="", before="", beside=False):  # noqa: PLR0913 - a task's parts
    """Pick on M5, then click the M6 link (or Open beside) with the mouse."""
    sel = "[data-beside]" if beside else "a[data-onward]"
    step = {"at": CLICK.format(sel=sel), "until": "/" if beside else "_function_desynth_overview",
            "then": WAIT if beside else SEEN}  # fmt: skip
    width, height, scheme = SIZES.get(name, SIZES["link"])
    task = {"name": name, "url": m5.as_uri(), "width": width, "height": height,
            "scheme": scheme, "act": PICK.format(n=n, press=press, before=before),
            "mouse": [step]}  # fmt: skip
    if name in SIZES:
        task["figure"] = ORA_FIG
    return {**task, "others": SEEN} if beside else task


def _root(root):
    """de and function on the synthetic design, its genes under made-up Ensembl IDs: a made-up
    gene's g00001 is no Ensembl gene, and M5 then draws no picker. Twelve more sets each hold
    one planted gene among others, so a short list's sets of smallest p hold a gene of it and
    fail FDR, as Rytis's 6 genes did on Step-1 (PK.3)."""
    write_gather_outputs(root, "desynth", design_counts(DE_GROUPS, batch_lfc=2.0), date.today(),
                         ensembl=True)  # fmt: skip
    gmt = root / load(DE_SPEC).resource("sets").dest
    write_gene_sets(gmt)
    filler = [f"G{i:05d}" for i in range(400, 700)]
    one = [f"ONE_{i:02d}\tsynthetic\t" + "\t".join([f"G{i + 1:05d}", *filler[i * 20:i * 20 + 30]])
           for i in range(12)]  # fmt: skip
    named = iter(HALLMARK)
    sets = re.sub(r"^(?:RANDOM|ONE)_\d\d\t", lambda _: f"HALLMARK_{next(named)}\t",
                  gmt.read_text() + "\n".join(one) + "\n", flags=re.MULTILINE)  # fmt: skip
    assert next(named, None) is None, "one Hallmark name for each other set"
    gmt.write_text(sets)
    for module in ("de", "function"):
        _ran(root, module, "--dataset", "desynth", "--dataset-file", str(DE_SPEC))
    (m5,) = (root / "outputs" / "05_de").glob("*_de_desynth_overview_v1.html")
    return m5


@pytest.fixture(scope="module")
def ran(tmp_path_factory):
    """Each way of sending the list, in one browser: (M5 page, downloads, results)."""
    m5 = once(tmp_path_factory, "handoff", _root)
    up = "fig.querySelector('[data-pick=up]').click();"
    (m6,) = (m5.parents[1] / "06_function").glob("*_function_desynth_overview_v1.html")
    built = [{"name": n, "url": m6.as_uri(), "width": w, "height": h, "scheme": c, "act": SEEN}
             for n, (w, h, c) in BUILT.items()]  # fmt: skip
    tasks = [*(_task(name, m5) for name in SIZES), _task("beside", m5, beside=True),
             _task("many", m5, 0, up), _task("refused", m5, before=REFUSED), *built]  # fmt: skip
    downloads = tmp_path_factory.mktemp("downloads")
    return m5, downloads, _browse(tasks, downloads)


def _listed(notes: str) -> list[str]:
    """The names the notes give after 'The list: '."""
    return re.search(r"The list: (.*)\.$", notes)[1].split(", ")


def _in_view(seen: dict) -> bool:
    return 0 <= seen["top"] < seen["window"] / 2 and seen["scrolled"] > 0


def _where(seen: dict) -> tuple:
    return seen["top"], seen["window"], seen["scrolled"]


def _expected(m5, names: list[str], tmp_path) -> pl.DataFrame:
    """``live.bars`` of the picked genes, as ``ctx function --genes-from`` tests the list."""
    root = m5.parents[2]
    spec = load(DE_SPEC)
    (path,) = (root / "outputs" / "05_de").glob(f"*_desynth_{spec.de.featured}_results_v1.parquet")
    sets = lite.read_gmt(root / spec.resource("sets").dest)
    listed = tmp_path / "list.txt"
    listed.write_text("\n".join(names) + "\n")
    return live.bars(live.listed(listed, pl.read_parquet(path), sets))


def _legend(series: str) -> str:
    family, _, failing = series.partition("_")
    words = text.raw("function", "live.families")
    key = "live.legend.failing" if failing else "live.legend.passing"
    return text.fill("function", key, family=words[family], alpha=live.ALPHA)


def test_a_clicked_link_shows_the_genes_on_m6(ran):
    """The click goes through the entry page to M6's overview with the six genes; M6 names
    each and brings its live ORA into the top half of the window.

    Mutations caught: the reveal left out (the panel 2,000 px down, as Rytis saw it); the
    names left out of the notes; the entry page dropping the fragment."""
    _, _, results = ran
    seen = results["link"]
    assert seen["errors"] == []
    picked, got = seen["acted"], seen["moused"][0]
    assert picked["count"].startswith("6 genes selected")
    assert got["url"].split("#")[0].endswith("_function_desynth_overview_v1.html")
    assert genes.decode(got["hash"]).genes and "&n=6&" in got["hash"]
    assert got["notes"].startswith("6 genes of the list have a symbol on this page")
    assert sorted(_listed(got["notes"])) == sorted(picked["names"])
    assert got["summary"].startswith("The list of the link: 6 genes")
    assert _in_view(got), _where(got)


@pytest.mark.parametrize("name", list(SIZES))
def test_six_genes_draw_the_sets_of_smallest_p_marked_as_failing(ran, name, tmp_path):
    """Rytis's 6 genes (PK.3): the panel draws the sets ``live.bars`` draws, in its order, each
    failing set grey and named so in the legend, the summary saying how many and at what FDR;
    the PDF and SVG are the view. The length is -log10 FDR, past a dashed line at FDR 0.05
    (Rytis 2026-10-21, D-102). At 1,280 and 375 px, light and dark, no sideways scroll.

    Mutations caught: the failing sets dropped (the empty chart Rytis saw); a failing set in
    its family's colour; the length -log10 p in the page's redraw; the dashed line left out;
    the grey sentence left out."""
    m5, downloads, results = ran
    seen = results[name]
    assert seen["errors"] == []
    got = seen["moused"][0]
    want = _expected(m5, seen["acted"]["names"], tmp_path)
    failing = want.filter(~pl.col("passes"))
    assert failing.height > 0, "the fixture's six genes must leave a family short of passing"
    drawn = {
        (b["name"], label): x for b in got["bars"] for label, x in zip(b["y"], b["x"], strict=True)
    }
    assert len(drawn) == want.height
    for row in want.iter_rows(named=True):
        x = drawn[(_legend(row["series"]), set_words(row["gene_set"]))]
        assert math.isclose(x, -math.log10(row["fdr"]), rel_tol=1.52e-10, abs_tol=1e-12), row
    assert got["categories"] == list(dict.fromkeys(set_words(s) for s in want["gene_set"]))
    for bar in got["bars"]:
        grey = "not passing" in bar["name"]
        assert (bar["colour"] == got["grey"]) == grey, bar
    assert got["legend"] is True
    line = f"FDR {live.ALPHA}"
    assert got["shapes"] == [{"x0": pytest.approx(-math.log10(live.ALPHA)),
                              "x1": pytest.approx(-math.log10(live.ALPHA)), "xref": "x",
                              "dash": "dash", "name": line}], got["shapes"]  # fmt: skip
    passing = int(re.search(r"(\d+) set results pass FDR", got["summary"])[1])
    sentence = text.fill("function", "live.grey", failing=failing.height, alpha=live.ALPHA,
                         top=live.TOP)  # fmt: skip
    assert sentence in got["summary"], got["summary"]
    assert passing == want["passes"].sum(), got["summary"]
    assert got["moved"] and got["scroll"] <= got["box"], (got["scroll"], got["box"])
    files = seen["files"]
    assert files["svg"]["name"].endswith("_view.svg"), files
    svg = (downloads / files["svg"]["guid"]).read_text()
    bars = re.findall(r'<g class="trace bars".*?</g></g>', svg, re.DOTALL)
    assert any(GREY_SVG.search(b) for b in bars), "no grey bar in the view's SVG"
    legend = [html.unescape(s) for s in re.findall(r'class="legendtext"[^>]*>([^<]*)', svg)]
    assert legend == [b["name"] for b in got["bars"]] + [line], legend
    assert (downloads / files["pdf"]["guid"]).read_bytes().startswith(b"%PDF")


def test_open_beside_shows_the_genes_in_the_window_it_opens(ran):
    """Open beside, clicked: the named window opens M6 with the list, names the genes and shows
    the panel. Mutation caught: the reveal left out."""
    _, _, results = ran
    seen = results["beside"]
    assert seen["errors"] == []
    (got,) = (v for url, v in seen["others"].items() if "_function_desynth_overview" in url)
    assert got["opener"]
    assert sorted(_listed(got["notes"])) == sorted(seen["acted"]["names"])
    assert _in_view(got), _where(got)


def test_a_long_list_names_its_first_genes_and_counts_the_rest(ran):
    """Every up call sent: the notes name ``genes.PAGE_ROWS`` genes, then the rest by count.
    Mutation caught: the cap ignored, a list of thousands of names in one line."""
    got = ran[2]["many"]["moused"][0]
    n = int(re.match(r"(\d+) genes of the list", got["notes"])[1])
    assert n > genes.PAGE_ROWS, n
    shown = re.search(r"The list: (.*) and (\d+) more\.$", got["notes"])
    assert len(shown[1].split(", ")) == genes.PAGE_ROWS
    assert int(shown[2]) == n - genes.PAGE_ROWS


def test_a_refused_address_still_sends_the_list(ran):
    """A page whose address cannot be rewritten (a browser's file:// rule) still links the
    list and M6 receives it. Mutation caught: replaceState left unguarded, which threw
    before the onward link and the windows beside were updated."""
    seen = ran[2]["refused"]
    assert seen["errors"] == []
    got = seen["moused"][0]
    assert "&n=6&" in seen["acted"]["href"] and "&n=6&" in got["hash"]
    assert sorted(_listed(got["notes"])) == sorted(seen["acted"]["names"])


def _seen(results: dict, name: str) -> dict:
    """What a task's reader saw of M6: after the click, or as opened."""
    seen = results[name]
    assert seen["errors"] == []
    return seen["moused"][0] if name in SIZES else seen["acted"]


@pytest.mark.parametrize("name", ["built_phone_light", "built_phone_dark", "link_phone_light",
                                  "link_phone_dark"])  # fmt: skip
def test_a_phone_keeps_every_set_label_on_screen(ran, name):
    """At 375 px, M6 as built and redrawn on the six genes (PK.4): every y tick label starts
    at or right of the screen's left edge, one label per bar's set, no page scrolls sideways.
    Before, Step-1's build view set its labels from -201 px and the six genes' from -18 px;
    on this fixture from -52.4 px and -58.5 px.

    Mutation caught: the redraw's phone layout left out (``ctxBars.narrow`` returning null),
    the six genes' labels then from the left of the screen."""
    got = _seen(ran[2], name)
    assert len(got["ticks"]) == len(got["categories"]) > 0, (got["ticks"], got["categories"])
    off = [(t["text"], round(t["left"], 1)) for t in got["ticks"] if t["left"] < 0]
    assert off == [], off
    assert any(t["lines"] > 1 for t in got["ticks"]), "no label wrapped: no phone layout drawn"
    assert got["scroll"] <= got["box"], (got["scroll"], got["box"])


@pytest.mark.parametrize("name", ["built_wide", "link", "link_wide_dark"])
def test_a_desktop_labels_every_set(ran, name):
    """At 1,280 px every set has its own tick label, on screen, and no two overlap: the live
    bars switch to the phone's layout below twice the widest label plus margins, as fixed bars
    do (Rytis 2026-10-21, D-105), so M6's 567 px plot wraps its labels as the GSEA bars beside
    it do. PK.4 asserted the opposite, one line at the build's height: its threshold, the
    widest label plus ``HELD_PAD_PX``, kept plotly's own ticks there, which drew 7 of 14.

    Mutation caught: that threshold back in either twin, every other label then skipped."""
    got = _seen(ran[2], name)
    squeeze = [re.sub(r"\s", "", t["text"]) for t in got["ticks"]]
    assert sorted(squeeze) == sorted(re.sub(r"\s", "", c) for c in got["categories"]), squeeze
    assert all(t["left"] >= 0 for t in got["ticks"]), got["ticks"]
    assert any(t["lines"] > 1 for t in got["ticks"]), "no label wrapped: no phone layout drawn"
    rows = sorted(got["ticks"], key=lambda t: t["top"])
    overlaps = [(a["text"], b["text"]) for a, b in pairwise(rows) if b["top"] < a["bottom"]]
    assert overlaps == [], overlaps
    assert got["scroll"] <= got["box"], (got["scroll"], got["box"])
