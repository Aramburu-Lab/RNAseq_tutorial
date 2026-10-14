"""The M5 gene picker in Chromium (GS.4, GS.6): a real lasso and box drag, at 1,280 and 375 px,
in both themes; the onward link and Open beside, both ways, through an entry page.

The drags are real mouse events (tests/browser.mjs, ``mouse``). Each test names the mutation
it catches. Skipped where no Chromium is found; set ``CTX_CHROMIUM``.
"""

from dataclasses import replace
from datetime import date

import polars as pl
import pytest
from synth import counts, demo_volcano
from test_export_browser import CHROMIUM, _browse

from ctx.io.naming import stem
from ctx.report import genes
from ctx.report.builder import (
    FigureText,
    Group,
    ModulePage,
    ParameterCard,
    Reproduction,
    SelfCheck,
    write,
)
from ctx.report.entry import write_entry
from ctx.viz.save import figure_block, write_figure
from ctx.viz.style import SELECTED, THEMES, contrast

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
DAY = date(2026, 10, 14)
SRC = "demo:case_control"
WAIT = "await new Promise(r => setTimeout(r, {ms}));"
# The volcano's plot box in the window, after scrolling it into view and choosing the drag.
SETUP = """(async () => {{
  {wait}
  const fig = document.querySelector('form.select').closest('.fig');
  fig.querySelector('input[name=drag][value={drag}]').click();
  fig.querySelector('.js-plotly-plot').scrollIntoView({{block: 'center', behavior: 'instant'}});
  await new Promise(r => setTimeout(r, 600));
  return null;
}})()"""
# A path over the plot's upper right: a box's two corners, or a lasso's polygon.
PATH = """(() => {{
  const d = document.querySelector('form.select').closest('.fig').querySelector('.js-plotly-plot');
  const r = d.getBoundingClientRect(), s = d._fullLayout._size;
  const at = (fx, fy) => [r.left + s.l + fx * s.w, r.top + s.t + fy * s.h];
  const corners = {shape}, xa = d._fullLayout.xaxis, ya = d._fullLayout.yaxis;
  // In data units: the page may scroll between the drag and the reading.
  window.ctxPath = corners.map(c => [xa.p2l(c[0] - r.left - s.l), ya.p2l(c[1] - r.top - s.t)]);
  const out = [['mousePressed', ...corners[0]]];
  for (let i = 1; i < corners.length; i++) {{
    for (let k = 1; k <= 8; k++) {{
      const a = corners[i - 1], b = corners[i];
      out.push(['mouseMoved', a[0] + (b[0] - a[0]) * k / 8, a[1] + (b[1] - a[1]) * k / 8]);
    }}
  }}
  out.push(['mouseReleased', ...corners[corners.length - 1]]);
  return out;
}})()"""
BOX = "[at(0.55, 0.05), at(0.98, 0.05), at(0.98, 0.6)]"
LASSO = "[at(0.55, 0.05), at(0.98, 0.05), at(0.98, 0.6), at(0.55, 0.6), at(0.55, 0.06)]"
# What the drag picked against what lies inside the path, by point, from the drawn traces.
READ = """(async () => {
  const EDGE = 4;
  await new Promise(r => setTimeout(r, 900));
  const fig = document.querySelector('form.select').closest('.fig');
  const d = fig.querySelector('.js-plotly-plot');
  const xa = d._fullLayout.xaxis, ya = d._fullLayout.yaxis;
  const px = (p) => [xa.l2p(p[0]), ya.l2p(p[1])];
  let path = window.ctxPath.map(px);
  if (path.length === 3) {
    const [a, , c] = path;
    path = [a, [c[0], a[1]], c, [a[0], c[1]]];
  }
  const inside = (p) => {
    let c = false;
    for (let i = 0, j = path.length - 1; i < path.length; j = i++) {
      const [xi, yi] = path[i], [xj, yj] = path[j];
      if ((yi > p[1]) !== (yj > p[1]) && p[0] < (xj - xi) * (p[1] - yi) / (yj - yi) + xi) c = !c;
    }
    return c;
  };
  // Within EDGE px of the path a point may fall either way: plotly tests its own drawn path.
  const near = (p) => path.some((a, i) => {
    const b = path[(i + 1) % path.length], dx = b[0] - a[0], dy = b[1] - a[1];
    const along = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / (dx * dx + dy * dy || 1);
    const t = Math.max(0, Math.min(1, along));
    return Math.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy) < EDGE;
  });
  const ring = d.data[d.data.length - 1], want = new Set(), edge = new Set(), got = new Set();
  d.data.slice(0, -1).forEach(t => (t.customdata || []).forEach((c, i) => {
    if (t.mode !== 'markers') return;
    const at = px([t.x[i], t.y[i]]);
    if (near(at)) edge.add(c[0]); else if (inside(at)) want.add(c[0]);
  }));
  ring.x.forEach((x, i) => got.add(`${x},${ring.y[i]}`));
  const points = new Set();
  d.data.slice(0, -1).forEach(t => (t.customdata || []).forEach((c, i) => {
    if (got.has(`${t.x[i]},${t.y[i]}`)) points.add(c[0]);
  }));
  const link = fig.querySelector('a[data-onward]');
  return {want: [...want].sort(), edge: [...edge].sort(), got: [...points].sort(),
          rings: ring.x.length,
          count: fig.querySelector('.select-count').textContent, href: link.getAttribute('href'),
          link: link.textContent, hash: location.hash,
          checked: fig.querySelectorAll('table.pick tbody input:checked').length,
          moved: !!fig.querySelector('a[data-export][role=button]'),
          ringColour: ring.marker.color, font: d._fullLayout.font.color,
          scroll: document.documentElement.scrollWidth <= document.documentElement.clientWidth};
})()"""


def _volcano():
    """The synthetic volcano with Ensembl IDs and symbols, and the picker of a dataset that
    sends lists to M6."""
    data = demo_volcano(counts())
    frame = data.frame.with_row_index("i").with_columns(
        gene_id=pl.format("ENSG{}", pl.col("i").cast(pl.String).str.zfill(11)).str.replace(
            "ENSG00000000000", "ENSG00000099999"),
        gene=pl.format("G{}", pl.col("i")),
    ).drop("i")  # fmt: skip
    columns = {**data.columns, "gene": ("none", "Synthetic symbol.")}
    hover = ("gene_id", "gene", "padj")
    select = genes.picker(
        SRC,
        hover,
        [("../06_function/index.html", "M6 Function")],
        {"id": "gene_id", "symbol": "gene", "lfc": "lfc", "padj": "padj"},
    )
    return replace(data, frame=frame, columns=columns, hover=hover, select=select)


def _page(folder, name):
    data = _volcano()
    folder.mkdir(parents=True, exist_ok=True)
    figure_stem = stem("de", f"{name}_volcano", 1, day=DAY)
    text = FigureText("Volcano", tab="Volcano")
    write_figure(data, folder, figure_stem, title="Volcano", data_releases={}, text=text)
    built = ModulePage(
        title=f"M5 Differential expression: {name}", purpose="Purpose.", input="Input.",
        method="Method.", parameters=(ParameterCard("--x", "A card.", "1", "Moves it."),),
        results=(Group("Volcano", (figure_block(data, figure_stem, text),)),),
        interpretation="Interpretation.", questions=(SelfCheck("A question.", "An answer."),),
        pitfalls=("A pitfall.",),
        reproduce=Reproduction("ctx demo", 1.0, 1, 1, "laptop", "none", {"python": "3.12"}),
    )  # fmt: skip
    return write(built, folder / f"{stem('de', f'{name}_overview', 1, day=DAY)}.html")


def _drag(name, page, width, drag, shape, dark=False):
    """One drag task. The browser's profile is shared by every task, so each sets its theme,
    then loads the page again: the theme is read before first paint."""
    theme = "black" if dark else "white"
    first = f"localStorage.setItem('ctx-theme', '{theme}'); location.reload(); null"
    setup = SETUP.format(wait=WAIT.format(ms=2500), drag=drag)
    mouse = [{"at": "(async () => { await " + setup + "; return []; })()"},
             {"at": PATH.format(shape=shape), "then": READ}]  # fmt: skip
    return {
        "name": name,
        "url": page.as_uri(),
        "width": width,
        "height": 900,
        "act": first,
        "mouse": mouse,
    }


@pytest.fixture(scope="module")
def ran(tmp_path_factory):
    root = tmp_path_factory.mktemp("pick") / "outputs"
    m5 = _page(root / "05_de", "demo")
    _page(root / "06_function", "demo")  # the receiving page, a picker too
    # A receiving module's overview, as M6 names it: the entry page links it.
    (root / "06_function" / f"{stem('function', 'demo_overview', 1, day=DAY)}.html").write_bytes(
        next((root / "06_function").glob("*_de_demo_overview_v1.html")).read_bytes()
    )
    write_entry(root / "06_function", "function")
    tasks = [_drag(f"{kind}_{w}{'_dark' if dark else ''}", m5, w, drag, shape, dark)
             for kind, drag, shape in (("box", "select", BOX), ("lasso", "lasso", LASSO))
             for w in (1280, 375) for dark in (False, True)]  # fmt: skip
    fig = "const f = document.querySelector('form.select').closest('.fig');"
    beside = {
        "name": "beside", "url": m5.as_uri(), "width": 1280, "height": 900,
        "act": "(async () => { " + WAIT.format(ms=2500)
               + " document.querySelector('[data-beside]').click(); " + WAIT.format(ms=3000)
               + " await " + SETUP.format(wait="", drag="select") + "; return null; })()",
        "mouse": [{"at": PATH.format(shape=BOX), "then": READ}],
        "others": "(async () => { " + WAIT.format(ms=1500) + fig
                  + " const before = f.querySelector('.select-count').textContent;"
                  + " f.querySelector('[data-pick=down]').click(); " + WAIT.format(ms=800)
                  + " return {before, opener: !!window.opener, href: location.href,"
                  + " problems: f.querySelector('.select-problems').textContent}; })()",
        "reset": "(async () => { " + WAIT.format(ms=1500) + fig
                 + " return {count: f.querySelector('.select-count').textContent,"
                 + " problems: f.querySelector('.select-problems').textContent}; })()",
    }  # fmt: skip
    return m5, _browse([*tasks, beside], tmp_path_factory.mktemp("downloads"))


@pytest.mark.parametrize("kind", ["box", "lasso"])
@pytest.mark.parametrize("width", [1280, 375])
@pytest.mark.parametrize("dark", [False, True])
def test_a_real_drag_picks_the_genes_inside_it(ran, kind, width, dark):
    """The genes ringed are the drawn points inside the drag's path; the count, the table,
    the fragment and the onward link follow; the export switches to the view; the page does
    not scroll sideways.

    Mutations caught: plotly_selected not wired; a point named by its curve's index instead
    of its gene; the onward link left without the fragment; the export switch left off."""
    _, results = ran
    seen = results[f"{kind}_{width}{'_dark' if dark else ''}"]
    assert seen["errors"] == []
    (got,) = seen["moused"][-1:]
    inside, edge, picked = set(got["want"]), set(got["edge"]), set(got["got"])
    assert inside and inside <= picked <= inside | edge, (len(picked), len(inside), len(edge))
    n = len(picked)
    assert got["rings"] == n and got["count"].startswith(f"{n} genes selected")
    assert got["href"].startswith("../06_function/index.html#v=1&src=demo:case_control&n=")
    assert got["hash"] == got["href"].split("index.html")[1]
    assert got["moved"] and got["scroll"]
    assert got["link"].startswith(f"Open these {n} genes in M6 Function")
    assert got["ringColour"] == SELECTED


@pytest.mark.parametrize("width", [1280, 375])
def test_the_dark_theme_draws_the_picked_view_in_its_own_ink(ran, width):
    """Mutation: the picker's relayouts dropping the theme's template."""
    _, results = ran
    fonts = [results[f"box_{width}{s}"]["moused"][-1]["font"] for s in ("", "_dark")]
    assert fonts[0] != fonts[1], fonts


def test_the_ring_reads_on_both_themes():
    """The selected ring at 3:1 or more against each theme's background (WCAG graphics)."""
    for theme in THEMES.values():
        assert contrast(SELECTED, theme.background) >= 3  # noqa: PLR2004 - 3:1


def test_open_beside_sends_both_ways_through_the_entry_page(ran):
    """The window Open beside opens, through the entry page, hears the box drag; a change
    there comes back as Replace.

    Mutation caught: no message sent. The opener survives the entry page's replace here."""
    m5, results = ran
    seen = results["beside"]
    assert seen["errors"] == []
    (other,) = (v for url, v in seen["others"].items() if "06_function" in url)
    picked = seen["moused"][0]
    n = len(picked["want"])
    assert other["opener"] and "06_function" in other["href"]
    assert other["before"].startswith(f"{n} genes selected")
    assert "The list changed on another page" in other["problems"]
    assert "The list changed on another page" in seen["reset"]["problems"]
    assert seen["reset"]["count"] != picked["count"]
