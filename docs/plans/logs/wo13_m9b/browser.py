# ruff: noqa: E501 - page scripts are quoted whole, a line each
"""Package M9b, E.1c2 and E.1d: a built M9 page in headless Chromium at 1,280 and 375 px, both themes.

1. Every packed card drawn; every box and F1 plot drew its points; each F1 view chosen in turn
   and drawn; the page wider than the window nowhere.
2. A pasted list of 150 tested genes: the live heatmap draws the first 100 in list order, and
   each cell's step equals the one Python reads from M5's exact log2 FC and call; the table
   counts all 150.
3. F4 opened alone: its first tissue unticked empties that column.
Run: python browser.py <M9 page> <out.json>
"""

import json
import sys
import tempfile
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "tests"))
from test_export_browser import _browse  # noqa: E402

from ctx.genecard import listing  # noqa: E402
from ctx.genecard.listing import unblock  # noqa: E402

DRIVE = """(async () => {
  if (%(dark)s && document.documentElement.getAttribute("data-theme") !== "black") document.getElementById("ctx-theme-toggle").click();
  for (const b of document.querySelectorAll("[data-unpack]")) b.click();
  for (let i = 0; i < 600; i++) {
    const left = Array.from(document.querySelectorAll(".packed-body")).filter(h => h.getAttribute("data-state") !== "drawn");
    if (!left.length) break;
    await new Promise(r => setTimeout(r, 100));
  }
  await new Promise(r => setTimeout(r, 2000));
  const figs = Array.from(document.querySelectorAll(".fig"));
  const boxes = figs.filter(f => (f.querySelector(".fig-title") || {}).textContent && f.querySelector(".fig-title").textContent.includes(", M5 normalised values"));
  const drew = boxes.map(f => { const d = f.querySelector(".js-plotly-plot"); return d && d.data ? d.data.reduce((n, t) => n + (t.y ? t.y.length : 0), 0) : -1; });
  const views = [];
  for (const f of boxes) {
    const form = f.querySelector("form.view");
    if (!form) continue;
    const sel = form.elements.view, seen = [];
    for (let i = 0; i < sel.options.length; i++) {
      sel.value = String(i); sel.dispatchEvent(new Event("change", {bubbles: true}));
      await new Promise(r => setTimeout(r, 150));
      const d = f.querySelector(".js-plotly-plot");
      seen.push({label: sel.options[i].textContent, points: d.data.reduce((n, t) => n + (t.y ? t.y.length : 0), 0),
                 lead: f.querySelector(".fig-lead").textContent.slice(0, 80)});
    }
    sel.value = "0"; sel.dispatchEvent(new Event("change", {bubbles: true}));
    views.push(seen);
  }
  const form = document.querySelector("form.genelist");
  form.elements.pasted.value = %(pasted)s;
  form.querySelector("[data-paste]").click();
  await new Promise(r => setTimeout(r, 2000));
  const panel = document.querySelector(".genelist-panel"), heat = panel.querySelector(".list-heat-plot");
  const W = document.documentElement.clientWidth;
  return {theme: document.documentElement.getAttribute("data-theme"), width: W,
          fits: document.documentElement.scrollWidth <= W, boxes: boxes.length, drew: drew, views: views,
          heat: panel.ctxHeat || null, pages: panel.querySelector(".pick-pages span").textContent,
          heat_width: heat.getBoundingClientRect().width, heat_min: heat.style.minWidth,
          heat_keys: (heat.data || []).slice(1).map(t => t.name),
          heat_hint: !panel.querySelector(".list-heat .scroll-hint").hidden,
          heat_ticks: heat.layout ? heat.layout.yaxis.ticktext.length : 0,
          caption: panel.querySelector(".list-heat-caption").textContent,
          summary: panel.querySelector(".list-summary").textContent};
})()"""
ALONE = """(async () => {
  await new Promise(r => setTimeout(r, 2000));
  const form = document.querySelector("form.subset");
  if (!form) return null;
  const fig = form.closest(".fig"), boxes = Array.from(form.querySelectorAll("input[name=kept]"));
  boxes[0].checked = false;
  form.dispatchEvent(new Event("change"));
  const d = fig.querySelector(".js-plotly-plot");
  return {dropped: boxes[0].value, tissues: boxes.length, column: d.data[0].z.map(r => r[0]),
          summary: fig.querySelector(".subset-summary").textContent,
          fits: document.documentElement.scrollWidth <= document.documentElement.clientWidth};
})()"""


def expected(page: Path, numbers: list[int]) -> list[list[int | None]]:
    """Each listed gene's cell per contrast, read from M5's tables as the evidence view reads them."""
    root = page.parents[2]
    from ctx.dataset import load  # noqa: PLC0415
    from ctx.genecard import cards  # noqa: PLC0415

    slug = page.name.split("_genecard_", 1)[1].rsplit("_overview", 1)[0]
    spec = load(Path(__file__).resolve().parents[4] / "conf" / "datasets" / f"{slug}.toml")
    grey = len(listing.STEPS)
    out = []
    for n in numbers:
        row = []
        for c in spec.contrasts:
            frame = cards.results(root, slug, c.name).filter(pl.col("number") == n)
            if frame.is_empty() or frame["padj"][0] is None:
                row.append(None)
                continue
            fc, call = frame["log2FoldChange"][0], frame["call"][0]
            if call in listing.CALLED:
                from ctx.viz.heatmap import diverging_levels  # noqa: PLC0415

                row.append(listing.STEPS.index(diverging_levels([fc])[0][0]))
            else:
                row.append(grey)
        out.append(row)
    return out


def main(page: Path, out: Path) -> None:
    """Drive the page and its first F4 page alone; write what each did and the checks."""
    html = page.read_text()
    block = json.loads(html.split('id="ctx-panel-1">', 1)[1].split("</script>", 1)[0])["block"]
    data = unblock(block)
    numbers = data["numbers"][:150]
    pasted = " ".join(f"ENSG{n:011d}" for n in numbers)
    alone = sorted(page.parent.glob(f"{page.name.split('_overview')[0]}_exons_*_v1.html"))[0]
    tasks = []
    for width in (1280, 375):
        for dark in (False, True):
            tasks.append(
                {
                    "name": f"page_{width}_{'dark' if dark else 'light'}",
                    "url": page.as_uri(),
                    "width": width,
                    "height": 900,
                    "scheme": "dark" if dark else "light",
                    "act": DRIVE
                    % {"dark": "true" if dark else "false", "pasted": json.dumps(pasted)},
                }
            )
        tasks.append(
            {
                "name": f"alone_{width}",
                "url": alone.as_uri(),
                "width": width,
                "height": 900,
                "act": ALONE,
            }
        )
    with tempfile.TemporaryDirectory() as downloads:
        results = _browse(tasks, Path(downloads), timeout=1800)
    want = expected(page, numbers[:100])
    report = {}
    for name, r in results.items():
        got = r["acted"]
        checks = {"errors": r["errors"]}
        if name.startswith("page"):
            checks |= {
                "theme": got["theme"],
                "fits": got["fits"],
                "boxes": got["boxes"],
                "every_box_drew": all(x > 0 for x in got["drew"]),
                "views": len(got["views"]),
                "view_options": sorted({len(v) for v in got["views"]}),
                "every_view_drew": all(s["points"] > 0 for v in got["views"] for s in v),
                "heat_rows": len(got["heat"]["genes"]),
                "heat_equal": got["heat"]["z"] == want,
                "heat_keys": got["heat_keys"],
                "heat_width": got["heat_width"],
                "heat_min": got["heat_min"],
                "heat_hint": got["heat_hint"],
                "pages": got["pages"],
                "caption": got["caption"],
                "summary": got["summary"],
                "first_view": got["views"][0] if got["views"] else None,
            }
        else:
            checks |= {
                "tissues": got["tissues"],
                "column_empty": all(c is None for c in got["column"]),
                "summary": got["summary"],
                "fits": got["fits"],
            }
        report[name] = checks
    out.write_text(json.dumps(report, indent=1, ensure_ascii=False))
    print(json.dumps(report, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main(Path(sys.argv[1]).absolute(), Path(sys.argv[2]))
