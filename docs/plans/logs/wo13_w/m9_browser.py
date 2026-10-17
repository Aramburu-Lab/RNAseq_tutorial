# ruff: noqa: E501 - page scripts and mutations are quoted whole, a line each
"""Package W, E.1a: a built M9 page, F3 and F4 against Python, and a list from M5's picker.

In headless Chromium, the list by link and by Open beside.

1. On the M9 page at 1,280 and 375 px in both themes: every packed card drawn, then in each
   F3 the first lineage left out and tau taken on nTPM, in each F4 the first tissue left out;
   the page's numbers are compared with ``figures.immune_subset`` and ``annotation.exon_use``
   on the rows the figure embeds, and those rows with the page's own tables.
2. On the M5 page: "Select the current up calls", then the genecard link as it stands and
   Open beside; the link is opened, and the window Open beside opened is read.
Run: python m9_browser.py <M9 page> <M5 page> <out.json>
"""

import json
import math
import sys
import tempfile
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "tests"))
from test_export_browser import _browse  # noqa: E402

from ctx.genecard import annotation, figures  # noqa: E402
from ctx.genecard.listing import unblock  # noqa: E402
from ctx.report.genes import decode  # noqa: E402
from ctx.viz.heatmap import diverging_levels  # noqa: E402

MOVE = """(async () => {
  if (%s) document.getElementById("ctx-theme-toggle").click();
  for (const b of document.querySelectorAll("[data-unpack]")) b.click();
  for (let i = 0; i < 300; i++) {
    const left = Array.from(document.querySelectorAll(".packed-body"))
      .filter(h => h.getAttribute("data-state") !== "drawn");
    if (!left.length) break;
    await new Promise(r => setTimeout(r, 100));
  }
  await new Promise(r => setTimeout(r, 1500));
  const out = [];
  for (const form of document.querySelectorAll("form.subset")) {
    const fig = form.closest(".fig"), d = fig.querySelector(".js-plotly-plot");
    const cfg = JSON.parse(document.getElementById(form.getAttribute("data-controls")).textContent);
    const t = JSON.parse(document.getElementById(form.getAttribute("data-rows")).textContent);
    const box = form.querySelector("input[name=kept]");
    box.checked = false;
    if (form.elements.scale) form.elements.scale.value = "linear";
    form.dispatchEvent(new Event("change"));
    out.push({twin: cfg.twin, title: fig.querySelector(".fig-title").textContent,
      dropped: box.value, options: cfg.options, exons: cfg.exons || null, rows: t,
      got: fig.ctxM9, lead: fig.querySelector(".fig-lead").textContent,
      summary: fig.querySelector(".subset-summary").textContent,
      visible: d.data.map(x => x.visible === undefined ? true : x.visible),
      categories: d.layout.yaxis.categoryarray ? d.layout.yaxis.categoryarray.slice() : null,
      z: cfg.twin === "exons" ? d.data[0].z : null,
      ticks: d.layout.yaxis2 ? d.layout.yaxis2.ticktext : null});
  }
  return {figures: out, theme: document.documentElement.getAttribute("data-theme"),
          width: document.documentElement.clientWidth,
          fits: document.documentElement.scrollWidth <= document.documentElement.clientWidth};
})()"""
PICK = """(async () => {
  await new Promise(r => setTimeout(r, 3000));
  document.querySelector("button[data-pick=up]").click();
  await new Promise(r => setTimeout(r, 500));
  const a = Array.from(document.querySelectorAll("a[data-onward]"))
    .filter(x => x.getAttribute("href").indexOf("11_genecard") >= 0)[0];
  const href = new URL(a.getAttribute("href"), location.href).href;
  a.parentElement.querySelector("[data-beside]").click();
  await new Promise(r => setTimeout(r, 8000));
  return {href: href, text: a.textContent};
})()"""
LISTED = """(async () => {
  /* The browser's first tab, about:blank, sits in the background, where timers run once a
     minute: read nothing there. */
  if (location.protocol !== "file:") return {url: location.href, list: null};
  let panel = null;
  for (let i = 0; i < 300 && !(panel = document.querySelector(".genelist-panel")); i++) {
    await new Promise(r => setTimeout(r, 100));
  }
  if (!panel) return {url: location.href, list: null};
  for (let i = 0; i < 300 && !panel.ctxList; i++) await new Promise(r => setTimeout(r, 100));
  return {url: location.href, list: panel.ctxList || null,
          summary: panel.querySelector(".list-summary").textContent,
          rows: Array.from(panel.querySelectorAll("tbody tr")).map(tr =>
            Array.from(tr.children).map(td => td.textContent)),
          pages: panel.querySelector(".pick-pages span").textContent,
          command: panel.querySelector(".list-command code").textContent.slice(0, 300),
          fits: document.documentElement.scrollWidth <= document.documentElement.clientWidth};
})()"""


def _rows(table: dict) -> list[dict]:
    """A figure's embedded rows, one dict each, as ctxRows makes them."""
    return [dict(zip(table["columns"], values, strict=True)) for values in table["data"]]


def _close(a: float | None, b: float | None) -> bool:
    if a is None or b is None:
        return a is b
    return abs(a - b) <= 1e-12 * max(1.0, abs(a))


def immune(fig: dict) -> list[str]:
    """What differs between the page's F3 and Python's on the same rows; empty if nothing."""
    rows = [(r["cell_type"], r["lineage"], r["ntpm"]) for r in _rows(fig["rows"])]
    kept = set(fig["options"]) - {fig["dropped"]}
    want = figures.immune_subset(rows, kept, log=False)
    bad = []
    if (fig["got"]["cells"], fig["got"]["top"]) != (want["cells"], want["top"]):
        bad.append(f"cells or top {fig['got']} {want}")
    if not _close(fig["got"]["tau"], want["tau"]):
        bad.append(f"tau {fig['got']['tau']} {want['tau']}")
    if fig["categories"] != [c for c, lineage, _ in rows if lineage in kept]:
        bad.append("categories")
    return bad


def exons(fig: dict) -> list[str]:
    """What differs between the page's F4 and Python's on the same rows; empty if nothing."""
    frame = pl.DataFrame(_rows(fig["rows"])).select("exon", "tissue", "share")
    kept = [t for t in fig["options"] if t != fig["dropped"]]
    use = annotation.exon_use(frame, kept)
    levels, _ = diverging_levels(list(use.steps))
    names = [s.level for s in diverging_levels([])[1]]
    want = dict(zip(use.frame.select("exon", "tissue").rows(), levels, strict=True))
    bad = []
    for i, e in enumerate(fig["exons"]):
        for j, t in enumerate(fig["options"]):
            cell = fig["z"][i][j]
            shown = names[cell] if cell is not None else None
            if shown != (None if t == fig["dropped"] else want[(e, t)]):
                bad.append(f"cell {e} {t}: {shown} {want.get((e, t))}")
    cv = {e: (None if v is None else v) for e, v in use.cv.items()}
    if any(not _close(fig["got"]["cv"].get(e), v) for e, v in cv.items()):
        bad.append("cv")
    pair = annotation.qpcr_pair(use.cv)
    got = fig["got"]["pair"]
    if (got is None) != (pair is None) or (pair and (got["first"], got["second"]) !=
                                           (pair.first, pair.second)):  # fmt: skip
        bad.append(f"pair {got} {pair}")
    worst = figures.least_even(use.cv)
    if pair and worst and f"{worst[1]:.2f}" not in fig["lead"]:
        bad.append(f"lead {fig['lead']}")
    if fig["ticks"] is not None and fig["ticks"] != [
            "" if use.cv[e] is None else f"{use.cv[e]:.2f}" for e in fig["exons"]]:  # fmt: skip
        bad.append("ticks")
    return bad


def main() -> None:
    """Both parts; the summary to stdout and the raw results to the out file."""
    # absolute(), not resolve(): a linked M5 folder resolves into the real root, whose links
    # reach the real M9 page.
    m9, m5, out = Path(sys.argv[1]).absolute(), Path(sys.argv[2]).absolute(), Path(sys.argv[3])
    tasks = [{"name": f"move_{w}_{'black' if dark else 'white'}", "url": m9.as_uri(), "width": w,
              "height": 900, "act": MOVE % ("true" if dark else "false")}
             for w in (1280, 375) for dark in (False, True)]  # fmt: skip
    tasks.append({"name": "pick", "url": m5.as_uri(), "width": 1280, "height": 900, "act": PICK,
                  "others": LISTED})  # fmt: skip
    with tempfile.TemporaryDirectory(prefix="w_m9_") as downloads:
        ran = _browse(tasks, Path(downloads))
        href = ran["pick"]["acted"]["href"]
        linked = _browse([{"name": "link", "url": href, "width": 375, "height": 812,
                           "until": "_overview_v1.html", "act": LISTED}],
                         Path(downloads))["link"]  # fmt: skip
    summary: dict = {"errors": {k: v["errors"] for k, v in ran.items() if v["errors"]}}
    for name, r in ran.items():
        if not name.startswith("move"):
            continue
        figs = r["acted"]["figures"]
        summary[name] = {"theme": r["acted"]["theme"], "width": r["acted"]["width"],
                         "fits": r["acted"]["fits"],
                         "f3": [len(immune(f)) for f in figs if f["twin"] == "immune"],
                         "f4": [len(exons(f)) for f in figs if f["twin"] == "exons"],
                         "problems": [p for f in figs
                                      for p in (immune(f) if f["twin"] == "immune" else exons(f))],
                         "tau_defined": sum(f["got"]["tau"] is not None for f in figs
                                            if f["twin"] == "immune")}  # fmt: skip
    sent = decode("#" + href.split("#", 1)[1])
    page = m9.read_text()
    block = json.loads(page.split('id="ctx-panel-1">', 1)[1].split("</script>", 1)[0])
    tested = set(unblock(block["block"])["numbers"])
    numbers = [abs(g) for g in sent.genes]
    summary["link"] = {"sent": len(numbers), "tested": sum(n in tested for n in numbers),
                       "carded": sum(c[0] in set(numbers) for c in block["cards"]),
                       "page": linked["acted"]["list"] and {k: (len(v) if isinstance(v, list) else v)
                                                            for k, v in linked["acted"]["list"].items()},
                       "summary": linked["acted"]["summary"], "pages": linked["acted"]["pages"],
                       "first_rows": linked["acted"]["rows"][:3], "fits": linked["acted"]["fits"],
                       "errors": linked["errors"], "url": linked["acted"]["url"][:120]}  # fmt: skip
    beside = ran["pick"].get("others", {})
    summary["beside"] = {url[:100]: (v.get("list") and {k: (len(x) if isinstance(x, list) else x)
                                                        for k, x in v["list"].items()},
                                     v.get("summary")) for url, v in beside.items()}  # fmt: skip
    out.write_text(json.dumps(ran, indent=1, default=str))
    print(json.dumps(summary, indent=1, default=lambda v: None if isinstance(v, float)
                     and not math.isfinite(v) else str(v)))  # fmt: skip


if __name__ == "__main__":
    main()
