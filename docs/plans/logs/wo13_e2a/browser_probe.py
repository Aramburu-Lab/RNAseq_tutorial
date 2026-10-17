"""Package wo13-e2a: a page's controls in headless Chromium at 1280 and 375 px, both themes.

For each figure whose title holds the given text, the probe moves one control (a select set
to a value, or a script), reads what the figure shows, checks the page does not scroll
sideways, presses the figure's PDF and SVG buttons and checks the files: the SVG holds no
raster image and carries the live line; the PDF opens and is the figure files' size.

    python browser_probe.py <page.html> <out.json> <title text> <control selector> <value>
           [<title text> <control selector> <value> ...]

<control selector> is a CSS selector inside the figure, such as 'form.view select'; a
<value> of "-" moves nothing and saves no file, for a page without the control.
"""

import json
import sys
from pathlib import Path

WT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(WT / "tests"))
from test_export_browser import CHROMIUM, _browse  # noqa: E402

page, out = Path(sys.argv[1]).resolve(), Path(sys.argv[2])
asks = [sys.argv[i : i + 3] for i in range(3, len(sys.argv), 3)]

ACT = """(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(2500);
  const want = %s, selector = %s, value = %s;
  const fig = Array.from(document.querySelectorAll('figure.fig')).find(f =>
    ((f.querySelector('.fig-title') || {}).textContent || '').includes(want));
  if (!fig) return {missing: want};
  /* A figure in a hidden tab: open its tab first, as a reader would. */
  const panel = fig.closest('[role=tabpanel]');
  if (panel && panel.hidden) {
    const tab = document.querySelector('[role=tab][aria-controls="' + panel.id + '"]');
    if (tab) tab.click();
    await wait(500);
  }
  fig.scrollIntoView();
  await wait(500);
  fig.setAttribute('data-probe', '1');
  const d = fig.querySelector('.js-plotly-plot');
  const read = () => ({
    traces: d.data.length, data: JSON.stringify(d.data.map(t => [t.x, t.y, t.name])).length,
    caption: (fig.querySelector('figcaption') || {}).textContent,
    lead: (fig.querySelector('.fig-lead') || {}).textContent,
    line: Array.from(fig.querySelectorAll('[aria-live]')).map(e => e.textContent).join(' | '),
    moved: !!fig.querySelector('a[data-export][role=button]'),
    pageWide: document.documentElement.scrollWidth, window: window.innerWidth,
    plotWide: d.getBoundingClientRect().width,
    boxWide: fig.querySelector('.fig-plot').getBoundingClientRect().width,
    theme: document.documentElement.getAttribute('data-theme'),
    paper: d._fullLayout.paper_bgcolor, font: d._fullLayout.font.color,
  });
  const out = {start: read()};
  const control = fig.querySelector(selector);
  if (!control) return {...out, noControl: selector};
  control.value = value;
  control.dispatchEvent(new Event('change', {bubbles: true}));
  await wait(800);
  out.after = read();
  return out;
})()"""


def task(name: str, width: int, scheme: str, ask: list[str]) -> dict:
    """One driver task: the page at a width and scheme, one control moved, its files saved."""
    want, selector, value = ask
    return {"name": name, "url": page.as_uri(), "width": width, "height": 900,
            "scheme": scheme, "act": ACT % (json.dumps(want), json.dumps(selector),
                                            json.dumps(value)),
            **({"figure": "[data-probe]"} if value != "-" else {})}  # fmt: skip


def check_files(found: dict, downloads: Path) -> dict:
    """What the saved files hold: raster images, text, the size of the PDF's page."""
    from pypdf import PdfReader  # noqa: PLC0415

    seen = {}
    for ext, saved in found.get("files", {}).items():
        path = downloads / saved["guid"]
        data = path.read_bytes() if path.exists() else b""
        if ext == "svg":
            text = data.decode("utf-8", "replace")
            seen["svg"] = {"bytes": len(data), "images": text.count("<image"),
                           "name": saved["name"],
                           "has_view_line": "This view shows" in text}  # fmt: skip
        else:
            reader = PdfReader(path)
            box = reader.pages[0].mediabox
            seen["pdf"] = {"bytes": len(data), "pages": len(reader.pages), "name": saved["name"],
                           "pt": [float(box.width), float(box.height)]}  # fmt: skip
    return seen


assert CHROMIUM, "no Chromium; set CTX_CHROMIUM"
downloads = WT / "outputs" / "scratch_roots" / "downloads" / (out.stem + "_downloads")
downloads.mkdir(parents=True, exist_ok=True)
tasks = [task(f"{i}:{w}:{s}", w, s, ask) for i, ask in enumerate(asks)
         for w in (1280, 375) for s in ("light", "dark")]  # fmt: skip
results = _browse(tasks, downloads)
report = {}
for name, found in results.items():
    acted = found.get("acted")
    report[name] = {"errors": found.get("errors"), "acted": acted,
                    "files": check_files(found, downloads)}  # fmt: skip
out.write_text(json.dumps(report, indent=1))
for name, r in report.items():
    a = r["acted"] or {}
    start, after = a.get("start", {}), a.get("after", {})
    print(
        name,
        "errors",
        len(r["errors"] or []),
        "moved",
        after.get("moved"),
        "wide",
        after.get("pageWide"),
        "/",
        after.get("window"),
        "plot",
        round(after.get("plotWide") or 0),
        "/",
        round(after.get("boxWide") or 0),
        "theme",
        after.get("theme"),
        "traces",
        start.get("traces"),
        "->",
        after.get("traces"),
        "files",
        r["files"],
    )
