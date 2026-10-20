"""PK: Rytis's steps on a kit's Step-1 M5 page, opened as file://, in one browser.

Ticks six genes in the picker's table, then clicks each onward link and each Open beside
button with a real mouse, and records what the page that opens receives.

Usage: python wo13_pk_repro.py <M5 overview page> <browser binary> <out.json>
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

DRIVER = Path(__file__).resolve().parents[3] / "tests" / "browser.mjs"
TICK = """(async () => {
  await new Promise(r => setTimeout(r, 6000));
  const fig = document.querySelector('form.select').closest('.fig');
  // The table redraws on every tick, so each box is found again.
  for (let i = 0; i < 6; i++) {
    fig.querySelectorAll('table.pick tbody input[type=checkbox]')[i].click();
  }
  await new Promise(r => setTimeout(r, 1500));
  const links = [...fig.querySelectorAll('a[data-onward]')];
  return {count: fig.querySelector('.select-count').textContent, hash: location.hash,
          links: links.map(a => [a.textContent, a.getAttribute('href'), a.target, a.rel])};
})()"""
# A real click at the centre of the element the selector names.
CLICK = """(() => {{
  const fig = document.querySelector('form.select').closest('.fig');
  const el = fig.querySelectorAll({sel})[{i}];
  el.scrollIntoView({{block: 'center', behavior: 'instant'}});
  const r = el.getBoundingClientRect(), x = r.left + r.width / 2, y = r.top + r.height / 2;
  return [['mousePressed', x, y], ['mouseReleased', x, y]];
}})()"""
# What the page that opened holds: its address, its fragment and what its receiver said.
READ = """(async () => {
  await new Promise(r => setTimeout(r, 6000));
  const text = (s) => [...document.querySelectorAll(s)].map(e => e.textContent.trim());
  let session = null, local = null;
  try { session = Object.keys(sessionStorage); } catch (e) { session = String(e); }
  try { local = Object.keys(localStorage); } catch (e) { local = String(e); }
  return {href: location.href.split('#')[0], hash: location.hash.slice(0, 120),
          hashLength: location.hash.length, opener: !!window.opener, origin: location.origin,
          history: history.length, session, local, ora: text('.ora-notes'),
          oraSummary: text('.ora-summary'), m9: text('.list-notes'),
          m9Summary: text('.list-summary'), top: scrollY, window: innerHeight,
          oraAt: [...document.querySelectorAll('form.ora')].map(f => f.closest('.fig')
            .getBoundingClientRect().top + scrollY),
          bars: [...document.querySelectorAll('form.ora')].map(f => f.closest('.fig')
            .querySelector('.js-plotly-plot').data.filter(t => t.type === 'bar')
            .reduce((n, t) => n + (t.y || []).length, 0)),
          m9At: [...document.querySelectorAll('.genelist-panel')].map(p =>
            p.getBoundingClientRect().top + scrollY),
          m9Rows: [...document.querySelectorAll('.genelist-panel tbody tr')].map(r =>
            r.cells[0].textContent),
          body: document.body.innerText.slice(0, 300)};
})()"""


def tasks(page: str) -> list[dict]:
    """The click tasks: each onward link and each Open beside button, in a fresh tab."""
    out = []
    for i, until in enumerate(("06_function/20", "11_genecard/20")):
        link = {"at": CLICK.format(sel="'a[data-onward]'", i=i), "until": until, "then": READ}
        out.append({"name": f"link_{i}", "url": page, "act": TICK, "mouse": [link]})
        beside = {"at": CLICK.format(sel="'[data-beside]'", i=i), "until": "/", "then": READ}
        out.append({"name": f"beside_{i}", "url": page, "act": TICK, "mouse": [beside],
                    "others": READ})  # fmt: skip
    return out


def main() -> None:
    """Run the tasks in the named browser and write what each page received."""
    m5, browser, out = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
    page = m5.absolute().as_uri()
    with tempfile.TemporaryDirectory() as downloads:
        plan = {
            "chromium": browser,
            "downloads": downloads,
            "tasks": tasks(page),
            "trace": True,
            "headed": os.environ.get("PK_HEADED") == "1",
        }
        node = shutil.which("node")
        assert node, "node not on PATH"
        done = subprocess.run([node, str(DRIVER)], input=json.dumps(plan), capture_output=True,
                              text=True, timeout=900, check=False, env=os.environ)  # fmt: skip
    if done.returncode:
        sys.exit(done.stderr[-4000:])
    results = json.loads(done.stdout)
    for r in results:
        r.pop("requests", None)
    out.write_text(json.dumps(results, indent=1))
    for r in results:
        print(r["name"], r["acted"]["count"], "errors:", r["errors"])
        for p in [*r["moused"], *r.get("others", {}).values()]:
            if "href" in p and "/outputs/" in p["href"]:
                print(
                    "  ",
                    p["href"].rsplit("/outputs/", 1)[1],
                    p["hashLength"],
                    p["opener"],
                    "top",
                    p["top"],
                    "window",
                    p["window"],
                    "ora",
                    p["oraAt"],
                    p["bars"],
                    p["ora"],
                    "m9",
                    p["m9At"],
                    p["m9Rows"],
                    p["m9Summary"],
                )


if __name__ == "__main__":
    main()
