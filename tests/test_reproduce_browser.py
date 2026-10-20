"""Y.1: a Reproduce panel prints its shared steps once, and each tab still copies and downloads
the whole script for that tab, byte for byte what main's panel wrote (``panel_golden``).

The page's script joins each tab's step 1 to the shared steps; this test reads what Chromium
copies and downloads, at 1,280 and 375 px in both themes. Skipped where no Chromium is found;
set ``CTX_CHROMIUM``. Each test names the mutation it catches.
"""

import base64
import hashlib
import re

import plotly.graph_objects as go
import pytest
from panel_golden import GOLDEN, WHERE, spec
from test_export_browser import CHROMIUM, _browse

from ctx.figures import panel
from ctx.report.builder import FigureBlock, FigurePage, write

# Every tab clicked, its copy button pressed with the clipboard recorded (the selection where a
# page cannot use it), its download link read; the shared block's state and the page's width.
ACT = """(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  let copied = null;
  Object.defineProperty(navigator, "clipboard", {configurable: true,
    value: {writeText: t => { copied = t; return Promise.resolve(); }}});
  document.querySelectorAll("details.repro").forEach(d => { d.open = true; });
  const tabs = {};
  for (const tab of document.querySelectorAll(".repro [role=tab]")) {
    tab.click(); await wait(50);
    const shown = document.getElementById(tab.getAttribute("aria-controls"));
    copied = null;
    shown.querySelector("button.copy").click(); await wait(150);
    const link = shown.querySelector("a[download]");
    tabs[tab.id.split("-").pop()] = {copied: copied === null ? window.getSelection().toString()
      : copied, href: link.href, name: link.getAttribute("download"),
      visible: link.getBoundingClientRect().width > 0,
      open: Array.from(document.querySelectorAll(".repro-tab")).filter(p => !p.hidden).length};
  }
  return {tabs, shared: document.querySelector(".repro-shared").hidden,
          wide: document.documentElement.scrollWidth, window: window.innerWidth,
          theme: document.documentElement.getAttribute("data-theme")};
})()"""


@pytest.fixture(scope="module")
def page(tmp_path_factory):
    """One figure page per place ctx can come from, each with its Reproduce panel."""
    folder = tmp_path_factory.mktemp("repro")
    pages = {}
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(panel, "_lock", lambda: None)
        for where in WHERE:
            block = FigureBlock(go.Figure(), "No test.", [], "a.csv",
                                reproduce=panel.block(spec(where)))  # fmt: skip
            pages[where] = write(FigurePage("Panel", block), folder / f"{where}.html")
    return pages


def test_the_shared_steps_are_printed_once(page, monkeypatch):
    """Each panel's page holds the shared steps once and each tab's step 1 once, and no tab's
    script as base64. Catches a tab that prints the whole block again, as F2's three did."""
    monkeypatch.setattr(panel, "_lock", lambda: None)
    for where, path in page.items():
        shown = path.read_text()
        assert shown.count("ctx figure desynth_m5_volcano --params params.toml") == 1, where
        assert shown.count("# 2. ") == 1 and 'href="data:text/x-shellscript' not in shown, where
        assert len(re.findall(r'<code id="ctx-repro-0-[a-z]+" data-shared=', shown)) == 3


def test_body_and_script_are_main_s_text(monkeypatch):
    """``panel.body`` and ``panel.script`` of every tab equal main's, by SHA-256. Catches the
    blank line between step 1 and step 2 lost, and the script head moved."""
    monkeypatch.setattr(panel, "_lock", lambda: None)
    for where in WHERE:
        for tab in panel.TABS:
            for kind, written in (("body", panel.body), ("script", panel.script)):
                text = written(spec(where), tab).encode()
                assert hashlib.sha256(text).hexdigest() == GOLDEN[f"{where}/{tab}/{kind}"]


@pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
def test_each_tab_copies_and_downloads_its_whole_script(page, tmp_path, monkeypatch):
    """In Chromium at 1,280 and 375 px, light and dark: each tab's copy is ``panel.body`` and
    its download ``panel.script`` with its file name, the shared block hidden, one tab shown,
    and the page no wider than the window. Catches the join's blank line, the script head left
    out of the download, a download holding another tab's block, and the shared steps shown
    twice with JavaScript."""
    monkeypatch.setattr(panel, "_lock", lambda: None)
    tasks = [{"name": f"{where} {width} {scheme}", "url": path.as_uri(), "width": width,
              "height": 900, "scheme": scheme, "act": ACT}
             for where, path in page.items() for width in (1280, 375)
             for scheme in ("light", "dark")]  # fmt: skip
    results = _browse(tasks, tmp_path)
    for name, seen in results.items():
        assert seen["errors"] == [], name
        acted, where = seen["acted"], name.split()[0]
        assert acted["shared"] is True and acted["wide"] <= acted["window"], name
        assert acted["theme"] == ("black" if name.endswith("dark") else "white"), name
        assert list(acted["tabs"]) == list(panel.TABS), name
        for tab, got in acted["tabs"].items():
            assert got["copied"] == panel.body(spec(where), tab), (name, tab)
            prefix = "data:text/x-shellscript;base64,"
            assert got["href"].startswith(prefix), (name, tab)
            script = base64.b64decode(got["href"][len(prefix) :]).decode()
            assert script == panel.script(spec(where), tab), (name, tab)
            assert got["name"] == f"reproduce_desynth_m5_volcano_desynth_{tab}.sh"
            assert got["visible"] and got["open"] == 1, (name, tab)


def test_without_javascript_every_step_reads_once(page):
    """The page as served, before any script: each tab's step 1, then the shared steps, the
    download links marked for JavaScript. Catches a no-script reader left without steps 2 to
    5, and a download link that does nothing shown to that reader."""
    shown = page["git"].read_text()
    assert shown.index('id="ctx-repro-0-pixi"') < shown.index('id="ctx-repro-0-shared"')
    assert len(re.findall(r'<p class="js-only"><a download=', shown)) == len(panel.TABS)
