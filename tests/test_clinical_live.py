"""M10's Cox screen on the live selection and its live captions, in Chromium and under node.

The screen's table follows a scored list, ordered as the build orders it; the captions of the
score's two figures are refilled for the list and restored at the build list; the screen's
strongest genes scored on the patients that chose them are flagged. Each test names the
mutation it catches.
"""

import json
import re
import shutil
import subprocess

import polars as pl
import pytest
from test_clinical_page import TEMPLATES, clinical_root, controls, outputs  # noqa: F401
from test_export_browser import CHROMIUM, _browse

from ctx.clinical.page import screen_rows
from ctx.report import text

ACT = """(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(3000);
  const form = document.querySelector('form.survival'), fig = form.closest('.fig');
  const table = document.querySelector('table.screen');
  const read = () => ({caption: fig.querySelector('figcaption').textContent,
    rows: Array.from(table.querySelectorAll('tbody tr')).map(r => r.cells[0].textContent),
    title: table.querySelector('caption').textContent,
    notes: fig.querySelector('.survival-notes').textContent});
  const out = {start: read()};
  const cfg = JSON.parse(document.getElementById(form.getAttribute('data-controls')).textContent);
  form.elements.genes.value = cfg.symbols.slice(3, 6).join(' ');
  form.querySelector('[type=submit]').click();
  await wait(1500); out.scored = read();
  form.querySelector('[data-screen]').click();
  await wait(1500); out.screen = read();
  form.querySelector('[data-reset]').click();
  await wait(1500); out.reset = read();
  out.symbols = cfg.symbols.slice(3, 6);
  return out;
})()"""


@pytest.fixture(scope="module")
def ran(clinical_root, tmp_path_factory):  # noqa: F811 - the fixture
    page = outputs(clinical_root, "overview_v1.html")
    return page, _browse([{"name": "m10", "url": page.as_uri(), "width": 1280, "height": 900,
                           "act": ACT}], tmp_path_factory.mktemp("downloads"))  # fmt: skip


@pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
def test_the_screen_and_the_captions_follow_the_scored_list(ran):
    """Mutations caught: the screen left at the build's rows; the caption left at the
    build's numbers; the build's caption not restored; the double-dip flag missing."""
    page, results = ran
    seen = results["m10"]
    assert seen["errors"] == []
    a = seen["acted"]
    words = text.raw("clinical", "live")
    assert sorted(a["scored"]["rows"]) == sorted(a["symbols"])
    assert a["scored"]["title"].startswith(f"The screen's rows for the {len(a['symbols'])} genes")
    assert a["scored"]["caption"] != a["start"]["caption"] and "[[" not in a["scored"]["caption"]
    assert words["double_dip"] in a["screen"]["notes"]
    assert words["double_dip"] not in a["scored"]["notes"]
    assert (
        a["reset"]["rows"] == a["start"]["rows"] and a["reset"]["caption"] == a["start"]["caption"]
    )


def test_the_screen_rows_sort_as_the_table_does(clinical_root):  # noqa: F811 - the fixture
    """ctxLive.screen: a list's rows by adjusted q, the refused last, ties by panel row, and
    the count of its genes with no row. Mutation: ties left in the list's order."""
    html = outputs(clinical_root, "overview_v1.html").read_text()
    cfg = controls(html, "survival")
    rows = [[0, 0.5, ["a"]], [1, None, ["b"]], [2, 0.1, ["c"]], [3, 0.5, ["d"]]]
    script = re.sub(r"\{% (end)?raw %\}", "", (TEMPLATES / "_live.js.j2").read_text())
    harness = (f"const r = ctxLive.screen({json.dumps(rows)}, [3, 1, 9, 0, 2]);"
               "console.log(JSON.stringify(r));")  # fmt: skip
    node = shutil.which("node")
    done = subprocess.run([node, "-e", "var ctxSurv = {};" + script + harness],
                          capture_output=True, text=True, check=True)  # fmt: skip
    out = json.loads(done.stdout)
    assert [r[0] for r in out["rows"]] == [2, 0, 3, 1] and out["absent"] == 1
    assert len(cfg["screen"]) == cfg["genes"]


def _strongest(rows: list, n: int) -> dict:
    """ctxLive.strongest under node."""
    script = re.sub(r"\{% (end)?raw %\}", "", (TEMPLATES / "_live.js.j2").read_text())
    harness = f"console.log(JSON.stringify(ctxLive.strongest({json.dumps(rows)}, {n})));"
    node = shutil.which("node")
    done = subprocess.run([node, "-e", "var ctxSurv = {};" + script + harness],
                          capture_output=True, text=True, check=True)  # fmt: skip
    return json.loads(done.stdout)


def test_a_screen_gene_takes_its_hazard_ratios_side_as_its_sign(clinical_root):  # noqa: F811
    """E.3: the screen's strongest genes entered the list each signed +1, so a protective gene
    subtracted from the harmful ones in the signed mean. Planted: gene B with an alone HR above
    1 and an adjusted HR below it, the one the screen ranks by; gene C at exactly 1. Mutations:
    every sign +1 again in page.screen_rows; the sign taken from the alone HR; the page
    script's signs all +1 again; the sign read from another position of the row."""
    table = pl.DataFrame({
        "gene_name": ["A", "B", "C", "D"], "hr_alone": [1.8, 1.4, 0.9, 2.0],
        "q_alone": [0.01, 0.02, 0.5, 0.01], "hr_adjusted": [1.6, 0.7, 1.0, None],
        "low_adjusted": [1.1, 0.5, 0.8, None], "high_adjusted": [2.2, 0.9, 1.2, None],
        "q_adjusted": [0.03, 0.01, 0.4, None],
    })  # fmt: skip
    rows = screen_rows(table)
    assert [r[3] for r in rows] == [1, -1, 1, None]
    assert _strongest(rows, 3) == {"genes": [1, 0, 2], "signs": [-1, 1, 1]}
    # On the built page every row's sign is its adjusted HR's side of 1, read from the table.
    cfg = controls(outputs(clinical_root, "overview_v1.html").read_text(), "survival")
    built = pl.read_csv(outputs(clinical_root, "cox_screen_v1.csv"), comment_prefix="#")
    signs = [None if q is None or hr is None else (1 if hr >= 1 else -1)
             for hr, q in zip(built["hr_adjusted"], built["q_adjusted"], strict=True)]  # fmt: skip
    assert [r[3] for r in cfg["screen"]] == signs
    assert set(signs) - {None} == {1, -1}, "the made-up root has genes on both sides"
