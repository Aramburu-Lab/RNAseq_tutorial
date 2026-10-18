"""M1 Figure 1's recount (work-order item 5, D-051): ``ctxRecount`` held to its one Python
definition, ``ctx.gather.composition.recount``, under node (D-012); the bars as built are that
definition's; the page recounts a contrast's runs, a colour and a count in Chromium.

Each test names the mutation it catches.
"""

import itertools
import json
import re
from dataclasses import replace
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from test_curation import FIXTURE, SPEC, _gather, _place_sdrf, sdrf_rows
from test_export_browser import CHROMIUM, _browse
from test_page_script import _node

from ctx.dataset import DatasetError, load
from ctx.de.deseq import select_samples
from ctx.gather import composition
from ctx.gather.fetchngs import read

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src" / "ctx" / "report" / "templates" / "_qc.js.j2"
# A contrast on the curation fixture's groups, and a paired one, so a filter, the two levels
# and the paired rule each decide which runs a contrast keeps.
CONTRASTS = """
[[contrasts]]
name = "a_case"
label = "Type A, case against control"
variable = "condition"
reference = "control"
tested = "case"
n = [1, 1]
filters = { cell_type = ["type_a"] }

[[contrasts]]
name = "paired"
variable = "condition"
reference = "control"
tested = "case"
n = [1, 1]
paired = true
"""
HARNESS = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
console.log(JSON.stringify(input.cases.map(([keep, colour, count]) => {
  const bars = ctxRecount.counts(input.rows, keep, colour, count);
  return {bars: bars, smallest: ctxRecount.smallest(bars)};
})));
"""


def _rows(seed: int = 3, n: int = 60) -> pl.DataFrame:
    """Runs with repeated individuals across groups and levels, as ``composition.runs`` gives."""
    rng = np.random.default_rng(seed)
    frame = pl.DataFrame({
        "group": [f"g{i}" for i in rng.integers(0, 6, n)],
        "donor": [f"d{i}" for i in rng.integers(0, 15, n)],
        "batch": [["short", "long"][i] for i in rng.integers(0, 2, n)],
        "sex": [["F", "M", ""][i] for i in rng.integers(0, 3, n)],
    })  # fmt: skip
    return frame.sort("group", maintain_order=True)


def test_the_page_twin_recounts_as_python_does():
    """Every count, colour and kept subset, the empty one too, and the smallest group.

    Mutations caught: individuals counted per run; a group's total summed over its levels;
    the tie broken by order instead of name."""
    rows = _rows()
    rng = np.random.default_rng(5)
    keeps = [None, sorted(rng.choice(rows.height, 25, replace=False).tolist()), [], [0]]
    cases = list(itertools.product(keeps, ("batch", "sex"), composition.COUNTS))
    got = _node({"rows": {c: rows[c].to_list() for c in rows.columns}, "cases": cases},
                HARNESS, SCRIPT)  # fmt: skip
    for (keep, colour, count), seen in zip(cases, got, strict=True):
        kept = rows if keep is None else rows[keep]
        bars = composition.recount(kept, colour, count)
        assert seen["bars"] == bars.to_dicts(), (keep, colour, count)
        expected = list(composition.smallest(bars)) if bars.height else [0, "", 0]
        assert seen["smallest"] == expected


def test_a_tie_for_the_smallest_group_goes_to_the_first_name():
    """Both twins, on groups whose order of appearance is not their names' order. Mutation
    caught: the page's tie broken by order."""
    rows = pl.DataFrame({"group": ["b", "a", "c", "c"], "donor": ["d1", "d2", "d3", "d4"],
                         "batch": "x"})  # fmt: skip
    bars = composition.recount(rows, "batch", "runs")
    assert composition.smallest(bars) == (3, "a", 1)
    got = _node({"rows": {c: rows[c].to_list() for c in rows.columns},
                 "cases": [[None, "batch", "runs"]]}, HARNESS, SCRIPT)  # fmt: skip
    assert got[0]["smallest"] == [3, "a", 1]


@pytest.fixture(scope="module")
def page(tmp_path_factory):
    """ctx gather on the curation fixture, its dataset file given two contrasts."""
    root = tmp_path_factory.mktemp("recount")
    spec = root / "synthetic.toml"
    spec.write_text(SPEC.read_text() + CONTRASTS)
    _place_sdrf(root, sdrf_rows(read(FIXTURE).runs))
    done = _gather(root, spec)
    assert done.returncode == 0, done.stderr[-3000:]
    out = root / "outputs" / "01_gather"
    (html,) = out.glob("*_synthetic_overview_v1.html")
    (sheet,) = out.glob("*_synthetic_samplesheet_v1.parquet")
    return html, pl.read_parquet(sheet), load(spec)


def _cfg(html: str) -> dict:
    figure = next(f for f in html.split('<figure class="fig"') if 'class="calls recount"' in f)
    return json.loads(re.search(r'id="ctx-controls-\d+">(.*?)</script>', figure, re.S).group(1))


def test_a_contrast_keeps_the_runs_m5_compares(page):
    """Mutation caught: a contrast's runs taken from its filters alone, the two levels and the
    pairing ignored."""
    html, sheet, spec = page
    cfg = _cfg(html.read_text())
    ordered = sheet.sort("cell_type", "condition", "run_accession", maintain_order=True)
    for contrast, held in zip(spec.contrasts, cfg["contrasts"], strict=True):
        used, _ = select_samples(sheet, contrast, spec.order)
        assert held["runs"] == [i for i, s in enumerate(ordered["sample"])
                                if s in set(used["sample"])]  # fmt: skip
    assert "Paired" in cfg["contrasts"][1]["line"]
    assert list(cfg["colours"]) == ["batch", "time_point"]


def test_the_bars_as_built_are_the_recount_of_every_run(page):
    """Mutation caught: the build's individuals counted as runs."""
    html, sheet, _ = page
    (csv,) = html.parent.glob("*_synthetic_composition_v1.csv")
    drawn = pl.read_csv(csv)
    rows = composition.runs(sheet, ["batch"])
    again = composition.recount(rows, "batch", "individuals")
    assert dict(zip(drawn["group"], drawn["n_individuals"], strict=True)) == dict(
        zip(again["group"], again["n"], strict=True))  # fmt: skip


ACT = """(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(1500);
  const form = document.querySelector('form.recount'), fig = form.closest('.fig');
  fig.setAttribute('data-probe', '1');
  const d = fig.querySelector('.js-plotly-plot');
  const read = () => ({traces: d.data.map(t => [t.name, Array.from(t.y), Array.from(t.x)]),
    title: d.layout.xaxis.title.text, lead: fig.querySelector('.fig-lead').textContent,
    line: fig.querySelector('.recount-summary').textContent,
    moved: !!fig.querySelector('a[data-export][role=button]')});
  const out = {start: read()};
  form.elements.contrast.value = '0'; form.elements.colour.value = 'time_point';
  form.elements.count.value = 'individuals';
  form.dispatchEvent(new Event('change', {bubbles: true}));
  await wait(600);
  /* A resize redraws every plot, and the title fitter with it, after the recount. */
  window.dispatchEvent(new Event('resize'));
  await wait(1500);
  out.after = read();
  return out;
})()"""


@pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
def test_the_page_draws_the_recount_and_exports_it(page, tmp_path):
    """In Chromium at 375 px: a contrast, a colour and a count give Python's bars, the lead
    names the smallest group in individuals, the files are drawn from the view. Mutation
    caught: the x title left as runs."""
    html, sheet, spec = page
    task = {"name": "r", "url": html.as_uri(), "width": 375, "height": 800, "act": ACT,
            "figure": "[data-probe]"}  # fmt: skip
    seen = _browse([task], tmp_path)["r"]
    assert seen["errors"] == []
    start, after = seen["acted"]["start"], seen["acted"]["after"]
    assert not start["moved"] and after["moved"] and after["title"] == "Individuals"
    rows = composition.runs(sheet, composition.colours(sheet, spec))
    keep = _cfg(html.read_text())["contrasts"][0]["runs"]
    bars = composition.recount(rows[keep], "time_point", "individuals")
    drawn = {(name, g): x for name, ys, xs in after["traces"] for g, x in zip(ys, xs, strict=True)}
    assert drawn == {(lv, g): n for g, lv, n in bars.select("group", "level", "n").rows()}
    groups, name, n = composition.smallest(bars)
    assert f"the smallest, {name}, holds {n} individual" in after["lead"]
    assert after["line"].startswith("The ") and after["line"].endswith(
        "Bars count individuals, coloured by time_point.")  # fmt: skip
    svg = (tmp_path / seen["files"]["svg"]["guid"]).read_text()
    words = " ".join(re.findall(r"<text[^>]*>([^<]*)</text>", svg))
    assert "<image" not in svg and "coloured by time_point" in words


def test_colours_are_the_ones_the_dataset_file_names(page):
    """M1 offers batch, then only the colours `[dataset] colours` names: Step-1's tissue and time
    point nest in cell type and condition (Rytis 2026-10-18). Mutations caught: every curated
    column offered again; an undeclared name passed through."""
    _, sheet, spec = page
    assert composition.colours(sheet, spec) == ["batch", "time_point"]
    unnamed = replace(spec, info={k: v for k, v in spec.info.items() if k != "colours"})
    assert composition.colours(sheet, unnamed) == ["batch"]
    wrong = replace(spec, info={**spec.info, "colours": ["tissue_of_nowhere"]})
    with pytest.raises(DatasetError, match="tissue_of_nowhere"):
        composition.colours(sheet, wrong)
