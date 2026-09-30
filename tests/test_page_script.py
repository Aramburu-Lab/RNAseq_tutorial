"""The page script's live call rule, run under node, against ctx.de.stats.classify (D-012).

The volcano page calls genes again in the browser for any cutoff the reader picks. That rule
is JavaScript (templates/_calls.js.j2); this test runs the same file and compares its counts
and cutoff line with the Python definition over a grid of every basis, cutoff and threshold.
"""

import dataclasses
import json
import math
import shutil
import subprocess
from collections import Counter
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from synth import demo_linked, design_counts, filtered_counts
from test_viz import filtered_bars

from ctx.dataset import Contrast
from ctx.de.deseq import fit
from ctx.de.run import _volcano
from ctx.de.stats import ADJUSTMENTS, adjust, classify, p_line
from ctx.report import text
from ctx.report.builder import _defaults
from ctx.viz.bar import significant
from ctx.viz.save import figure_block
from ctx.viz.scatter import CUTOFF_PRESETS, MIN_CUTOFF, _js_number

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src" / "ctx" / "report" / "templates" / "_calls.js.j2"
# The presets as the page's controls carry them, from their one definition.
CFG = {"gate": "padj", "lfc": "log2FoldChange", "p": "pvalue", "alpha": 0.05,
       "min_cutoff": MIN_CUTOFF, "max_cutoff": max(CUTOFF_PRESETS), "n_tested": 400,
       "presets": [_js_number(c) for c in CUTOFF_PRESETS]}  # fmt: skip
# The tested genes _rows makes; the untested one is not counted.
BASES = (*ADJUSTMENTS, "pvalue")
CUTOFFS = (0.1, 0.05, 0.01, 0.001)
THRESHOLDS = (0.0, math.log2(1.2), 1.0, 2.0)


def _rows() -> pl.DataFrame:
    """Random tested genes, plus the boundary cases the rule must get exactly right."""
    rng = np.random.default_rng(7)
    n = CFG["n_tested"]
    p = np.concatenate([rng.uniform(0, 1e-4, 40), rng.uniform(0, 1, n - 40)])
    lfc = rng.normal(0, 1.5, n)
    lfc[:6] = [1.0, 1.0 - 1e-9, -1.0, -1.0 + 1e-9, 2.0, 0.0]  # at and just inside |LFC| = 1
    p[:6] = 1e-6
    p[6] = 0.0  # pydeseq2 returns 0 for the strongest effects
    frame = pl.DataFrame({"gene": [f"g{i}" for i in range(n)], "log2FoldChange": lfc, "pvalue": p})
    frame = frame.with_columns(**{column: pl.Series(adjust(p, column)) for column in ADJUSTMENTS})
    # A padj exactly at the default cutoff: strictly below is called, so this is not.
    frame = frame.with_columns(
        padj=pl.when(pl.col("gene") == "g10").then(pl.lit(0.05)).otherwise(pl.col("padj"))
    )
    untested = pl.DataFrame(
        {"gene": ["untested"], "log2FoldChange": [5.0], "pvalue": [1e-9]}
    ).with_columns(**{column: pl.lit(None, pl.Float64) for column in ADJUSTMENTS})
    return pl.concat([frame, untested])


def _python(rows: pl.DataFrame, basis: str, alpha: float, t: float) -> dict[str, float | None]:
    """The oracle: classify on the chosen p, over tested genes only; p_line for the line."""
    chosen = pl.when(pl.col("padj").is_not_null()).then(pl.col(basis))
    calls = rows.select(classify(chosen, pl.col("log2FoldChange"), alpha, t)).to_series()
    passing = rows.filter(pl.col("padj").is_not_null() & (pl.col(basis) < alpha)).height
    return {
        "up": int((calls == "up").sum()),
        "down": int((calls == "down").sum()),
        "line": p_line(basis, alpha, CFG["n_tested"], passing),
    }


def _node(payload: dict, harness: str | None = None, script: Path = SCRIPT) -> dict:
    """Run a template's logic under node on one payload; returns what the harness prints."""
    node = shutil.which("node")
    assert node, "nodejs is pinned in the ctx-dev environment; run tests with pixi"
    rule = script.read_text().replace("{% raw %}", "").replace("{% endraw %}", "")
    harness = (
        harness
        or """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
const out = {
  grid: input.grid.map(([basis, alpha, t]) => {
    const found = ctxCalls.call(input.rows, input.cfg, t, basis, alpha);
    return { up: found.parts.up.length, down: found.parts.down.length, line: found.line };
  }),
  thresholds: input.thresholds.map(([v, unit]) => ctxCalls.threshold(v, unit)),
  cutoffs: input.cutoffs.map((v) => ctxCalls.cutoff(v, input.cfg)),
  walks: [-1, 1].map((dir) => {
    let v = String(input.cfg.alpha); const seen = [];
    for (let i = 0; i < 6; i++) {
      v = String(ctxCalls.step(v, dir, input.cfg)); seen.push(Number(v));
    }
    return seen;
  }),
  between: [ctxCalls.step("0.03", -1, input.cfg), ctxCalls.step("0,03", 1, input.cfg)],
};
console.log(JSON.stringify(out));
"""
    )
    done = subprocess.run(
        [node, "-e", rule + harness], input=json.dumps(payload), capture_output=True, text=True,
        check=True,
    )  # fmt: skip
    return json.loads(done.stdout)


@pytest.fixture(scope="module")
def ran():
    rows = _rows()
    grid = [(b, a, t) for b in BASES for a in CUTOFFS for t in THRESHOLDS]
    payload = {
        "rows": rows.to_dicts(),
        "cfg": CFG,
        "grid": grid,
        "thresholds": [["20", "pct"], ["1", "log2"], ["-3", "log2"], ["", "log2"], ["100", "pct"]],
        "cutoffs": ["0", "", "0.5", "0.00001", "abc", "0.01", "0,01"],
    }
    return rows, grid, _node(payload)


def test_counts_and_line_match_classify(ran):
    rows, grid, out = ran
    assert len(out["grid"]) == len(grid) == 80
    for (basis, alpha, t), page in zip(grid, out["grid"], strict=True):
        # The BY constant is summed pairwise by numpy and in order by the page: last bits differ.
        want = _python(rows, basis, alpha, t)
        assert page == want | {"line": pytest.approx(want["line"], rel=1e-12)}, (basis, alpha, t)


def test_the_grid_exercises_the_boundaries(ran):
    # Without these the comparison above could pass while the edges were wrong.
    rows, _, out = ran
    at_default = _python(rows, "padj", 0.05, 1.0)
    assert at_default["up"] > 0 and at_default["down"] > 0
    assert _python(rows, "padj_bonferroni", 0.001, 2.0) != _python(rows, "padj", 0.1, 0.0)
    raw = _python(rows, "pvalue", 0.1, 0.0)
    everything = rows.filter(pl.col("pvalue") < 0.1).height
    assert raw["up"] + raw["down"] == everything - 1  # the untested gene is never called


def test_threshold_and_cutoff_parsing(ran):
    *_, out = ran
    assert out["thresholds"] == pytest.approx([math.log2(1.2), 1.0, 0.0, 0.0, 1.0])
    # 0 and values below the minimum are raised to it, empty and junk give the default,
    # values above the maximum are lowered to it.
    assert out["cutoffs"] == [0.0001, 0.05, 0.1, 0.0001, 0.05, 0.01, 0.01]


def test_arrows_walk_every_preset(ran):
    """The arrows reach every preset from the default and hold at both ends.

    A number box's own arrows went from 0.05 to 0.0001 and then nowhere: they step by 1 from
    the minimum, which the range 0.0001 to 0.1 cannot hold.
    """
    *_, out = ran
    assert out["walks"] == [[0.01, 0.001, 0.0001, 0.0001, 0.0001, 0.0001],
                            [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]]  # fmt: skip
    assert out["between"] == [0.01, 0.05]


LOAD = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
const cfg = input.cfg;
const found = ctxCalls.call(input.rows, cfg, cfg.lfc_min, cfg.gate, cfg.alpha);
console.log(JSON.stringify({ line: found.line, up: found.parts.up.length,
                             down: found.parts.down.length }));
"""


def test_static_line_is_the_line_the_page_draws_on_load():
    # The PDF, the SVG and the no-JavaScript view show Calls.p_line; the script redraws the
    # line on load from the embedded rows and cfg.n_tested. A producer counting m or k on
    # other genes than the page does would ship two lines, and only this test would see it.
    # Filtered genes, so the tested set is not every gene and m has one right value.
    synthetic = design_counts(condition_lfc=1.0)
    counts = filtered_counts(synthetic)
    result = fit(counts, synthetic.samples,
                 Contrast("c", "condition", "control", "case", (6, 6)), {}, threads=1)  # fmt: skip
    assert result.results["padj"].null_count() > 0
    data = _volcano(result, counts.select("gene_id", gene_name=pl.col("gene_id")))
    block = figure_block(data, "x")
    rows = [{k: None if isinstance(v, float) and not math.isfinite(v) else v
             for k, v in r.items()} for r in block.rows]  # fmt: skip
    out = _node({"rows": rows, "cfg": block.controls}, LOAD)
    calls = data.frame["call"].value_counts()
    assert out == {
        "line": pytest.approx(data.calls.p_line, rel=1e-12),
        "up": dict(calls.rows())["up"],
        "down": dict(calls.rows())["down"],
    }


SWITCH = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
const up = (t) => ctxCalls.call(input.rows, input.cfg, t, "padj", 0.05).parts.up.map((r) => r.gene);
const state = { exact: null, shown: null, unit: "log2" };
const steps = [];
let t = ctxCalls.applied(state, "1", "log2");
steps.push({ step: "log2 1", t, up: up(t) });
// The select moved to % without a change event, as when a browser restores the form.
t = ctxCalls.applied(state, "34", "pct");
steps.push({ step: "34 % typed", t, up: up(t) });
let shown = ctxCalls.switchUnit(state, "34", "log2");
t = ctxCalls.applied(state, shown, "log2");
steps.push({ step: "switched to log2", shown, t, up: up(t) });
shown = ctxCalls.switchUnit(state, shown, "pct");
t = ctxCalls.applied(state, shown, "pct");
steps.push({ step: "switched back to %", shown, t, up: up(t) });
shown = ctxCalls.switchUnit(state, shown, "log2");
t = ctxCalls.applied(state, "0.420", "log2");
steps.push({ step: "typed 0.420", shown, t, up: up(t) });
console.log(JSON.stringify(steps));
"""


def test_unit_switch_keeps_the_exact_threshold():
    # One gene between the rounded and the exact threshold: log2 1.34 = 0.4222 shows as 0.42,
    # and a gene at 0.4211 is called only if the rounded number were applied.
    exact = math.log2(1.34)
    rows = [
        {"gene": "edge", "log2FoldChange": 0.4211, "pvalue": 1e-4, "padj": 0.001},
        {"gene": "strong", "log2FoldChange": 2.0, "pvalue": 1e-6, "padj": 0.0001},
    ]
    steps = _node({"rows": rows, "cfg": CFG}, SWITCH)
    got = {s["step"]: s for s in steps}
    assert got["log2 1"]["up"] == ["strong"]
    assert got["34 % typed"]["t"] == pytest.approx(exact)
    assert got["34 % typed"]["up"] == ["strong"]
    assert got["switched to log2"]["shown"] == "0.42"
    assert got["switched to log2"]["t"] == pytest.approx(exact)
    assert got["switched to log2"]["up"] == ["strong"], "the rounded 0.42 was applied"
    assert got["switched back to %"]["shown"] == "34.0"
    assert got["switched back to %"]["t"] == pytest.approx(exact)
    # A number the reader types replaces the kept one.
    assert got["typed 0.420"]["t"] == 0.42
    assert got["typed 0.420"]["up"] == ["edge", "strong"]


BRUSH = ROOT / "src" / "ctx" / "report" / "templates" / "_brush.js.j2"
PAINT = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
const cfg = input.cfg, rows = input.rows;
let painted = ctxBrush.start(rows.length);
const out = [];
for (const [indices, event, k] of input.steps) {
  painted = ctxBrush.paint(painted, event ? ctxBrush.indices(event) : indices, k);
  out.push({ summary: ctxBrush.summary(painted, rows, cfg),
             markers: ctxBrush.markers(painted, cfg), codes: ctxBrush.codes(painted) });
}
painted = ctxBrush.clear(painted);
out.push({ summary: ctxBrush.summary(painted, rows, cfg), codes: ctxBrush.codes(painted) });
console.log(JSON.stringify(out));
"""


def _brush_oracle(painted: list[int], rows: list[dict], cfg: dict) -> str:
    """The live line in Python: per paint in palette order, n and the detail, most first."""
    lines = []
    for k in range(len(cfg["colours"])):
        members = [str(rows[i][cfg["detail"]]) for i, p in enumerate(painted) if p == k]
        if members:
            counts = sorted(Counter(members).items(), key=lambda kv: (-kv[1], kv[0]))
            detail = ", ".join(f"{level} {n:,}" for level, n in counts)
            lines.append(cfg["summary"].format(n=f"{len(members):,}", colour=cfg["names"][k],
                                               detail=detail))  # fmt: skip
    return "; ".join(lines) or cfg["empty"]


def test_brush_paints_both_views_and_counts_as_python_does():
    data = demo_linked()
    block = figure_block(data, "x")
    rows, cfg = list(block.rows), dict(block.controls)
    selected = {"points": [{"curveNumber": 0, "pointIndex": i, "pointNumber": i}
                           for i in (10, 11, 12)]}  # fmt: skip
    tapped = {"points": [{"curveNumber": 1, "pointNumber": i} for i in (20, 21, 3)]}
    steps = [[[0, 1, 2, 3], None, 0], [[2, 3, 4], None, 5], [[], selected, 2], [[], tapped, 0]]
    out = _node({"rows": rows, "cfg": cfg, "steps": steps}, PAINT, BRUSH)
    painted = [-1] * len(rows)
    for (indices, event, k), page in zip(steps, out, strict=False):
        picked = (
            [p.get("pointNumber", p.get("pointIndex")) for p in event["points"]]
            if event
            else indices
        )
        for i in picked:
            painted[i] = k
        assert page["summary"] == _brush_oracle(painted, rows, cfg)
        assert page["codes"] == [p + 1 for p in painted]
        assert page["markers"] == [cfg["reference"] if p < 0 else cfg["colours"][p]
                                   for p in painted]  # fmt: skip
    # A repaint moves a sample: 2 and 3 went from orange to vermillion, then 3 back.
    assert painted[:5] == [0, 0, 5, 0, 5]
    assert out[-1] == {"summary": cfg["empty"], "codes": [0] * len(rows)}
    assert "; " in out[-2]["summary"]


FILTER = ROOT / "src" / "ctx" / "report" / "templates" / "_filter.js.j2"
PICK = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
const cfg = input.cfg;
console.log(JSON.stringify({
  views: input.choices.map((chosen) => {
    const view = ctxFilter.view(input.rows, cfg, chosen);
    return { view, summary: ctxFilter.summary(view, cfg, input.template, input.none) };
  }),
  numbers: input.numbers.map((v) => ctxFilter.number(v)),
}));
"""


def _filter_oracle(data, chosen: dict[str, str]) -> dict:
    """What the page must draw for one choice, from BarData itself: its rows, its lengths()."""
    frame = data.frame
    for column, level in chosen.items():
        if level:
            frame = frame.filter(pl.col(column).cast(pl.String) == level)
    if frame.is_empty():
        return {"categories": [], "bars": 0, "median": None, "traces": [[] for _ in data.colours]}
    lengths = dataclasses.replace(data, frame=frame).lengths()
    parts = [frame.filter(pl.col(data.series) == level) for level in data.colours]
    return {
        "categories": frame[data.category].unique(maintain_order=True).cast(pl.String).to_list(),
        "bars": len(lengths),
        "median": float(np.median(lengths)),
        "traces": [p[data.category].cast(pl.String).to_list() for p in parts],
        "values": [p[data.value].to_list() for p in parts],
        "hover": [p.select(data.hover).rows() for p in parts],
    }


@pytest.mark.parametrize(("kind", "horizontal"),
                         [("runs", False), ("groups", False), ("grouped", True)])  # fmt: skip
def test_filter_draws_what_bardata_defines(kind, horizontal):
    """ctxFilter against BarData: the rows kept, the bars per trace, the categories, the
    median over lengths() (a total per category when stacked) and the live line."""
    data = dataclasses.replace(filtered_bars(kind), horizontal=horizontal)
    cfg = data.controls()
    levels = {c: ["", *v] for c, v in cfg["filters"].items()}
    choices = [{"cell_type": t, "condition": c}
               for t in levels["cell_type"] for c in levels["condition"]]  # fmt: skip
    choices.append({"cell_type": "NDN", "condition": "nowhere"})  # a choice no bar matches
    template, none = text.raw("ui", "filter.summary"), text.word("ui", "filter.none")
    numbers = [0.0, 0.5, 2.675, 7 / 3, 26.116325, 45.25, 99.96, 999.5, 1234.5, 0.0001234,
               1_000_000.5, 123_456_789.0]  # fmt: skip
    payload = {"rows": data.page_table().to_dicts(), "cfg": cfg, "choices": choices,
               "template": template, "none": none, "numbers": numbers}  # fmt: skip
    out = _node(payload, PICK, FILTER)
    assert out["numbers"] == [significant(v) for v in numbers]
    along, at = ("x", "y") if horizontal else ("y", "x")
    for chosen, page in zip(choices, out["views"], strict=True):
        want, view = _filter_oracle(data, chosen), page["view"]
        assert view["categories"] == want["categories"], chosen
        assert (view["bars"], view["median"]) == (want["bars"], pytest.approx(want["median"]))
        if not want["bars"]:
            assert page["summary"] == none
            continue
        assert view[at] == want["traces"] and view[along] == want["values"], chosen
        assert view["customdata"] == [[list(r) for r in h] for h in want["hover"]]
        filled = {"shown": f"{want['bars']:,}", "total": f"{cfg['total']:,}",
                  "median": f"{significant(want['median'])} {cfg['unit']}"}  # fmt: skip
        assert page["summary"] == template.format(**filled)
    # "All" in every select is what the page shows before any change, and without JavaScript.
    assert out["views"][0]["summary"] == _defaults(cfg)["summary"]
    assert any(not p["view"]["bars"] for p in out["views"])
    assert len({p["summary"] for p in out["views"]}) > len(choices) // 2
