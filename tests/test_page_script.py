"""The page script's live call rule, run under node, against ctx.de.stats.classify (D-012).

The volcano page calls genes again in the browser for any cutoff the reader picks. That rule
is JavaScript (templates/_calls.js.j2); this test runs the same file and compares its counts
and cutoff line with the Python definition over a grid of every basis, cutoff and threshold.
"""

import json
import math
import shutil
import subprocess
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from ctx.de.stats import ADJUSTMENTS, adjust, classify, p_line

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src" / "ctx" / "report" / "templates" / "_calls.js.j2"
CFG = {"gate": "padj", "lfc": "log2FoldChange", "p": "pvalue", "alpha": 0.05,
       "min_cutoff": 0.0001, "max_cutoff": 0.1, "n_tested": 400}  # fmt: skip
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


def _node(payload: dict, harness: str | None = None) -> dict:
    """Run the template's rule under node on one payload; returns what the harness prints."""
    node = shutil.which("node")
    assert node, "nodejs is pinned in the ctx-dev environment; run tests with pixi"
    rule = SCRIPT.read_text().replace("{% raw %}", "").replace("{% endraw %}", "")
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
        "cutoffs": ["0", "", "0.5", "0.00001", "abc", "0.01"],
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
    assert out["cutoffs"] == [0.0001, 0.05, 0.1, 0.0001, 0.05, 0.01]


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
