"""The page script's live call rule, run under node, against ctx.de.stats.classify (D-012).

The volcano page calls genes again in the browser for any cutoff the reader picks. That rule
is JavaScript (templates/_calls.js.j2); this test runs the same file and compares its counts
and cutoff line with the Python definition over a grid of every basis, cutoff and threshold.
The bar filter, its lead's numbers and the hover rule are held to their Python definitions the
same way.
"""

import dataclasses
import itertools
import json
import math
import re
import shutil
import subprocess
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import polars as pl
import pytest
from synth import demo_linked, demo_volcano, design_counts, filtered_counts
from test_viz import filtered_bars

from ctx.dataset import Contrast
from ctx.de.deseq import fit
from ctx.de.run import volcano
from ctx.de.stats import ADJUSTMENTS, adjust, classify, p_line
from ctx.report import text
from ctx.report.builder import FigureText
from ctx.report.figures import defaults as _defaults
from ctx.report.figures import finite
from ctx.viz.axes import hover_rows, hover_value
from ctx.viz.bar import significant
from ctx.viz.save import figure_block
from ctx.viz.scatter import CUTOFF_PRESETS, MIN_CUTOFF, _js_number
from ctx.viz.style import THEMES

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


def _node(payload: dict, harness: str | None = None, script: Path | tuple[Path, ...] = SCRIPT
          ) -> dict:  # fmt: skip
    """Run templates' logic under node on one payload; returns what the harness prints.

    Several templates run in the order given, as the page includes them.
    """
    node = shutil.which("node")
    assert node, "nodejs is pinned in the ctx-dev environment; run tests with pixi"
    scripts = script if isinstance(script, tuple) else (script,)
    rule = "".join(s.read_text() for s in scripts).replace("{% raw %}", "").replace(
        "{% endraw %}", "")  # fmt: skip
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
    data = volcano(result, counts.select("gene_id", gene_name=pl.col("gene_id")))
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


HOVER = ROOT / "src" / "ctx" / "report" / "templates" / "_hover.js.j2"
# The filter as the page includes it: the hover rule first.
FILTER = (HOVER, ROOT / "src" / "ctx" / "report" / "templates" / "_filter.js.j2")
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
        "hover": [hover_rows(p, data.hover) for p in parts],
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


def test_a_filter_line_without_a_median_and_its_selects_in_words():
    """Bars stacked to 100 %: the live line counts the bars shown and prints no median, before
    and after a choice, and each select is labelled in words, its column kept for the script.

    Mutations caught: ``_defaults`` ignoring summary_median (a median of 100 % printed); the
    select labelled by its column name ("cell_type")."""
    data = dataclasses.replace(filtered_bars("groups"), summary_median=False,
                               labels={"cell_type": "Cell type"})  # fmt: skip
    cfg = data.controls()
    shown = _defaults(cfg)
    assert shown["summary"] == f"{cfg['total']:,} of {cfg['total']:,} bars shown."
    assert "median" not in shown["template"]
    assert ("cell_type", cfg["filters"]["cell_type"], "Cell type") in shown["filters"]
    assert ("condition", cfg["filters"]["condition"], "condition") in shown["filters"]
    choices = [{"cell_type": "", "condition": ""}, {"cell_type": "NDN", "condition": ""}]
    payload = {"rows": data.page_table().to_dicts(), "cfg": cfg, "choices": choices,
               "template": shown["template"], "none": shown["none"], "numbers": []}  # fmt: skip
    out = _node(payload, PICK, FILTER)
    assert out["views"][0]["summary"] == shown["summary"]
    assert re.fullmatch(r"\d+ of \d+ bars shown\.", out["views"][1]["summary"])


HOVERS = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
console.log(JSON.stringify({
  values: input.values.map(ctxHover.value), nan: ctxHover.value(NaN),
  views: input.choices.map((chosen) => ctxFilter.view(input.rows, input.cfg, chosen).customdata),
  aggregate: ctxFilter.aggregate(input.rows, input.rows.map((_, i) => i), input.aggregates),
  none: ctxFilter.aggregate(input.rows, [0], input.aggregates),
}));
"""


def test_a_missing_hover_value_prints_nothing_as_python_writes_it():
    """Package C: under a number format plotly prints a missing value as 0, so Python writes it
    as "" (``axes.hover_value``). The page builds customdata again from its rows, where a
    missing number is null, whenever a filter or a live call moves; it follows the same rule.
    A lead's numbers pass a missing value by, in both languages, and with none left give no
    number. Mutation: the page's customdata of raw values (null printed as 0 again)."""
    from ctx.viz import aggregates  # noqa: PLC0415

    data = filtered_bars()
    fdr = [None if i % 4 == 0 else i / 100 for i in range(data.frame.height)]
    data = dataclasses.replace(
        data, frame=data.frame.with_columns(fdr=pl.Series(fdr, dtype=pl.Float64)),
        hover=("cell_type", "fdr"),
        columns={**data.columns, "fdr": ("none", "An FDR.")})  # fmt: skip
    values = [None, 0, 0.0, -1.5, 12, "", "LDN"]
    choices = [{"cell_type": "", "condition": ""}, {"cell_type": "LDN", "condition": ""}]
    rows = [finite(r) for r in data.page_table().to_dicts()]
    specs = (aggregates.Aggregate("m", "median", "fdr"),
             aggregates.Aggregate("low", "min", "fdr", row=("bar", "bar")))  # fmt: skip
    out = _node({"values": values, "rows": rows, "cfg": data.controls(), "choices": choices,
                 "aggregates": [dataclasses.asdict(s) for s in specs]}, HOVERS, FILTER)  # fmt: skip
    assert out["values"] == [hover_value(v) for v in values] and out["nan"] == hover_value(math.nan)
    for chosen, page in zip(choices, out["views"], strict=True):
        assert page == [[list(r) for r in h] for h in _filter_oracle(data, chosen)["hover"]]
    assert any("" in r for trace in out["views"][0] for r in trace)
    assert out["aggregate"] == aggregates.fill(rows, specs) == {"m": "0.200", "low": "0.0100",
                                                                "bar": "r01"}  # fmt: skip
    assert out["none"] is None and aggregates.fill(rows[:1], specs) is None


def _depth(median: int, n: int, seed: int):
    """M1's depth bars and lead over ``n`` runs (odd) whose median is ``median`` reads.

    Monocyte runs are steady state only, so one filter level shows no run.
    """
    from ctx.gather import harmonise  # noqa: PLC0415

    rng = np.random.default_rng(seed)
    half = n // 2
    reads = [*rng.integers(median // 5, median, half), median,
             *rng.integers(median + 1, 3 * median, half)]  # fmt: skip
    types = [("NDN", "LDN", "Monocyte")[i % 3] for i in rng.permutation(n)]
    sheet = pl.DataFrame({
        "run_accession": [f"XRR{i:04d}" for i in range(n)], "read_count": reads,
        "batch": ["short" if i % 2 else "long" for i in range(n)],
        "read_length": [76.0 if i % 2 else 101.0 for i in range(n)], "cell_type": types,
        "condition": ["steady_state" if t == "Monocyte" or i % 2 else "gcsf"
                      for i, t in enumerate(types)],
    })  # fmt: skip
    bands, _ = harmonise._depth_standard()
    data = harmonise._library_sizes(
        sheet, SimpleNamespace(batch_bins={"short": (76,), "long": (101,)}), bands
    )
    return data, harmonise.library_lead(bands[0])


def _samples() -> pl.DataFrame:
    """Ten samples over three cell types and two conditions, for M3's bars; Mono is b only,
    so one filter level shows no sample."""
    return pl.DataFrame({"sample": [f"S{i}" for i in range(10)],
                         "cell_type": ["NDN", "NDN", "LDN", "LDN", "LDN", "Mono", "Mono", "NDN",
                                       "LDN", "Mono"],
                         "condition": ["a", "b", "a", "b", "a", "b", "b", "b", "a",
                                       "b"]})  # fmt: skip


def _mapping():
    """M3's mapping bars and lead: S3 and S6 tie at the lowest uniquely mapped share, so a tie
    keeps the first row."""
    from ctx.process import Aligner, blocks, figures  # noqa: PLC0415
    from ctx.process import tables as tab  # noqa: PLC0415

    unique = [900, 850, 700, 310, 820, 760, 310, 905, 640, 500]
    reads = {f"S{i}": {"unique": u, "multi": 40 + i, "unmapped": 60 - i, "input": u + 100}
             for i, u in enumerate(unique)}  # fmt: skip
    reads["S6"] = dict(reads["S3"])
    frame = tab.outcomes_long(reads, Aligner.HISAT2).join(
        _samples(), left_on="label", right_on="sample", maintain_order="left")  # fmt: skip
    data = figures.outcome_bars(frame, Aligner.HISAT2, "mapping", ("cell_type", "condition"),
                                "Made up; no test.")  # fmt: skip
    return data, blocks.mapping_lead()


def _strand():
    """M3's strand bars and lead: shares below 0.1, near 0.5 and near 1, and the calls tied
    two to two within LDN."""
    from ctx.process import blocks, figures  # noqa: PLC0415

    shares = [0.0213, 0.0456, 0.4921, 0.0312, 0.5034, 0.9611, 0.0278, 0.4877, 0.0453, 0.0391]
    calls = ["reverse" if s < 0.1 else "forward" if s > 0.9 else "unstranded" for s in shares]
    rows = [{"sample": f"S{i}", "setting": s, "reads": 1000 + 10 * i, "percent": 50.0 + i,
             "no_feature": 100, "ambiguous": 10, "forward_share": share, "call": call}
            for i, (share, call) in enumerate(zip(shares, calls, strict=True))
            for s in ("s0", "s1", "s2")]  # fmt: skip
    frame = pl.DataFrame(rows).join(_samples(), on="sample", maintain_order="left")
    data = figures.strand_bars(frame, ("cell_type", "condition"), "Made up; no test.")
    return data, blocks.strand_lead()


# The made-up M1 page printed its median as "0.5 million" in the lead and "0.494 million
# reads" in the live line (kit v6); Step-1's is 26.1 million (work order item 2, D1 and D2).
LEADS = {"depth_made_up": lambda: _depth(494_123, 33, 3),
         "depth_step1": lambda: _depth(26_123_456, 41, 5), "mapping": _mapping,
         "strand": _strand}  # fmt: skip
LEAD = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
const cfg = input.cfg, lead = cfg.lead;
console.log(JSON.stringify(input.choices.map((chosen) => {
  const view = ctxFilter.view(input.rows, cfg, chosen);
  const values = view.bars ? ctxFilter.aggregate(input.rows, view.shown, lead.aggregates) : null;
  return {
    lead: values === null ? null
      : ctxFilter.fill(lead.template, Object.assign({}, lead.fixed, values)),
    summary: ctxFilter.summary(view, cfg, input.template, input.none, ctxFilter.moved(chosen)),
  };
})));
"""


def _python_lead(lead, rows: list[dict], chosen: dict[str, str]) -> str | None:
    """The oracle: ``ctx.viz.aggregates`` over the rows the choice keeps, or None for none."""
    from ctx.viz import aggregates  # noqa: PLC0415

    shown = [r for r in rows if all(not v or str(r[c]) == v for c, v in chosen.items())]
    values = aggregates.fill(shown, lead.aggregates) if shown else None
    return None if values is None else lead.template.format_map({**lead.fixed, **values})


@pytest.mark.parametrize("name", list(LEADS))
def test_a_lead_follows_every_filter_level_as_python_fills_it(name):
    """D1: the M1 depth bars' and M3's mapping and strand bars' leads, refilled by the page at
    every level of every filter, equal ``ctx.viz.aggregates`` over the same rows, a level with
    no bar or no number to take giving none; at "All" the lead is the one the page prints and
    the live line is the default's, and once a filter is set the live line says what the
    caption counts. D2: at "All" the made-up M1 page's median prints 0.494 and Step-1's 26.1.
    Mutations caught: the median taken over every row; the depth printed as ".1f" (0.5)."""
    data, lead = LEADS[name]()
    block = figure_block(data, "x", FigureText(title="Demo", reading=("Read.",), follows=lead))
    rows, cfg = [finite(r) for r in block.rows], dict(block.controls)
    levels = itertools.product(*(["", *v] for v in cfg["filters"].values()))
    choices = [dict(zip(cfg["filters"], combo, strict=True)) for combo in levels]
    words = _defaults(cfg)
    out = _node({"rows": rows, "cfg": cfg, "choices": choices, "template": words["template"],
                 "none": words["none"]}, LEAD, FILTER)  # fmt: skip
    for chosen, page in zip(choices, out, strict=True):
        assert page["lead"] == _python_lead(lead, rows, chosen), chosen
        if page["lead"] is not None:
            assert page["summary"].endswith(" " + lead.note) == any(chosen.values()), chosen
    assert out[0]["lead"] == block.text.lead and out[0]["summary"] == words["summary"]
    assert any(p["lead"] is None for p in out)
    assert len({p["lead"] for p in out}) > len(choices) // 3
    depth = {"depth_made_up": "0.494 million", "depth_step1": "26.1 million"}.get(name)
    assert depth is None or f"The median run has {depth} reads." in block.text.lead


def test_a_lead_names_its_first_row_on_a_tie_and_holds_three_digits():
    """The planted ties and shares the leads above are compared on: S3 and S6 tie lowest, and
    S3, first, is named; a forward share of 0.0213 prints 0.0213, not 0.021 held to three
    decimals nor 0.0210; LDN's calls tie two to two and the first seen is named."""
    data, lead = _mapping()
    rows = [finite(r) for r in data.page_table().to_dicts()]
    assert lead.text(rows) == "Median uniquely mapped: 87.9 % of input reads; lowest 75.6 % in S3."
    data, lead = _strand()
    rows = [finite(r) for r in data.page_table().to_dicts()]
    assert lead.text(rows) == "Forward share 0.0213 to 0.961: 6 of 10 samples are called reverse."
    ldn = [r for r in rows if r["cell_type"] == "LDN"]
    assert lead.text(ldn) == "Forward share 0.0312 to 0.503: 2 of 4 samples are called unstranded."


EXPORT = ROOT / "src" / "ctx" / "report" / "templates" / "_export.js.j2"
WRAP = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
console.log(JSON.stringify(input.cases.map(([text, ex]) =>
  ctxExport.wrap(text, ctxExport.fits(ex, input.metrics)))));
"""
# Step-1's PCA caption, one of the course's 99 (2026-10-07): at 85 mm the advances without
# their kerning end its third line a word early.
KERNED = (
    "PERMANOVA on Euclidean distances for condition: R\u00b2 = 0.40, pseudo-F = 8.28, "
    "p = 0.001 (999 permutations; n = 55 samples in 5 groups). Blind VST, the 500 most variable "
    "genes after the transform; PC1 33%, PC2 9% of their variance; n = 55 samples."
)


def _caption_lines(figure) -> list[str]:
    """The caption a figure file draws, line by line: the text place_caption sets at the
    insets."""
    from ctx.viz.style import CAPTION_BOTTOM, CAPTION_INSET  # noqa: PLC0415

    (caption,) = [t for t in figure.texts if t.get_position() == (CAPTION_INSET, CAPTION_BOTTOM)]
    return caption.get_text().split("\n")


def test_export_wraps_its_caption_as_the_figure_files_do():
    """The view's caption breaks where the figure files' caption breaks: the page measures a
    line with the numbers ctx.viz.metrics measures the files' captions with (the export's
    "metrics"), greedy, never splitting a word, a hyphenated one included. Checked against the
    caption each demo figure's own static() draws, and against wrap_to at both widths on a
    caption whose breaks the kerning moves and on edge cases. Mutations: the character-count
    wrap the first export used; a line closed one word early; the kerning left out of the
    page's measure; a split at a hyphen."""
    from ctx.report import export  # noqa: PLC0415
    from ctx.viz.metrics import wrap_to  # noqa: PLC0415
    from ctx.viz.style import CAPTION_INSET, SIZES, WIDTHS_MM  # noqa: PLC0415

    demos = [demo_volcano(design_counts(condition_lfc=1.0)), demo_linked(), filtered_bars()]
    cases, expected = [], []
    for data in demos:
        ex = figure_block(data, "x").controls["export"]
        cases.append([ex["caption"], ex])
        expected.append(_caption_lines(data.static(THEMES["white"])))
    texts = [KERNED, "a  b\tc " + "x" * 12 + " d", "", "one-word-wide-as-two two"]
    for column, words in [(c, t) for c in sorted(WIDTHS_MM) for t in texts[:1]] + [
            ("tiny", t) for t in texts[1:]]:  # fmt: skip
        width_pt = 30 / (1 - 2 * CAPTION_INSET) if column == "tiny" else (
            WIDTHS_MM[column] / 25.4 * 72)  # fmt: skip
        cases.append([words, {"width_pt": width_pt, "inset": CAPTION_INSET,
                              "size": SIZES["caption"]}])  # fmt: skip
        expected.append(wrap_to(words, width_pt * (1 - 2 * CAPTION_INSET), SIZES["caption"]))
    extra = frozenset("".join(words for words, _ in cases)) - export.package_characters()
    metrics = json.loads(str(export.page_assets(extra)["config"]))["metrics"]
    out = _node({"cases": cases, "metrics": metrics}, WRAP, EXPORT)
    assert out == expected
    # The cases hold what they test: the linked view's caption takes two lines, the kerned
    # caption breaks otherwise at 85 mm than at 170 mm, and the edge cases.
    assert [len(lines) for lines in expected[:3]] == [1, 2, 1]
    assert expected[3] != expected[4] and expected[5] == ["a b c", "x" * 12, "d"]
    assert expected[6] == [] and expected[7] == ["one-word-wide-as-two", "two"]


CAPTION = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
console.log(JSON.stringify(input.cases.map(([ex, view]) =>
  ctxExport.captionLines(ex, view, (line) => line.length <= 1000))));
"""


def test_export_caption_is_caption_live_line_and_note_each_closed():
    """The file's caption is the figure files' caption, the live line and the note, in that
    order, each a sentence: the linked view's live line ends with no stop, and the note ran
    on from it ("vermillion: type_a 3 Drawn in the browser"). Mutations: the note left out;
    the stop never added."""
    ex = {"caption": "No test.", "note": "Drawn in the browser."}
    cases = [[ex, "3 painted vermillion: type_a 3"], [ex, "Up: 4 genes."], [ex, ""],
             [{**ex, "caption": "A caption without a stop"}, "  "]]  # fmt: skip
    out = [lines[0] for lines in _node({"cases": cases}, CAPTION, EXPORT)]
    assert out == [
        "No test. 3 painted vermillion: type_a 3. Drawn in the browser.",
        "No test. Up: 4 genes. Drawn in the browser.",
        "No test. Drawn in the browser.",
        "A caption without a stop. Drawn in the browser.",
    ]


FIGURE = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
console.log(JSON.stringify(input.cases.map(([parts, cfg]) =>
  ctxExport.figure(parts, cfg, {journal: true}, 200, input.dashes))));
"""


def _figure(cases: list) -> list:
    """ctxExport.figure() under node on each (parts, controls), with the page's dash patterns."""
    from ctx.viz.style import plotly_dashes  # noqa: PLC0415

    return _node({"cases": cases, "dashes": plotly_dashes()}, FIGURE, EXPORT)


def _plain(figure) -> dict:
    """A plotly figure as the page holds it: arrays as lists."""
    return json.loads(
        json.dumps(
            figure.to_plotly_json(),
            default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o),
        )
    )


def _filled(figure, data, cfg) -> dict:
    """A live-calls figure as the page draws it at the default rule: its empty point traces
    filled from the rows it embeds, by the default call, as liveCalls fills them (GS.3)."""
    plain, rows = _plain(figure), data.on_page()
    for level, index in cfg["traces"].items():
        part = rows.filter(pl.col(data.colour) == level)
        plain["data"][index]["x"] = part[cfg["lfc"]].to_list()
        plain["data"][index]["y"] = [-math.log10(max(p, cfg["p_floor"])) for p in part[cfg["p"]]]
    return plain


def test_export_draws_marks_at_the_figure_files_sizes_from_any_zoom():
    """Every point MARKER_SIZE pt, every line its LINE_WIDTHS pt, from the module page and
    from a page opened alone at twice the size, each plot's zoom read off its layout's
    meta.zoom, which the page script sets to 1 below the floor; the journal template and the
    files' width; a linked view's parallel sets under the scatter. Mutation: the page's zoom
    not divided out."""
    from ctx.viz.style import LINE_WIDTHS, MARKER_SIZE, zoom  # noqa: PLC0415

    volcano = figure_block(demo_volcano(design_counts(condition_lfc=1.0)), "x")
    linked = figure_block(demo_linked(), "y")
    cases = []
    for scale in (1.0, 2.0):
        for block in (volcano, linked):
            parts = [zoom(_plain(f), scale) for f in (block.figure, block.below) if f is not None]
            for part in parts if scale != 1 else ():
                part["layout"]["meta"] = {**(part["layout"].get("meta") or {}), "zoom": scale}
            cases.append([parts, block.controls])
    out = _figure(cases)
    for (_, cfg), drawn in zip(cases, out, strict=True):
        layout, ex = drawn["layout"], cfg["export"]
        assert layout["template"] == {"journal": True}
        assert layout.get("margin") == ({"l": ex["sides"]["l"], "r": ex["sides"]["r"]}
                                        if "sides" in ex else None)  # fmt: skip
        assert (layout["width"], layout["height"]) == (ex["width_pt"], 200)
        sizes = {t["marker"]["size"] for t in drawn["data"] if t.get("mode") == "markers"}
        assert sizes and all(s == pytest.approx(MARKER_SIZE) for s in sizes), sizes
        widths = {t["line"]["width"] for t in drawn["data"] if t.get("mode") == "lines"}
        assert all(w == pytest.approx(LINE_WIDTHS["line"]) for w in widths), widths
    linked_out = out[1]
    assert [t["type"] for t in linked_out["data"]] == ["scatter", "parcats"]
    assert linked_out["data"][1]["domain"]["y"] == linked.controls["below"]
    assert linked_out["layout"]["yaxis"]["domain"] == linked.controls["above"]


def test_export_draws_the_x_title_as_built():
    """A title the page wrapped for a phone, whose titles guard keeps the title on one line in
    layout.meta.xtitle, is drawn whole, as the figure files set it; with none kept the title
    stands; meta never reaches the file. The volcano's title is built with its unit on a line
    of its own (axes.screen_title, M3.6), which the guard reads as one line. Mutation: the
    restore removed, which broke the 170 mm files' title where a phone's plot broke it
    (review, 2026-10-02)."""
    data = demo_volcano(design_counts(condition_lfc=1.0))
    block = figure_block(data, "x")
    wrapped = _plain(block.figure)
    built = wrapped["layout"]["xaxis"]["title"]["text"]
    files = data.static(THEMES["white"]).axes[0].get_xlabel()
    assert built == "Fold change<br>(log2 fold change)" and files == built.replace("<br>", " ")
    wrapped["layout"]["xaxis"]["title"]["text"] = files.replace(" ", "<br>", 1)
    kept = {**wrapped, "layout": {**wrapped["layout"], "meta": {"xtitle": files}}}
    out = _figure([[[kept], block.controls], [[wrapped], block.controls]])
    assert [o["layout"]["xaxis"]["title"]["text"] for o in out] == [
        files, files.replace(" ", "<br>", 1)]  # fmt: skip
    assert all("meta" not in o["layout"] for o in out)


TITLE = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
console.log(JSON.stringify(input.cases.map(([text, width]) =>
  ctxExport.wrap(ctxExport.words(text), (line) => line.length <= width))));
"""


def test_a_wrapped_axis_title_keeps_its_unit_whole():
    """The page's titles guard wraps an x title from its words with each parenthesised span
    held whole, as the files hold a unit (annotate.fit_title) and the page breaks before it
    (axes.screen_title, M3.6): split at its spaces, "(log2(count + 1))" left "1))" alone on a
    line. A span wider than a line takes its own. Mutation: words() splitting at every space."""
    cases = [["Mean of both (log2(count + 1))", 16], ["Fold change (log2 fold change)", 22],
             ["PC1, 42% of variance (VST units)", 12], ["Reads (%)", 4], ["", 5]]  # fmt: skip
    out = _node({"cases": cases}, TITLE, EXPORT)
    assert out == [
        ["Mean of both", "(log2(count + 1))"],
        ["Fold change", "(log2 fold change)"],
        ["PC1, 42% of", "variance", "(VST units)"],
        ["Reads", "(%)"],
        [],
    ]


COMPOSE = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
const svg = ctxExport.compose(input.plot, input.lines, input.ex);
console.log(JSON.stringify({svg, pdf: ctxExport.forPdf(svg, "CtxSans"),
                            height: ctxExport.captionHeight(input.lines.length, input.ex)}));
"""


def test_export_svg_is_sized_in_points_with_its_caption_inside():
    """The file is the figure files' size in pt, one user unit to the point; the caption's
    lines sit inside its insets, escaped; the PDF's copy sets every text in the face it
    embeds; plotly's hard-coded band edges go. Mutations: the size left in plotly's px, the
    caption's last line below the page."""
    from xml.etree import ElementTree  # noqa: PLC0415

    ex = figure_block(demo_volcano(design_counts(condition_lfc=1.0)), "x").controls["export"]
    plot = ('<svg xmlns="http://www.w3.org/2000/svg" class="main-svg" width="481" '
            'height="230" style="" viewBox="0 0 481 230">'
            "<text style=\"font-family: Arial, 'Liberation Sans', sans-serif; font-size: 9px\">"
            '1</text><path class="path" fill="rgb(1, 2, 3)" stroke="lightgray" '
            'stroke-width="0.2"/></svg>')  # fmt: skip
    lines = ["Up in case: 3 < 4 & so on.", "A second line."]
    out = _node({"plot": plot, "lines": lines, "ex": ex}, COMPOSE, EXPORT)
    root = ElementTree.fromstring(out["svg"])
    w, h = ex["width_pt"], ex["height_pt"]
    assert (root.get("width"), root.get("height")) == (f"{w}pt", f"{h}pt")
    assert root.get("viewBox") == f"0 0 {w} {h}"
    texts = [t for t in root.iter("{http://www.w3.org/2000/svg}text")][1:]
    assert [t.text for t in texts] == lines
    baselines = [float(t.get("y")) for t in texts]
    assert baselines[-1] <= h * (1 - ex["bottom"]) + 1e-9
    assert baselines[1] - baselines[0] == pytest.approx(ex["size"] * ex["spacing"])
    assert float(texts[0].get("x")) == pytest.approx(ex["inset"] * w)
    assert out["height"] == pytest.approx(2 * ex["size"] * ex["spacing"] + ex["bottom"] * h)
    assert "lightgray" not in out["svg"] and 'stroke="none"' in out["svg"]
    assert set(re.findall(r"font-family:\s*([^;\"]+)", out["pdf"])) == {"CtxSans"}


MOVED = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
const cfg = input.cfg;
console.log(JSON.stringify({
  brush: [ctxBrush.moved(ctxBrush.start(3)), ctxBrush.moved([-1, 2, -1]),
          ctxBrush.moved(ctxBrush.clear([0, 1, 2]))],
  filter: [ctxFilter.moved({a: "", b: ""}), ctxFilter.moved({a: "", b: "x"}), ctxFilter.moved({})],
  calls: [ctxCalls.moved(cfg, cfg.lfc_min, cfg.gate, cfg.alpha),
          ctxCalls.moved(cfg, Math.log2(Math.pow(2, cfg.lfc_min)), cfg.gate, cfg.alpha),
          ctxCalls.moved(cfg, 2, cfg.gate, cfg.alpha),
          ctxCalls.moved(cfg, cfg.lfc_min, "padj_holm", cfg.alpha),
          ctxCalls.moved(cfg, cfg.lfc_min, cfg.gate, 0.01)],
}));
"""


def test_a_view_is_moved_only_away_from_the_defaults():
    """Nothing painted, every filter on all and the default rule are the figure files' view;
    any paint, filter, threshold, basis or cutoff is not. Mutation: the cutoff left out of
    ctxCalls.moved, so a stricter cutoff still linked the default file."""
    cfg = figure_block(demo_volcano(design_counts(condition_lfc=1.0)), "x").controls
    scripts = [SCRIPT, BRUSH, *FILTER]
    rule = "\n".join(p.read_text().replace("{% raw %}", "").replace("{% endraw %}", "")
                     for p in scripts)  # fmt: skip
    node = shutil.which("node")
    done = subprocess.run([node, "-e", rule + MOVED], input=json.dumps({"cfg": cfg}),
                          capture_output=True, text=True, check=True)  # fmt: skip
    out = json.loads(done.stdout)
    assert out == {"brush": [False, True, False], "filter": [False, True, False],
                   "calls": [False, False, True, True, True]}  # fmt: skip


def test_export_draws_what_the_figure_files_draw_and_the_page_lacks(monkeypatch):
    """A thinned volcano's left-out points back in their trace; its labels as annotations at
    the offsets static() placed them, with a solid leader; plotly's named dashes as the
    figure files' lengths; legend keys opaque beside points at the page's opacity.
    Mutations: the thinned points dropped; the labels dropped; "dash" left as plotly's own
    9 on 9 off; the legend key left at the points' opacity."""
    from ctx.viz import scatter  # noqa: PLC0415
    from ctx.viz.style import DASH_PATTERNS, POINT_ALPHA, SIZES  # noqa: PLC0415

    monkeypatch.setattr(scatter, "THIN_MAX", 300)
    data = demo_volcano(design_counts(condition_lfc=1.0))
    block = figure_block(data, "x")
    (drawn,) = _figure([[[_filled(block.figure, data, block.controls)], block.controls]])
    traces, cfg = drawn["data"], block.controls
    points = [t for t in traces if t.get("mode") == "markers" and t["x"] != [None]]
    assert sum(len(t["x"]) for t in points) == data.frame.height
    (ns,) = [t for t in points if t["name"] == "ns"]
    assert len(ns["x"]) == data.frame.filter(pl.col("call") == "ns").height
    keys = [t for t in traces if t.get("x") == [None]]
    assert [k["name"] for k in keys] == [t["name"] for t in points]
    assert all(k["marker"]["opacity"] == 1 and k["legendgroup"] == t["legendgroup"]
               and t["showlegend"] is False and t["marker"]["opacity"] == POINT_ALPHA
               for k, t in zip(keys, points, strict=True))  # fmt: skip
    dashed = {t["line"]["dash"] for t in traces if t.get("mode") == "lines"}
    assert dashed == {",".join(f"{v}px" for v in DASH_PATTERNS["dashed"])}
    labels = drawn["layout"]["annotations"]
    placed = cfg["export"]["labels"]
    assert [(a["text"], a["x"], a["y"], a["ax"], -a["ay"], a["xanchor"]) for a in labels] == [
        (p["text"], p["x"], p["y"], p["dx"], p["dy"], p["anchor"]) for p in placed]  # fmt: skip
    assert placed, "the static file labels some genes"
    assert {a["font"]["size"] for a in labels} == {SIZES["annotation"]}
    assert {a["axref"] for a in labels} == {"pixel"} and {a["arrowhead"] for a in labels} == {0}
