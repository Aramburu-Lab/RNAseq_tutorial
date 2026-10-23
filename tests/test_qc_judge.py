"""M2's judge panel (work-order item 5): ``ctxJudge`` held to ``ctx.qc.checkpoint.judge``
under node (D-012), over every k step and strata grouping the panel offers, on values that
reach each rule: small strata, a stratum with no spread, both sides, a keep metric, a banded
metric, a log10 metric, a strand call and a metric judged within its read-length group; and
the numbers its reasons print, as Python prints them.

Each test names the mutation it catches.
"""

import itertools
import math
from dataclasses import replace
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from test_export_browser import CHROMIUM, _browse
from test_page_script import _node
from test_qc_views import page  # noqa: F401 - the fixture: ctx qc on the rnaseq reports

from ctx.qc import checkpoint, judge
from ctx.qc import registry as reg
from ctx.report import text

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src" / "ctx" / "report" / "templates" / "_qc.js.j2"
# One metric of each rule the twin decides, from the registry the run reads.
METRICS = ("library_size", "genes_detected", "pca_distance", "percent_gc",
           "percent_duplicates", "tg_total_reads", "salmon_strand_bias", "read_length")  # fmt: skip
GROUPINGS = (("cell_type",), (), ("cell_type", "condition"), ("batch",))
HARNESS = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
console.log(JSON.stringify(input.cases.map((c) => {
  const labels = ctxJudge.labels(input.units, "sample", c.cols, input.words.none);
  const found = ctxJudge.judge(c.values, labels, c.metric, c.k, input.cfg, c.absolute, c.calls,
                               input.words);
  return {bands: found.bands, reasons: found.reasons};
})));
"""


def _values(seed: int = 11, n: int = 40) -> pl.DataFrame:
    """Samples in strata of 1, 2 and many, one stratum at one value, missing values and NaN."""
    rng = np.random.default_rng(seed)
    cell = ["a"] * 20 + ["b"] * 14 + ["c"] * 3 + ["d"] * 2 + ["e"]
    frame = pl.DataFrame({
        "sample": [f"S{i:03d}" for i in range(n)], "cell_type": cell,
        "condition": [["x", "y"][i % 2] for i in range(n)],
        "batch": [["long", "short"][i % 3 == 0] for i in range(n)],
        checkpoint.LENGTH_GROUP: [["75 nt", "101 nt"][i % 4 == 0] for i in range(n)],
        "library_size": rng.lognormal(16, 0.6, n), "genes_detected": rng.normal(15000, 900, n),
        "pca_distance": np.abs(rng.normal(0, 3, n)), "percent_gc": rng.normal(50, 2, n),
        "percent_duplicates": rng.uniform(10, 90, n), "tg_total_reads": rng.uniform(2e4, 2e5, n),
        "salmon_strand_bias": rng.choice([0.02, 0.5, 0.97, 0.7, 0.51], n),
        "read_length": rng.choice([75.0, 76.0, 101.0], n),
    })  # fmt: skip
    frame = frame.with_columns(
        pl.when(pl.col("sample") == "S001").then(1e3).otherwise(pl.col("library_size"))
        .alias("library_size"),
        # In the bad band and far under its stratum: both rules flag it, the MAD's reason last.
        pl.when(pl.col("sample") == "S002").then(500.0).otherwise(pl.col("tg_total_reads"))
        .alias("tg_total_reads"),
        pl.when(pl.col("cell_type") == "c").then(12000.0).otherwise(pl.col("genes_detected"))
        .alias("genes_detected"),
        pl.when(pl.col("sample").is_in(["S004", "S009"])).then(None)
        .otherwise(pl.col("pca_distance")).alias("pca_distance"),
        pl.when(pl.col("sample") == "S005").then(float("nan"))
        .otherwise(pl.col("percent_gc")).alias("percent_gc"),
    )  # fmt: skip
    # Each metric's values under the column the registry reads it from.
    registry = reg.load(reg.resolve(None))
    return frame.rename({m: registry.metrics[m].column for m in METRICS})


def _words() -> dict:
    return {**{f: text.raw("qc", f"reasons.{f}") for f in ("low", "high", "strand")},
            "small": text.raw("qc", "notes.small_stratum"),
            "spread": text.word("qc", "notes.no_spread"),
            "none": text.word("qc", "strata.none")}  # fmt: skip


def test_the_page_judges_as_checkpoint_judge_does():
    """Mutations caught: the MAD unscaled; the strand call's tie broken by count order; the
    bad side ignored; a keep metric flagged; the absolute reason winning over the MAD's."""
    registry = reg.load(reg.resolve(None))
    frame = _values()
    samples = frame.select("sample", "cell_type", "condition", "batch", checkpoint.LENGTH_GROUP)
    grouped = ["cell_type", "condition", "batch", checkpoint.LENGTH_GROUP]
    units = judge.unit_table(frame.select("sample"), samples, "sample", grouped)
    index = {s: i for i, s in enumerate(frame["sample"])}
    cases, expected = [], []
    for k, strata in itertools.product(judge.K_STEPS, GROUPINGS):
        bands, flags = checkpoint.judge(frame.drop(samples.columns[1:]), samples,
                                        replace(registry, mad_multiplier=k), strata)  # fmt: skip
        for metric_id in METRICS:
            metric = registry.metrics[metric_id]
            cut = frame.select("sample", metric.column)
            absolute, calls = judge.decided(cut, metric)
            within = metric.judged_within == checkpoint.LENGTH_GROUP
            values = [None if v is None or math.isnan(v) else v for v in cut[metric.column]]
            cases.append({
                "cols": [checkpoint.LENGTH_GROUP] if within else list(strata), "values": values,
                "metric": judge._metric(metric), "k": k,  # noqa: SLF001
                "absolute": {str(i): r for i, r in absolute.items()}, "calls": calls,
            })  # fmt: skip
            mine = flags.filter(pl.col("metric") == metric.id).select("sample", "reason")
            expected.append((bands.filter(pl.col("metric") == metric.id).drop("metric").to_dicts(),
                             {str(index[s]): r for s, r in mine.rows()}))  # fmt: skip
    got = _node({"units": {c: units[c].to_list() for c in units.columns}, "cases": cases,
                 "cfg": {"min_samples": registry.min_samples, "scaled": registry.mad_scaled},
                 "words": _words()}, HARNESS, SCRIPT)  # fmt: skip
    seen = set()
    for case, (bands, reasons), found in zip(cases, expected, got, strict=True):
        where = (case["metric"]["id"], case["k"], tuple(case["cols"]))
        assert found["reasons"] == reasons, where
        assert len(found["bands"]) == len(bands), where
        for mine, theirs in zip(found["bands"], bands, strict=True):
            for key, value in theirs.items():
                if isinstance(value, float):
                    assert mine[key] == pytest.approx(value, rel=1e-12), (where, key)
                else:
                    assert mine[key] == value, (where, key)
        seen |= {"flagged" if reasons else "clean"} | {b["note"][:8] for b in bands if b["note"]}
    assert {"flagged", "clean", "no band:"} <= seen


VALUES = (
    0.125,
    0.0125,
    2.5,
    1234.5,
    1e-5,
    123456.789,
    -0.0004567,
    99.95,
    100.0,
    1e21,
    0.5,
    5e-324,
    1.0000000000000002,
    999.9996,
    9.995,
    -1234567.5,
    0.00012345,
    1e-4,
)


def test_numbers_print_as_python_prints_them():
    """Every value at every digit count a reason uses, ties included (Python rounds half to
    even, JavaScript's own formatting rounds them up). Mutation caught: ties rounded up."""
    harness = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
console.log(JSON.stringify(input.values.map((v) => input.digits.map((d) => ctxJudge.show(v, d)))));
"""
    digits = list(range(3, 16))
    got = _node({"values": list(VALUES), "digits": digits}, harness, SCRIPT)
    for value, printed in zip(VALUES, got, strict=True):
        assert printed == [checkpoint._show(value, d) for d in digits], value  # noqa: SLF001


PANEL = """(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(1500);
  const forms = Array.from(document.querySelectorAll('form.judge'));
  const fig = forms[0].closest('.fig');
  fig.setAttribute('data-probe', '1');
  const read = () => {
    const d = fig.querySelector('.js-plotly-plot'), pts = d.data[d.data.length - 1];
    const symbols = Array.from(pts.marker.symbol);
    return {title: fig.querySelector('.fig-title').textContent,
            flagged: Array.from(pts.customdata).filter((c, i) => symbols[i] === 'triangle-up')
              .map(c => c[1]),
            ks: forms.map(f => f.elements.k.value),
            moved: !!fig.querySelector('a[data-export][role=button]'),
            line: fig.querySelector('.judge-summary').textContent};
  };
  const set = async (k) => {
    forms[0].elements.k.value = k;
    forms[0].dispatchEvent(new Event('change', {bubbles: true}));
    await wait(600);
  };
  await set('2');
  const out = {moved: read()};
  await set('3');
  out.back = read();
  await set('2');
  return out;
})()"""


@pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
def test_the_panel_rejudges_every_strip_and_comes_back(page, tmp_path):  # noqa: F811
    """In Chromium at 1280 px: k 2 moves every panel of the scope, the strip's triangles are
    the samples checkpoint.judge flags at k 2, the line names it, and k 3 draws the build.
    Mutations caught: the panels left unsynchronised; the build not restored."""
    _, out, sheet, html = page
    seen = _browse([{"name": "j", "url": html.as_uri(), "width": 1280, "height": 900,
                     "act": PANEL, "figure": "[data-probe]"}], tmp_path)["j"]  # fmt: skip
    assert seen["errors"] == []
    moved, back = seen["acted"]["moved"], seen["acted"]["back"]
    assert set(moved["ks"]) == {"2"} and set(back["ks"]) == {"3"}
    assert moved["moved"] and not back["moved"] and back["line"] == ""
    registry = reg.load(reg.resolve(None))
    metric = next(m for m in registry.metrics.values() if m.title == moved["title"])
    values = pl.read_csv(next(out.glob("*_desynth_metrics_v1.csv")))
    samples = checkpoint.length_groups(sheet, declared=False).unique(
        "sample", keep="first", maintain_order=True)  # fmt: skip
    _, flags = checkpoint.judge(values, samples, replace(registry, mad_multiplier=2.0),
                                ("cell_type",))  # fmt: skip
    expected = set(flags.filter(pl.col("metric") == metric.id)["sample"])
    assert sorted(moved["flagged"]) == sorted(expected)
    assert moved["line"].startswith("At k 2 within cell_type: ")
    assert seen["files"]["svg"]["name"].endswith("_view.svg")
