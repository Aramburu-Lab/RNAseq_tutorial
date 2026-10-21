"""M10.G: one gene's median split on the ER figures. ctx.clinical.gene.split against R's
survival (tests/fixtures/survival/gene_split.json, written by make_gene_split.R under R 4.4.2
and survival 3.7.0), and its page twin ctxSurv.gene against it, under node, on the R cohort and
on the made-up GDC root's page (tests/clinical_root.py), four genes and four endpoints. Each test
names the mutation it catches.

Tolerances: the Kaplan-Meier estimate to 1e-12 relative and its bands to 1e-10, as
tests/test_survival.py holds them; the log-rank chi-squared to 1e-12 between Python and the
twin, 1e-10 against R; its p to 1e-6 between Python and the twin, whose chi-squared tail is an
erfc approximation of fractional error under 1.2e-7, and 1e-10 against R. Group membership and
at-risk counts are exact.
"""

import base64
import gzip
import json
import re
import shutil
import subprocess
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
import pytest
from test_clinical_page import clinical_root, controls, outputs  # noqa: F401

from ctx.clinical import gene, panel
from ctx.dataset import ENDPOINTS
from ctx.report.builder import unpacked
from ctx.viz.style import KM

ROOT = Path(__file__).resolve().parents[1]
TEMPLATES = ROOT / "src/ctx/report/templates"
R = json.loads((ROOT / "tests/fixtures/survival/gene_split.json").read_text())
A = json.loads((ROOT / "tests/fixtures/survival/fixture.json").read_text())["a"]
COLUMNS = {"time": "time", "at_risk": "n_risk", "events": "n_event", "censored": "n_censor",
           "survival": "surv", "lower": "lower", "upper": "upper"}  # fmt: skip


def nan(values: list[Any]) -> np.ndarray:
    """R's NA and node's NaN, both written as text, as NaN."""
    return np.array([np.nan if v in ("NA", "NaN", None) else v for v in values], float)


def cohort() -> pl.DataFrame:
    """Cohort A as a patient table with one endpoint, its times read as months."""
    return pl.DataFrame({"OS": [r["event"] for r in A],
                         "OS_months": [float(r["time"]) for r in A]})  # fmt: skip


def rows(found: gene.GeneSplit, side: str) -> pl.DataFrame:
    """One half's Kaplan-Meier steps as the figure's contract holds them."""
    return found.figure.frame.filter(pl.col("group") == side)


@pytest.mark.parametrize("name", ["x1", "x2", "tied", "floor"])
def test_the_split_equals_rs(name):
    """Python's split against R's, on a continuous gene, on two with ties at the median, and on
    one whose median is its lowest value. Mutations caught: the split at > rather than >= (x2
    and tied put 108 and 102 of 160 above, R's rule), the log-rank's variance without the
    ties factor, an at-risk count read from the last time before a tick, and a test where one
    side is empty."""
    r = R[name]
    found = gene.split(cohort(), "OS", np.array(r["values"], float))
    assert found.upper.tolist() == r["upper"]
    for side in gene.NAMES:
        got = rows(found, side)
        if r[side] is None:
            assert got.height == 0
            continue
        for ours, theirs in COLUMNS.items():
            if ours in ("lower", "upper"):
                # R leaves the log-log band NA where the estimate is still 1, a half whose
                # first patients are censored at time 0; kaplan_meier draws it at 1 there.
                values = got[ours].cast(pl.Float64).fill_null(np.nan).to_numpy()
                theirs_ = nan(r[side][theirs])
                one = np.isnan(theirs_) & (got["survival"].to_numpy() == 1)
                assert (values[one] == 1).all()
                np.testing.assert_allclose(values[~one], theirs_[~one], rtol=1e-10)
            else:
                np.testing.assert_allclose(got[ours].to_numpy(), r[side][theirs], rtol=1e-12)
        at = replace(found.figure, risk_times=tuple(float(t) for t in R["ticks"]))
        assert at.at_risk(side) == r[side]["risk"]
    if r["logrank"] is None:
        assert found.logrank is None
    else:
        chi, df, p = found.logrank
        assert df == 1
        assert chi == pytest.approx(r["logrank"]["chisq"], rel=1e-10)
        assert p == pytest.approx(r["logrank"]["p"], rel=1e-10)


HARNESS = r"""
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
const out = input.cases.map(c => {
  const found = ctxSurv.gene(c.cfg, c.values, c.endpoint);
  let last = 0;
  found.groups.forEach(g => {
    if (g.patients) last = Math.max(last, g.curve.time[g.curve.time.length - 1]);
  });
  const step = ctxLive.riskStep(last, c.cfg.risk_columns), ticks = [];
  for (let t = 0; t <= last; t += step) ticks.push(t);
  return {upper: found.upper, logrank: found.logrank, ticks,
          groups: found.groups.map(g => ({patients: g.patients, events: g.events, curve: g.curve,
                                         risk: ctxSurv.atRisk(g.curve, ticks)}))};
});
console.log(JSON.stringify(out, (k, v) => (typeof v === "number" && !isFinite(v) ? String(v) : v)));
"""


def twin(cases: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """ctxSurv.gene on each case under node, with the ticks and counts the figure draws."""
    binary = shutil.which("node")
    assert binary, "nodejs is pinned in the ctx-dev environment; run tests with pixi"
    script = "".join(re.sub(r"\{%-? (?:end)?raw -?%\}", "", (TEMPLATES / n).read_text())
                     for n in ("_surv.js.j2", "_live.js.j2"))  # fmt: skip
    done = subprocess.run([binary, "-e", script + HARNESS], input=json.dumps({"cases": cases}),
                          capture_output=True, text=True, check=True)  # fmt: skip
    return json.loads(done.stdout)


def same(found: gene.GeneSplit, got: dict[str, Any]) -> None:
    """The twin's split, curves, log-rank and at-risk rows equal Python's."""
    assert got["upper"] == found.upper.tolist()
    assert got["ticks"] == list(found.figure.risk_times)
    for side, g in zip(gene.NAMES, got["groups"], strict=True):
        ours = rows(found, side)
        assert g["patients"] == int((found.upper == gene.NAMES.index(side)).sum())
        if not ours.height:
            assert g["patients"] == 0
            continue
        curve = g["curve"]
        assert curve["time"] == ours["time"].to_list()
        assert curve["at_risk"] == ours["at_risk"].to_list()
        assert curve["censored"] == ours["censored"].to_list()
        np.testing.assert_allclose(curve["survival"], ours["survival"].to_numpy(), rtol=1e-12)
        for bound in ("lower", "upper"):
            np.testing.assert_allclose(
                nan(curve[bound]),
                ours[bound].cast(pl.Float64).fill_null(np.nan).to_numpy(),
                rtol=1e-10,
            )
        assert g["risk"] == found.figure.at_risk(side)
    if found.logrank is None:
        assert got["logrank"] is None
    else:
        assert got["logrank"]["chi2"] == pytest.approx(found.logrank[0], rel=1e-12)
        assert got["logrank"]["p"] == pytest.approx(found.logrank[2], rel=1e-6)


def test_the_twin_equals_python_on_the_r_cohort():
    """ctxSurv.gene on cohort A, days read as months. Mutations caught: the twin's split at >
    (x2, tied), its log-rank variance without the ties factor, its at-risk count off by one
    time, a side drawn when empty (floor)."""
    names = ("x1", "x2", "tied", "floor")
    cfg = {"patients": len(A), "days_per_month": 1, "z": 1.959963984540054,
           "risk_columns": KM["risk_columns"],
           "endpoints": {"OS": {"days": [float(r["time"]) for r in A],
                                "event": [r["event"] for r in A]}}}  # fmt: skip
    got = twin([{"cfg": cfg, "values": R[n]["values"], "endpoint": "OS"} for n in names])
    for name, result in zip(names, got, strict=True):
        found = gene.split(cohort(), "OS", np.array(R[name]["values"], float))
        same(found, result)


def page_patients(cfg: dict[str, Any]) -> pl.DataFrame:
    """The page's own endpoints as a patient table: each flag, and its days in months."""
    columns: dict[str, list[Any]] = {}
    for e in ENDPOINTS:
        data = cfg["endpoints"][e]
        columns[e] = data["event"]
        columns[f"{e}_months"] = [None if d is None else d / cfg["days_per_month"]
                                  for d in data["days"]]  # fmt: skip
    return pl.DataFrame(columns, schema={k: pl.Float64 if k.endswith("_months") else pl.Int64
                                         for k in columns})  # fmt: skip


def test_the_twin_equals_python_on_the_page(clinical_root):  # noqa: F811 - the fixture
    """The made-up root's page: four of its panel genes on all four endpoints, from the
    page's own block and endpoints, so both languages read the bytes the browser reads.
    Mutations caught as above, and the twin reading another gene's plane."""
    cfg = controls(unpacked(outputs(clinical_root, "overview_v1.html").read_text()), "survival")
    values = panel.decode(cfg["block"], cfg["genes"], cfg["patients"])
    raw = gzip.decompress(base64.b64decode(cfg["block"]))
    assert len(raw) == 2 * cfg["genes"] * cfg["patients"]
    patients = page_patients(cfg)
    picked = [0, 1, cfg["genes"] // 2, cfg["genes"] - 1]
    lite = {k: cfg[k] for k in ("patients", "days_per_month", "z", "risk_columns", "endpoints")}
    cases = [{"cfg": lite, "values": values[g].tolist(), "endpoint": e}
             for g in picked for e in ENDPOINTS]  # fmt: skip
    got = twin(cases)
    for case, result in zip(cases, got, strict=True):
        found = gene.split(patients, case["endpoint"], np.array(case["values"]))
        same(found, result)
    assert any(r["logrank"] for r in got)


def test_the_screen_rows_carry_the_screens_own_p_values(clinical_root):  # noqa: F811
    """The gene box prints each model's p from the screen's rows: they are the screen table's
    p values as the page prints numbers. Mutation caught: the p cells written as none."""
    cfg = controls(unpacked(outputs(clinical_root, "overview_v1.html").read_text()), "survival")
    table = pl.read_csv(outputs(clinical_root, "cox_screen_v1.csv"), comment_prefix="#")
    printed = [[format(r["p_alone"], ".2g"), format(r["p_adjusted"], ".2g")]
               for r in table.iter_rows(named=True)]  # fmt: skip
    assert [r[4] for r in cfg["screen"]] == printed
