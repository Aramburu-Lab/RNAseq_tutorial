"""``ctx clinical`` end to end on a made-up GDC root (tests/gdc_fixture.py): the page, its tables,
and the live panel's default state against the build's own CSVs (GS plan, test 7). Each test
names the mutation it catches."""

import base64
import gzip
import json
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from clinical_root import PATIENTS, ctx, make

from ctx.report.builder import MAX_BYTES, plotly_js

ROOT = Path(__file__).resolve().parents[1]
TEMPLATES = ROOT / "src" / "ctx" / "report" / "templates"


@pytest.fixture(scope="module")
def clinical_root(tmp_path_factory):
    """A gathered made-up GDC root with ctx clinical run on it; returns the root."""
    root, spec = make(tmp_path_factory.mktemp("m10"))
    ctx("--root", str(root), "clinical", "--dataset", "fake_brca", "--dataset-file", str(spec))
    return root


def test_genes_from_replaces_the_build_list(tmp_path):
    """Catches --genes-from read and then ignored: the page scores the file's list."""
    root, spec = make(tmp_path)
    listed = tmp_path / "mine.txt"
    listed.write_text("# source: x:y\n# count: 2\nG2\nG7\n")
    ctx("--root", str(root), "clinical", "--dataset", "fake_brca", "--dataset-file", str(spec),
        "--genes-from", str(listed))  # fmt: skip
    html = outputs(root, "overview_v1.html").read_text()
    cfg = controls(html, "survival")
    assert cfg["defaults"]["genes"] == "G2\nG7" and cfg["source"] == "mine.txt"


def outputs(root: Path, what: str) -> Path:
    """The one clinical output of this kind."""
    (path,) = (root / "outputs" / "12_clinical").glob(f"*_clinical_fake_brca_{what}")
    return path


def controls(html: str, kind: str) -> dict:
    """The first figure's page controls of this kind."""
    for body in re.findall(r'id="ctx-controls-\d+">(.*?)</script>', html, re.S):
        found = json.loads(body)
        if found["kind"] == kind:
            return found
    raise AssertionError(f"no {kind} controls")


def test_the_page_is_offline_under_the_cap_and_names_its_sources(clinical_root):
    """Catches a page over 15 MB, a script fetched from the network, and the GDC
    acknowledgement or the CDR's licence left out."""
    html = outputs(clinical_root, "overview_v1.html").read_text()
    assert len(html.encode()) <= MAX_BYTES
    assert html.count(plotly_js()) == 1
    assert not re.search(r"<script[^>]+\bsrc\s*=", html, re.IGNORECASE)
    assert "TCGA Research Network" in html and "CC BY-NC-ND 4.0" in html
    assert "ctxLive" in html and "ctxSurv" in html


def test_every_patient_once_and_every_table_has_its_readme(clinical_root):
    """Catches a patient with two tumours (the second aliquot of none here, but the normals of
    patients 1 to 5) and a table without its README."""
    patients = pl.read_parquet(outputs(clinical_root, "patients_v1.parquet"))
    assert patients.height == PATIENTS and patients["patient"].n_unique() == PATIENTS
    assert patients["sample"].str.contains("-01A-").all()
    for table in (clinical_root / "outputs" / "12_clinical").glob("*.parquet"):
        assert table.with_name(table.name.replace(".parquet", "_README.md")).is_file()


HARNESS = r"""
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
const cfg = input.cfg;
const panel = ctxSurv.planes(Uint8Array.from(input.bytes), cfg.genes, cfg.patients);
const signs = cfg.build.map(() => 1);
const found = ctxSurv.analyse(cfg, panel, cfg.build, signs, cfg.endpoint, cfg.split);
console.log(JSON.stringify({
  alone: found.alone, adjusted: found.adjusted, corrected: found.scan.corrected,
  groups: found.groups.map(g => ({name: g.name, patients: g.patients, events: g.events,
                                  median: g.median[0], rmst: g.rmst})),
  summary: ctxLive.fill(cfg.words.summary, ctxLive.values(cfg, found,
    {endpoint: cfg.endpoint, split: cfg.split, genes: cfg.build, source: cfg.source}, cfg.none))
}));
"""


def node_default(html: str) -> dict:
    """The panel's default state as the page computes it, under node."""
    binary = shutil.which("node")
    assert binary, "nodejs is pinned in the ctx-dev environment; run tests with pixi"
    cfg = controls(html, "survival")
    cfg["split_words"] = dict(cfg["defaults"]["splits"])
    raw = gzip.decompress(base64.b64decode(cfg.pop("block")))
    script = "".join(re.sub(r"\{%-? (?:end)?raw -?%\}", "", (TEMPLATES / n).read_text())
                     for n in ("_surv.js.j2", "_live.js.j2"))  # fmt: skip
    done = subprocess.run([binary, "-e", script + HARNESS],
                          input=json.dumps({"cfg": cfg, "bytes": list(raw)}),
                          capture_output=True, text=True, check=True)  # fmt: skip
    return json.loads(done.stdout)


def test_the_live_default_equals_the_builds_tables(clinical_root):
    """GS plan test 7: with no list given, the page's own computation of the build list gives
    the build's Cox and group tables and its summary line. Catches a default list, endpoint or
    split that differs from the build's, the block decoded in the wrong plane order, and the
    adjusted model fitted over rows a covariate leaves empty."""
    html = outputs(clinical_root, "overview_v1.html").read_text()
    page = node_default(html)
    models = pl.read_csv(outputs(clinical_root, "cox_score_v1.csv"))
    for name in ("alone", "adjusted"):
        row = models.filter(pl.col("model") == name).row(0, named=True)
        got = page[name]
        assert (got["patients"], got["events"]) == (row["patients"], row["events"])
        assert got["hr"] == pytest.approx(row["hr"], rel=1e-6)
        assert got["p"] == pytest.approx(row["p"], rel=1e-5)
        assert got["c_index"] == pytest.approx(row["c_index"], rel=1e-9)
    groups = pl.read_csv(outputs(clinical_root, "groups_score_v1.csv"))
    for got, row in zip(page["groups"], groups.iter_rows(named=True), strict=True):
        assert (got["name"], got["patients"], got["events"]) == (row["group"], row["patients"],
                                                                 row["events"])  # fmt: skip
        assert got["rmst"] == pytest.approx(row["rmst"], rel=1e-12)
    summary = re.search(r'<p class="survival-summary"[^>]*>(.*?)</p>', html, re.S)
    assert summary is not None
    assert page["summary"] == summary.group(1).replace("&#39;", "'")


def test_the_risk_rows_count_patients_followed_at_each_time():
    """At-risk counts are those of the first distinct time at or after each count time, 0 past
    the last; the band holds each bound until the next time. Catches the count taken after the
    time (t > at), and a band joining bounds diagonally."""
    from ctx.viz.km import COLUMNS, KmData  # noqa: PLC0415

    frame = pl.DataFrame([("a", 1.0, 10, 1, 0, 0.9, 0.7, 0.98), ("a", 2.0, 9, 2, 1, 0.7, 0.5, 0.85),
                          ("a", 5.0, 6, 1, 2, 0.58, 0.4, 0.75)], schema=list(COLUMNS),
                         orient="row")  # fmt: skip
    data = KmData(frame, {c: ("none", c) for c in COLUMNS}, ("a",), {"a": "#0072B2"}, "t", "s",
                  "risk", (0.0, 2.0, 4.0, 6.0), "caption")  # fmt: skip
    assert data.at_risk("a") == [10, 9, 6, 0]
    x, y = data.band("a")
    assert x[:4] == [0.0, 1.0, 1.0, 2.0] and y[:4] == [1.0, 1.0, 0.98, 0.98]


def test_the_risk_step_is_the_same_in_both_languages():
    """At most seven at-risk columns of whole years, as the page redraws them. Catches the twin
    rounding the step another way, or a step that is not whole years."""
    from ctx.clinical.figures import risk_step  # noqa: PLR0402, PLC0415
    from ctx.viz.style import KM  # noqa: PLC0415

    lasts = [5.0, 59.9, 72.0, 72.1, 283.4, 1000.0]
    script = re.sub(r"\{%-? (?:end)?raw -?%\}", "", (TEMPLATES / "_live.js.j2").read_text())
    probe = (
        f"console.log(JSON.stringify({lasts}.map(l => ctxLive.riskStep(l, {KM['risk_columns']}))));"
    )
    done = subprocess.run([shutil.which("node") or "node", "-e", script + probe],
                          capture_output=True, text=True, check=True)  # fmt: skip
    assert json.loads(done.stdout) == [risk_step(v) for v in lasts]
    assert [risk_step(v) for v in (72.0, 72.1, 283.4)] == [12.0, 24.0, 48.0]


def test_the_block_holds_the_builds_stored_values(clinical_root):
    """The panel's block decodes to the stored log2 CPM of the build list's genes, in order.
    Catches values scaled twice, or the planes swapped."""
    from ctx.clinical import panel  # noqa: PLC0415

    html = outputs(clinical_root, "overview_v1.html").read_text()
    cfg = controls(html, "survival")
    values = panel.decode(cfg["block"], cfg["genes"], cfg["patients"])
    box = pl.read_csv(outputs(clinical_root, "expression_er_v1.csv"))
    first = box.filter(pl.col("feature") == cfg["symbols"][0])
    patients = pl.read_parquet(outputs(clinical_root, "patients_v1.parquet"))["patient"].to_list()
    stored = dict(zip(patients, values[0].tolist(), strict=True))
    for sample, value in first.select("sample", "value").rows():
        assert round(value * 1000) == stored[sample]
    assert np.asarray(values).dtype == np.int16
