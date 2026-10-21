"""``ctx clinical`` end to end on a made-up GDC root (tests/gdc_fixture.py): the page, its tables,
and the live panel's default state against the build's own CSVs (GS plan, test 7). Each test
names the mutation it catches."""

import base64
import gzip
import json
import re
import shutil
import subprocess
import tomllib
import types
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from clinical_root import CLINICAL, PATIENTS, ctx, make
from gdc_fixture import PROJECT, SPEC

from ctx.clinical import ClinicalError
from ctx.clinical.run import build_list
from ctx.dataset import load
from ctx.report import genes
from ctx.report.builder import plotly_js
from ctx.report.cap import MAX_BYTES

ROOT = Path(__file__).resolve().parents[1]
TEMPLATES = ROOT / "src" / "ctx" / "report" / "templates"


@pytest.fixture(scope="session")
def clinical_root(clinical_run):
    """A gathered made-up GDC root with ctx clinical run on it; returns the root. The session's
    run, which test_cli and test_outputs_contract also read."""
    return clinical_run[0]


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


def test_a_registered_gene_the_table_lacks_stops_the_run(tmp_path):
    """D-063's KIRC list names NDUFA4L2, which GENCODE 50 calls COXFA4L2: a name the table
    lacks stops the run, named. Mutation caught: the name dropped in silence, the page scoring
    three genes where the dataset file registered four."""
    ids = [f"ENSG{n:011d}" for n in (1, 3, 185633)]
    table = genes.table(ids, ["G1", "G3", "NDUFA4L2"])
    spec_file = tmp_path / "fake_brca.toml"
    spec_file.write_text(SPEC.format(protocol="unstranded", cdr="")
                         + CLINICAL.replace('"G5"', '"NDUFA4L2"'))  # fmt: skip
    rows, _ = build_list(load(spec_file), None, table)
    assert rows == [0, 1, 2]
    spec_file.write_text(SPEC.format(protocol="unstranded", cdr="")
                         + CLINICAL.replace('"G5"', '"COXFA4L2"'))  # fmt: skip
    with pytest.raises(ClinicalError, match=r"\['COXFA4L2'\] name no expressed gene"):
        build_list(load(spec_file), None, table)


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
    assert "TCGA Research Network: https://www.cancer.gov/tcga." in html
    assert "CC BY-NC-ND 4.0" in html
    assert f"Data of {PROJECT} were retrieved from the NCI Genomic Data Commons" in html
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


CONTENT = ROOT / "src" / "ctx" / "report" / "content"
# Label words of the dataset files that name no cancer type alone ("clear" is ui.toml's button).
NOT_A_CANCER = {"tcga", "and", "cell", "clear"}


def cancer_words() -> set[str]:
    """Every word naming a cancer type in a dataset file M10 builds: its CDR type and the words
    of its label, read from conf/datasets so that no cohort identifier is typed here."""
    words = set()
    for path in (ROOT / "conf" / "datasets").glob("*.toml"):
        raw = tomllib.loads(path.read_text())
        if "clinical" in raw and raw.get("gdc", {}).get("cdr_type"):
            words.add(raw["gdc"]["cdr_type"])
            words |= {w for w in raw["dataset"]["label"].split() if w.lower() not in NOT_A_CANCER}
    return words


def strings(node) -> list[str]:
    """Every string of a parsed content file."""
    if isinstance(node, str):
        return [node]
    items = node.values() if isinstance(node, dict) else node if isinstance(node, list) else ()
    return [s for item in items for s in strings(item)]


def test_no_shared_page_text_names_a_cancer_type():
    """R1.2: a content file's text prints on every cohort's page, so none names a cancer type;
    a cohort's own words come from its dataset file. Catches the Question asking why OS is
    advised against "for breast cancer", and any other cohort's name typed into a shared key."""
    words = cancer_words()
    assert len(words) >= 8  # noqa: PLR2004 - five cohorts' types and their labels' words
    named = [(path.name, w) for path in sorted(CONTENT.glob("*.toml"))
             for s in strings(tomllib.loads(path.read_text())) for w in words
             if re.search(rf"(?<![\w-]){re.escape(w)}(?![\w-])", s,
                          0 if w.isupper() else re.I)]  # fmt: skip
    assert named == []


def note_context(marks: dict | None, endpoint: str = "PFI", default: str = "PFI"):
    """A page context holding only what ``endpoint_note`` reads."""
    from ctx.dataset import ClinicalSpec  # noqa: PLC0415

    published = {"table3": {"source": "s", "OS": [9, 90], "PFI": [7, 70], "DFI": [5, 50],
                            "DSS": [3, 30]}}  # fmt: skip
    if marks is not None:
        published["table3_use"] = {"source": "s", **marks}
    clinical = ClinicalSpec(("G1",), default, published)
    gdc = types.SimpleNamespace(cdr_type="XXXX", project="FAKE-XXXX")
    spec = types.SimpleNamespace(clinical=clinical, gdc=gdc, path=Path("fake.toml"))
    return types.SimpleNamespace(spec=spec, endpoint=endpoint)


MARKS = {"OS": "caution", "PFI": "recommended", "DFI": "not_recommended", "DSS": "caution",
         "caution": "a reason"}  # fmt: skip


def test_the_endpoint_note_is_the_dataset_files_table_3():
    """The endpoint, who chose it, and Table 3's mark and counts for it, from the dataset file;
    then its cautions and refusals. Catches the counts of another endpoint, a caution or
    refusal left out, the --endpoint choice called the dataset file's, and a mark misspelt."""
    from ctx.clinical.page import endpoint_note  # noqa: PLC0415

    first, second = endpoint_note(note_context(MARKS))
    assert first == (
        "The models below fit PFI, the dataset file's endpoint. Liu et al. 2018 "
        "Table 3 marks PFI recommended for XXXX, with 7 events and 70 censored."
    )
    assert second == ('Table 3 marks OS and DSS for caution on XXXX: "a reason". '
                      "Table 3 does not recommend DFI for XXXX.")  # fmt: skip
    first, _ = endpoint_note(note_context(MARKS, endpoint="OS"))
    assert (
        first.startswith(
            "The models below fit OS, chosen by --endpoint over the dataset file's PFI."
        )
        and "OS recommended with caution" in first
    )
    assert first.endswith("with 9 events and 90 censored.")
    every = {e: "recommended" for e in ("OS", "PFI", "DFI", "DSS")}
    assert endpoint_note(note_context(every))[1] == "Table 3 recommends every endpoint for XXXX."
    for marks in (None, {**every, "OS": ""}):
        (only,) = endpoint_note(note_context(marks))
        assert only.endswith("holds no copy of Liu et al. 2018 Table 3's use marks for XXXX.")
    with pytest.raises(ClinicalError, match="'OS': 'yes'"):
        endpoint_note(note_context({**every, "OS": "yes"}))
    with pytest.raises(ClinicalError, match="Explanation/Caution"):
        endpoint_note(note_context({**every, "OS": "caution"}))


def sex_context(sexes: list, variables: tuple[str, ...]):
    """A page context holding only what ``sex_flag`` reads."""
    from ctx.dataset import ClinicalSpec  # noqa: PLC0415

    spec = types.SimpleNamespace(clinical=ClinicalSpec(("G1",), "OS", {}, variables))
    clinical = types.SimpleNamespace(patients=pl.DataFrame({"sex": sexes}))
    return types.SimpleNamespace(spec=spec, clinical=clinical)


def test_the_sex_sentence_counts_the_few_where_sex_is_no_variable():
    """Only a cohort that does not declare sex flags the less frequent one, by its own count.
    Catches the sentence on a cohort of both sexes, the count of the frequent sex, and men
    counted where women are the few."""
    from ctx.clinical.page import sex_flag  # noqa: PLC0415

    few_men = ["female"] * 5 + ["male"] * 2 + [None]
    assert sex_flag(sex_context(few_men, ("stage", "age"))) == (
        "2 men are kept and flagged by sex.",)  # fmt: skip
    few_women = ["male"] * 4 + ["female"]
    assert sex_flag(sex_context(few_women, ("stage",))) == ("1 women are kept and flagged by sex.",)
    assert sex_flag(sex_context(few_men, ("stage", "sex"))) == ()
    assert sex_flag(sex_context(["female"] * 3, ("stage",))) == ()
    assert sex_flag(sex_context(["female", "male"], ("stage",))) == ()


def test_the_page_prints_the_cohorts_own_endpoint_and_no_sex_flag(clinical_root):
    """End to end: the made-up cohort declares no Table 3 marks and has women only, so the page
    says it holds no copy and prints no sex sentence. Catches the sentences built and dropped."""
    html = outputs(clinical_root, "overview_v1.html").read_text().replace("&#39;", "'")
    assert "The models below fit PFI, the dataset file's endpoint." in html
    assert "holds no copy of Liu et al. 2018 Table 3's use marks for FAKE-BRCA." in html
    assert "flagged by sex" not in html


def test_a_variable_with_nothing_missing_says_so():
    """The missing-values table's Reasons cell: a variable no patient lacks reads "nothing
    missing", one with gaps lists each reason most first, and reasons that do not count every
    gap stop the page. Catches the page's word for a number that does not exist ("not
    estimable") printed as a reason, and a gap with no reason passed in silence."""
    from ctx.clinical.page import reasons  # noqa: PLC0415

    patients = pl.DataFrame({
        "age": [61, 70, 55], "age_missing": [None, None, None],
        "stage": [1, None, None], "stage_missing": [None, "Stage X", "discrepancy"],
        "sex": ["female", None, None], "sex_missing": [None, "not reported", None],
    })  # fmt: skip
    assert reasons(patients, "age") == "nothing missing"
    assert reasons(patients, "stage") == "Stage X 1; discrepancy 1"
    with pytest.raises(ClinicalError, match="sex: 2 missing, reasons for 1"):
        reasons(patients, "sex")
