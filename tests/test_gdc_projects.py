"""Package S's gather rules on made-up GDC folders: two projects in one sheet, a case excluded by
the dataset file, and PDC's biospecimens; each test names the mutation it catches."""

import json
import re
from dataclasses import replace
from pathlib import Path

import polars as pl
import pytest
from clinical_root import ctx
from gdc_fixture import PROJECT, RELEASE, SPEC, aliquot, build

from ctx.dataset import DatasetError, load
from ctx.gather import gdc, gdc_run

OTHER = "FAKE-LUSC"


def two_projects(tmp_path: Path) -> Path:
    """A dataset file naming both made-up projects, their tumours compared."""
    text = SPEC.format(protocol="unstranded", cdr="").replace(
        f'project = "{PROJECT}"', f'projects = ["{PROJECT}", "{OTHER}"]'
    )
    text += (
        '\n[[contrasts]]\nname = "other_first"\nvariable = "project"\n'
        f'reference = "{PROJECT}"\ntested = "{OTHER}"\n'
        'filters = { condition = ["tumour"] }\nn = [3, 2]\n'
    )
    path = tmp_path / "two.toml"
    path.write_text(text)
    return path


def folders(tmp_path: Path) -> tuple[Path, Path]:
    """Two folders: three tumours and a normal of one project, two tumours of the other."""
    first = build(tmp_path / "a", [aliquot(1), aliquot(2), aliquot(3), aliquot(3, code="11")])
    second = build(tmp_path / "b", [aliquot(7), aliquot(8)], seed=5, project=OTHER)
    return first.path, second.path


def test_two_projects_make_one_sheet_with_a_project_column(tmp_path):
    """Catches: the second folder read with the first project's spec (its queries refused); the
    folders taken in the order given rather than the file's; run accessions repeated across the
    folders; the counts matrix out of step with the sheet."""
    a, b = folders(tmp_path)
    spec = load(two_projects(tmp_path))
    found = gdc_run.read_projects([b, a], spec)
    assert [one.project for one, _ in found] == [PROJECT, OTHER]
    combined = gdc_run.combine([r for _, r in found])
    sheet = combined.sheet
    assert sheet["project"].to_list() == [PROJECT] * 4 + [OTHER] * 2
    assert sheet["run_accession"].n_unique() == sheet.height == combined.counts.matrix.shape[1]
    assert combined.folder.files["sample"].to_list() == sheet["sample"].to_list()
    single = gdc_run.read_all(b, replace(spec, gdc=spec.gdc.one(OTHER)))
    assert (combined.counts.matrix[:, 4:] == single.counts.matrix).all()


def test_a_missing_or_doubled_project_is_refused(tmp_path):
    """Catches: one folder of a two-project dataset gathered alone; two folders of one project."""
    a, b = folders(tmp_path)
    spec = load(two_projects(tmp_path))
    with pytest.raises(gdc.GdcError, match="no --gdc-dir"):
        gdc_run.read_projects([a], spec)
    with pytest.raises(gdc.GdcError, match="both"):
        gdc_run.read_projects([a, a, b], spec)


def test_two_project_gather_writes_the_sheet_m5_fits_on(tmp_path):
    """End to end through the CLI with a repeatable --gdc-dir. Catches: typer taking only the
    last --gdc-dir; the project column dropped from the written sheet."""
    a, b = folders(tmp_path)
    root = tmp_path / "root"
    ctx(
        "--root",
        str(root),
        "gather",
        "--dataset",
        "fake_brca",
        "--gdc-dir",
        str(a),
        "--gdc-dir",
        str(b),
        "--dataset-file",
        str(two_projects(tmp_path)),
    )
    sheet = pl.read_parquet(next((root / "outputs" / "01_gather").glob("*_samplesheet_v1.parquet")))
    chosen = gdc.chosen(sheet, "Primary Tumor")
    assert chosen.group_by("project").len().sort("project").rows() == [(PROJECT, 3), (OTHER, 2)]
    # S.4: one download-checks table per project, each named by it; two under one title read
    # as one table printed twice. Catches: the project-blind title.
    (page,) = (root / "outputs" / "01_gather").glob("*_overview_v1.html")
    titles = re.findall(r"<caption>([^<]*)</caption>", page.read_text())
    checks = [t for t in titles if t.startswith("Download checks")]
    assert checks == [f"Download checks, {PROJECT}", f"Download checks, {OTHER}"]
    assert len(set(titles)) == len(titles)


def test_schema_refuses_project_and_projects_or_one_listed(tmp_path):
    """Catches: both keys accepted; a one-project list read as a combined dataset; a CDR declared
    for two projects (the CDR reader takes one type)."""
    text = two_projects(tmp_path).read_text()
    for bad, said in (
        (text.replace("[gdc]\n", f'[gdc]\nproject = "{PROJECT}"\n'), "not both"),
        (text.replace(f', "{OTHER}"]', "]"), "two or more"),
        (text.replace("[gdc]\n", '[gdc]\ncdr_file = "c.xlsx"\ncdr_type = "X"\n'), "one project"),
    ):
        path = tmp_path / "bad.toml"
        path.write_text(bad)
        with pytest.raises(DatasetError, match=said):
            load(path)


def test_an_excluded_case_drops_first_and_m5_never_sees_it(tmp_path):
    """Catches: the exclusion applied after the tie-break (the case keeps a chosen file); an
    exclusion without its reason accepted; a redaction reason overwriting the exclusion."""
    build(tmp_path / "gdc", [aliquot(1), aliquot(2), aliquot(2, code="11"), aliquot(2, vial="B")])
    case = aliquot(2).case
    text = SPEC.format(protocol="unstranded", cdr="").replace(
        "[guard]", f'exclude_cases = ["{case}"]\nexclude_reason = "made up."\n\n[guard]'
    )
    path = tmp_path / "x.toml"
    path.write_text(text)
    sheet = gdc_run.read_all(tmp_path / "gdc", load(path)).sheet
    mine = sheet.filter(pl.col("case") == case)
    assert mine["drop_reason"].to_list() == ["excluded"] * 3 and not mine["chosen"].any()
    assert gdc.chosen(sheet)["case"].to_list() == [aliquot(1).case]
    path.write_text(text.replace('exclude_reason = "made up."\n', ""))
    with pytest.raises(DatasetError, match="exclude_reason"):
        load(path)


def biospecimen_folder(tmp_path: Path, rows: list[dict[str, str]], logged: bool = True) -> Path:
    """A one-tumour folder with a PDC biospecimen record, logged or not."""
    folder = build(tmp_path / "gdc", [aliquot(1)]).path
    (folder / "pdc.json").write_text(json.dumps({"data": {"biospecimenPerStudy": rows}}))
    if logged:
        with (folder / gdc.RETRIEVAL).open("a") as log:
            log.write(f"2099-01-01T00:00:00Z\t{RELEASE}\tpdc.json\tmade up\n")
    return folder


ROW = {
    "aliquot_submitter_id": "CPTX1",
    "case_submitter_id": "FAKE-A1-0001",
    "sample_type": "Primary Tumor",
}


def test_biospecimens_read_and_checked(tmp_path):
    """Catches: an unlogged record read (no release behind it); an aliquot listed twice kept."""
    text = SPEC.format(protocol="unstranded", cdr="").replace(
        "[guard]", 'biospecimens = "pdc.json"\n\n[guard]'
    )
    (tmp_path / "x.toml").write_text(text)
    gdc_spec = load(tmp_path / "x.toml").gdc
    assert gdc_spec is not None
    good = biospecimen_folder(tmp_path / "a", [ROW, {**ROW, "aliquot_submitter_id": "CPTX0"}])
    found = gdc.biospecimens(good, gdc_spec)
    assert found["aliquot"].to_list() == ["CPTX0", "CPTX1"]
    with pytest.raises(gdc.GdcError, match="no row"):
        gdc.biospecimens(biospecimen_folder(tmp_path / "b", [ROW], logged=False), gdc_spec)
    with pytest.raises(gdc.GdcError, match="twice"):
        gdc.biospecimens(biospecimen_folder(tmp_path / "c", [ROW, ROW]), gdc_spec)
