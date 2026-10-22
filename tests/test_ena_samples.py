"""scripts/python/ena_samples.py on fixture responses, through a stub opener; no network."""

import csv
import importlib.util
import io
import re
import sys
import urllib.error
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "python" / "ena_samples.py"
FIXTURES = ROOT / "tests" / "fixtures" / "ena"
STUDY = "PRJNA000001"


@pytest.fixture(scope="module")
def ena():
    spec = importlib.util.spec_from_file_location("ena_samples", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    # dataclasses resolves the module's annotations through sys.modules.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class Stub:
    """Answers the file report and the XML from the fixtures; ``failures`` 5xx first."""

    def __init__(self, failures=0, code=503, report=None):
        self.failures, self.code, self.urls = failures, code, []
        self.report = report if report is not None else (FIXTURES / "read_run.tsv").read_bytes()

    def __call__(self, url, timeout):
        self.urls.append(url)
        if self.failures:
            self.failures -= 1
            raise urllib.error.HTTPError(url, self.code, "fault", None, None)
        if "filereport" in url:
            return self.report
        return (FIXTURES / "samples.xml").read_bytes()


def read(path):
    return list(csv.DictReader(path.open(newline=""), delimiter="\t"))


def test_parse_keeps_every_tag_raw_and_in_order(ena):
    samples = ena.parse_sample_xml((FIXTURES / "samples.xml").read_bytes())
    assert [s.accession for s in samples] == ["SAMN00000001", "SAMN00000002", "SAMN00000003"]
    first = samples[0]
    assert (first.secondary, first.alias, first.title) == ("SRS0000001", "GSM0000001",
                                                           "P1-T (RNA-seq)")  # fmt: skip
    assert (first.scientific_name, first.broker_name) == ("Homo sapiens", "NCBI")
    assert first.center_name == "Example Lab, Example University"
    assert [t for t, _ in first.attributes] == ["recurrence", "patient", "tissue", "Sex",
                                                "ENA-FIRST-PUBLIC", "ENA-LAST-UPDATE"]  # fmt: skip
    assert dict(samples[2].attributes)["stage"] == "IA\t(pT1a)"
    assert dict(samples[2].attributes)["recurrence"] == "Progression (<2y)"


def test_wide_and_long_tables(ena, tmp_path):
    """Mutation caught: parse_sample_xml skipping its last SAMPLE_ATTRIBUTE drops
    ENA-LAST-UPDATE from every sample, and the header and long-table counts both fail."""
    stub = Stub()
    line = ena.run(STUDY, tmp_path, None, opener=stub)
    wide = read(tmp_path / f"{STUDY}_samples.tsv")
    assert list(wide[0]) == [*ena.FIXED, "recurrence", "patient", "tissue", "Sex",
                             "ENA-FIRST-PUBLIC", "ENA-LAST-UPDATE", "stage"]  # fmt: skip
    assert [r["runs"] for r in wide] == ["SRR0000011;SRR0000012", "SRR0000021", "SRR0000031"]
    assert {r["layout"] for r in wide} == {"PAIRED"}
    # A sample lacking a tag reads empty, never shifted into the next column.
    assert [r["recurrence"] for r in wide] == ["No Progression (>5y)", "", "Progression (<2y)"]
    assert [r["Sex"] for r in wide] == ["female", "female", ""]
    assert [r["stage"] for r in wide] == ["", "", "IA\t(pT1a)"]
    assert wide[2]["ENA-LAST-UPDATE"] == "2023-09-10"
    long = read(tmp_path / f"{STUDY}_attributes_long.tsv")
    assert len(long) == 6 + 5 + 6
    assert long[-1] == {"sample_accession": "SAMN00000003", "tag": "ENA-LAST-UPDATE",
                        "value": "2023-09-10"}  # fmt: skip
    text = (tmp_path / f"{STUDY}_samples.tsv").read_text()
    assert text.count('"') == 2  # only the value holding a tab is quoted
    assert line.startswith(f"{STUDY}: 3 samples, 4 runs, 2 requests (0 retries)")
    assert "recurrence 2, patient 3, tissue 3, Sex 2" in line and "stage 1" in line
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        f"{STUDY}_attributes_long.tsv", f"{STUDY}_samples.tsv"]  # fmt: skip


def test_xml_requests_batch(ena, monkeypatch, tmp_path):
    monkeypatch.setattr(ena, "BATCH", 2)
    stub = Stub()
    ena.run(STUDY, tmp_path, None, opener=stub)
    xml = [u.rsplit("/", 1)[1] for u in stub.urls if "browser" in u]
    assert xml == ["SAMN00000001,SAMN00000002", "SAMN00000003"]


def test_retries_a_5xx(ena, monkeypatch, tmp_path):
    monkeypatch.setattr(ena.time, "sleep", lambda _: None)
    stub = Stub(failures=2)
    line = ena.run(STUDY, tmp_path, None, opener=stub)
    assert "4 requests (2 retries)" in line
    stub = Stub(failures=3)
    with pytest.raises(ena.EnaError, match=rf"{STUDY}: no answer after 3 tries"):
        ena.run(STUDY, tmp_path / "x", None, opener=stub)


def test_a_4xx_is_not_retried(ena, tmp_path):
    stub = Stub(failures=1, code=400)
    with pytest.raises(ena.EnaError, match=f"{STUDY}: ENA refused"):
        ena.run(STUDY, tmp_path, None, opener=stub)
    assert len(stub.urls) == 1


def test_refuses_an_unresolved_or_malformed_accession(ena, tmp_path):
    header_only = b"run_accession\tsample_accession\tlibrary_layout\n"
    with pytest.raises(ena.EnaError, match=f"{STUDY}: ENA resolves no runs"):
        ena.run(STUDY, tmp_path, None, opener=Stub(report=header_only))
    with pytest.raises(ena.EnaError, match="SAMN00000001: not an ENA study"):
        ena.run("SAMN00000001", tmp_path, None, opener=Stub())
    assert list(tmp_path.iterdir()) == []


def test_a_missing_record_writes_nothing(ena, tmp_path):
    report = (FIXTURES / "read_run.tsv").read_bytes() + (
        b"SRR0000041\tSAMN00000004\tSRS0000004\tGSM0000004\tP3-T\tSRX0000041\tPAIRED\n"
    )
    with pytest.raises(ena.EnaError, match="no record for 1 samples"):
        ena.run(STUDY, tmp_path, None, opener=Stub(report=report))
    assert list(tmp_path.iterdir()) == []


def test_write_is_atomic(ena, tmp_path, monkeypatch):
    first, second = tmp_path / "a.tsv", tmp_path / "b.tsv"

    def fail(part, path):
        raise OSError("disk full")

    monkeypatch.setattr(ena.os, "replace", fail)
    with pytest.raises(OSError):
        ena.write_atomic({first: b"a", second: b"b"})
    assert list(tmp_path.iterdir()) == []


# labdocs lint's rule for a deliverable's stem, the one these names must pass.
DELIVERABLE_STEM = re.compile(r"^\d{4}-\d{2}-\d{2}_[a-z0-9]+(_[a-z0-9]+)*_v\d+$")


def test_deliverable_names(ena, tmp_path):
    """The accession alone by default, lowercased, as the deliverable rule allows no capital."""
    wide, long = ena.output_paths(STUDY, tmp_path, 1)
    today = ena.date.today().isoformat()
    assert wide.name == f"{today}_ena_samples_{STUDY.lower()}_v1.tsv"
    assert long.name == f"{today}_ena_samples_{STUDY.lower()}_attributes_long_v1.tsv"
    assert all(DELIVERABLE_STEM.match(p.stem) for p in (wide, long))
    assert ena.output_paths(STUDY, tmp_path, None)[0].name == f"{STUDY}_samples.tsv"


def test_a_label_names_the_files(ena, tmp_path, monkeypatch):
    """--label puts the label before the accession, through main, in both tables' names."""
    # main passes no opener, so run's default is the one to stub.
    monkeypatch.setattr(ena.run, "__defaults__", (Stub(), None))
    code = ena.main([STUDY, "--out-dir", str(tmp_path), "--version", "1",
                     "--label", f"{STUDY}=lung_tissue_luad"])  # fmt: skip
    assert code == 0
    today = ena.date.today().isoformat()
    stem = f"{today}_ena_samples_lung_tissue_luad_{STUDY.lower()}"
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        f"{stem}_attributes_long_v1.tsv", f"{stem}_v1.tsv"]  # fmt: skip
    assert all(DELIVERABLE_STEM.match(p.stem) for p in tmp_path.iterdir())
    assert ena.output_paths(STUDY, tmp_path, None, "blood")[1].name == (
        f"blood_{STUDY}_attributes_long.tsv")  # fmt: skip


@pytest.mark.parametrize("label", [
    "Lung_Tissue", "lung tissue", "lung-tissue", "lung_(luad)", "_lung", "lung__tissue", "",
])  # fmt: skip
def test_a_label_outside_the_naming_rule_is_refused(ena, tmp_path, label, capsys):
    out = tmp_path / "out"
    assert ena.main([STUDY, "--out-dir", str(out), "--label", f"{STUDY}={label}"]) == 1
    assert "lowercase words joined by underscores" in capsys.readouterr().err
    assert not out.exists()


@pytest.mark.parametrize("pair", ["lung", "PRJNA000002=lung"])
def test_a_label_names_a_study_given(ena, pair):
    with pytest.raises(ena.EnaError, match="<study>=<label>"):
        ena.parse_labels([pair], [STUDY])


def test_a_study_is_labelled_once(ena):
    with pytest.raises(ena.EnaError, match="twice"):
        ena.parse_labels([f"{STUDY}=a", f"{STUDY}=b"], [STUDY])


def test_tsv_quotes_only_when_it_must(ena):
    body = ena.tsv_bytes(["a", "b"], [["x y", "p\tq"]]).decode()
    assert body == 'a\tb\nx y\t"p\tq"\n'
    assert list(csv.reader(io.StringIO(body), delimiter="\t"))[1] == ["x y", "p\tq"]


def test_a_tag_named_as_a_fixed_column_keeps_both(ena):
    sample = ena.Sample(
        "SAMEA0000001",
        scientific_name="Homo sapiens",
        attributes=[("scientific_name", "Homo sapiens"), ("sex", "female")],
    )
    header, rows = ena.wide_rows([sample])
    assert header[len(ena.FIXED) :] == ["attribute_scientific_name", "sex"]
    assert len(set(header)) == len(header)
    assert ena.long_rows([sample])[0] == ["SAMEA0000001", "scientific_name", "Homo sapiens"]


def test_a_repeated_tag_is_joined_and_counted_once(ena):
    sample = ena.Sample(
        "SAMEA0000001", attributes=[("organism", "Homo sapiens"), ("organism", "Homo sapiens")]
    )
    header, rows = ena.wide_rows([sample])
    assert rows[0][header.index("organism")] == "Homo sapiens; Homo sapiens"
    assert "organism 1" in ena.summary("PRJEB000001", [sample], ena.Requests())
