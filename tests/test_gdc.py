"""ctx.gather.gdc on made-up GDC folders (tests/gdc_fixture.py), each test with the mutation it
catches; and, opt-in with -m realdata, the four TCGA projects on this Mac against the plan's
numbers (docs/plans/2026-10-07_m10_real_data_plan_proposal.md, section 1.1), typed here."""

import csv
import importlib.util
import io
import os
import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest
from gdc_fixture import PROJECT, aliquot, build, counts_text, spec_file

from ctx.dataset import DatasetError, load
from ctx.gather import gdc, gdc_run
from ctx.qc import registry
from ctx.report.builder import unpacked

ROOT = Path(__file__).resolve().parents[1]
CTX = Path(sys.executable).parent / "ctx"
XLSX = importlib.util.find_spec("openpyxl") is not None
NO_XLSX = "openpyxl ==3.1.5 is pinned in pixi.toml but not installed here; run pixi install"


def metric() -> registry.Metric:
    """The registry's strand metric, as gather reads it."""
    return registry.load(registry.resolve(None)).metrics[gdc.STRAND_METRIC]


def gathered(tmp_path: Path, items, protocol: str = "unstranded", **kw) -> gdc_run.Read:
    """Build a folder of these aliquots and read it as gather does."""
    build(tmp_path / "gdc", items, **kw)
    return gdc_run.read_all(tmp_path / "gdc", load(spec_file(tmp_path / "x.toml", protocol)))


def test_reads_a_folder_by_column_name_whatever_the_order_and_line_ends(tmp_path):
    """GDC fixes neither column order nor line ends: the fixture shuffles the metadata columns
    and ends lines with CRLF. Catches: a positional read (every field lands in another column),
    and a carriage return left on the last value (md5 then never matches the manifest)."""
    items = [aliquot(1), aliquot(2, code="11"), aliquot(3)]
    plain = gathered(tmp_path / "a", items, crlf=False, shuffle=False)
    shuffled = gathered(tmp_path / "b", items)
    assert plain.sheet.equals(shuffled.sheet)
    assert shuffled.sheet["sample"].to_list() == sorted(i.sample for i in items)
    assert shuffled.sheet["condition"].to_list() == ["tumour", "normal", "tumour"]
    assert not any("\r" in v for v in shuffled.sheet["md5sum"])


def test_another_gene_model_or_gene_order_or_par_y_count_is_refused(tmp_path):
    """Catches: a file of another GENCODE taken (the first line unchecked), two genes swapped in
    one file (misaligned rows), and a _PAR_Y row with counts dropped silently."""
    folder = build(tmp_path / "gdc", [aliquot(1), aliquot(2)])
    spec = load(spec_file(tmp_path / "x.toml"))
    second = folder.path / "00000001-0000-4000-8000-000000000000"
    path = next(second.glob("*.tsv"))
    good = path.read_text()
    rows = folder.counts[aliquot(2).sample]
    damaged = {
        "first line": counts_text(rows, model="GENCODE v38"),
        "gene order": counts_text([rows[1], rows[0], *rows[2:]]),
        "a _PAR_Y row": counts_text([*rows[:-1], (rows[-1][0], 3, 1, 2)]),
    }
    for label, body in damaged.items():
        path.write_text(body)
        with pytest.raises(gdc.GdcError, match=r"first line|gene order|_PAR_Y"):
            gdc.read_counts(gdc.read_folder(folder.path, spec.gdc), spec.gdc)
        assert label
    path.write_text(good)
    found = gdc.read_counts(gdc.read_folder(folder.path, spec.gdc), spec.gdc)
    # Stripped after the drop, the 12 IDs are unique; stripped before, the _PAR_Y copy repeats.
    assert found.genes.height == 12 and found.par_y == 1
    assert found.genes["gene_id"][-1] == "ENSG00000000012"


def test_library_size_sums_genes_only_and_takes_the_protocols_column(tmp_path):
    """Catches: the four N_ rows summed into library size (11 + 22 + 33 + 44 more), and a
    stranded cohort counted from the unstranded column."""
    items = [aliquot(1, strand="reverse"), aliquot(2, strand="reverse")]
    written = build(tmp_path / "gdc", items).counts
    found = gdc_run.read_all(tmp_path / "gdc", load(spec_file(tmp_path / "x.toml",
                                                              "stranded_second")))  # fmt: skip
    for j, item in enumerate(items):
        rows = written[item.sample][:-1]
        assert found.counts.matrix[:, j].tolist() == [s for _, _, _, s in rows]
        assert found.sheet["library_size"][j] == sum(s for _, _, _, s in rows)
    assert gdc.chosen(found.sheet).height == 2


def test_a_share_at_the_unstranded_edge_is_called_in_exact_arithmetic(tmp_path):
    """9 of 20 reads on the first strand sit on the band's edge, |2 x 0.45 - 1| = 0.1, which is
    not inside it (< 0.1): undetermined, dropped. One read more in 10**18 is inside it exactly,
    though its float rounds to 0.45. Catches the share called as a float (both undetermined)
    in place of the two sums as a Fraction."""
    near = (450_000_000_000_000_001, 549_999_999_999_999_999)
    found = gathered(tmp_path, [aliquot(1, strand=(9, 11)), aliquot(2, strand=near)])
    calls = found.sheet.select("stranded_first", "strand_call", "drop_reason").rows()
    assert calls == [(9, "undetermined", "protocol"), (near[0], "unstranded", "")]


def test_the_tie_break_runs_after_every_drop(tmp_path):
    """Case 1's FFPE aliquot and case 2's stranded aliquot each sort first. Catches: the first
    aliquot taken before the drops (case 1 and 2 then lose their tumour), and a drop applied
    after it."""
    items = [
        aliquot(1, vial="A", preservation="FFPE"), aliquot(1, vial="B"),
        aliquot(2, vial="A", strand="reverse"), aliquot(2, vial="B"), aliquot(2, vial="C"),
        aliquot(3, vial="A", preservation="FFPE"),
    ]  # fmt: skip
    sheet = gathered(tmp_path, items).sheet
    taken = gdc.chosen(sheet)
    assert sorted(taken["sample"]) == [aliquot(1, vial="B").sample, aliquot(2, vial="B").sample]
    reasons = dict(sheet.select("sample", "drop_reason").rows())
    assert reasons[aliquot(1).sample] == "ffpe" and reasons[aliquot(2).sample] == "protocol"
    assert reasons[aliquot(2, vial="C").sample] == "not_first"
    assert reasons[aliquot(3).sample] == "ffpe"


def test_a_redaction_drops_a_file_and_a_notification_drops_nothing(tmp_path):
    """Catches: a notification excluding a patient, and a redaction ignored."""
    a1, a2 = aliquot(1), aliquot(2)
    notes = [
        {"case_submitter_id": a1.case, "entity_type": "case", "entity_submitter_id": a1.case,
         "classification": "Notification", "category": "Neoadjuvant therapy"},
        {"case_submitter_id": a2.case, "entity_type": "aliquot", "entity_submitter_id": a2.sample,
         "classification": "Redaction", "category": "Genotype mismatch"},
        {"case_submitter_id": a1.case, "entity_type": "aggregated_somatic_mutation",
         "entity_submitter_id": "x", "classification": "Redaction", "category": "General"},
    ]  # fmt: skip
    found = gathered(tmp_path, [a1, a2], annotations=notes)
    assert gdc.chosen(found.sheet)["sample"].to_list() == [a1.sample]
    actions = dict(found.annotations.select("category", "action").rows())
    assert actions == {"Neoadjuvant therapy": "logged", "Genotype mismatch": "drops its files"}


def test_an_md5_mismatch_a_missing_file_or_another_release_stops_before_writing(tmp_path):
    """Catches: files never checksummed, a folder from another release read as this one, and a
    file the retrieval log does not name."""
    folder = build(tmp_path / "gdc", [aliquot(1), aliquot(2)])
    spec = load(spec_file(tmp_path / "x.toml"))
    path = next((folder.path / "00000000-0000-4000-8000-000000000000").glob("*.tsv"))
    path.write_text(path.read_text().replace("\t1.0\t", "\t2.0\t", 1))
    with pytest.raises(gdc.GdcError, match="md5 mismatch"):
        gdc.verify(gdc.read_folder(folder.path, spec.gdc))
    path.unlink()
    with pytest.raises(gdc.GdcError, match="missing"):
        gdc.verify(gdc.read_folder(folder.path, spec.gdc))
    status = folder.path / gdc.STATUS
    status.write_text(status.read_text().replace("99.0", "98.0"))
    with pytest.raises(gdc.GdcError, match=r"gdc_status.json: 'Data Release 98\.0"):
        gdc.read_folder(folder.path, spec.gdc)
    status.write_text(status.read_text().replace("98.0", "99.0"))
    log = folder.path / gdc.RETRIEVAL
    log.write_text(log.read_text().replace("99.0", "98.0", 1))
    with pytest.raises(gdc.GdcError, match=r"retrieved from 'Data Release 98\.0"):
        gdc.read_folder(folder.path, spec.gdc)
    log.write_text(log.read_text().replace("98.0", "99.0"))
    log.write_text("\n".join(line for line in log.read_text().splitlines()
                             if "00000001-" not in line) + "\n")  # fmt: skip
    with pytest.raises(gdc.GdcError, match="no row for 00000001-"):
        gdc.read_folder(folder.path, spec.gdc)


def test_a_file_of_two_cases_or_another_project_is_refused(tmp_path):
    """Catches: a second case's column ignored, and a query for another project read as this."""
    folder = build(tmp_path / "gdc", [aliquot(1)])
    spec = load(spec_file(tmp_path / "x.toml"))
    meta = folder.path / gdc.METADATA
    text = meta.read_text()
    meta.write_text(text.replace("cases.0.case_id", "cases.1.case_id", 1), newline="")
    with pytest.raises(gdc.GdcError, match="second entity"):
        gdc.read_folder(folder.path, spec.gdc)
    meta.write_text(text, newline="")
    query = folder.path / f"{gdc.METADATA}{gdc.QUERY}"
    query.write_text(query.read_text().replace(PROJECT, "OTHER-X"))
    with pytest.raises(gdc.GdcError, match="project_id"):
        gdc.read_folder(folder.path, spec.gdc)


def test_an_aliquot_listed_under_two_samples_of_its_case_is_one_file(tmp_path):
    """CPTAC-3 lists an aliquot cut from pooled tumour pieces once per piece (9 of 277 ccRCC
    files), and GDC's TSV gives every file the widest file's columns. Catches: a second entity
    refused by its header alone, and a second listing accepted when it names another aliquot or
    sample type."""
    folder = build(tmp_path / "gdc", [aliquot(1), aliquot(2)], crlf=False, shuffle=False)
    spec = load(spec_file(tmp_path / "x.toml"))
    meta = folder.path / gdc.METADATA
    rows = list(csv.DictReader(io.StringIO(meta.read_text()), delimiter="\t"))
    first, second = "cases.0.samples.0.", "cases.0.samples.1."
    fields = ("sample_type", "preservation_method", "submitter_id",
              "portions.0.analytes.0.aliquots.0.submitter_id")  # fmt: skip

    def listed(again: bool = True, **changed: str) -> pl.DataFrame:
        """Read the folder with the first file listed again under a second sample, or not."""
        for i, row in enumerate(rows):
            for f in fields:
                value = changed.get(f.split(".")[0], row[first + f])
                row[second + f] = value if again and i == 0 else ""
        out = io.StringIO()
        writer = csv.DictWriter(out, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        meta.write_text(out.getvalue(), newline="")
        return gdc.read_folder(folder.path, spec.gdc).files

    assert listed(again=False).height == 2
    piece = rows[0][first + "submitter_id"] + "-2"
    files = listed(submitter_id=piece)
    assert files.height == 2
    assert piece not in files["sample_id"].to_list()
    for changed in ({"portions": "CPT0000000001"}, {"sample_type": "Solid Tissue Normal"}):
        with pytest.raises(gdc.GdcError, match="second entity in 1 files"):
            listed(submitter_id=piece, **changed)


def test_barcode_parts_come_from_a_tcga_shaped_aliquot_of_the_same_case_only():
    """Catches: parts taken from an aliquot of another case, or from an identifier of another
    shape, where a CPTAC aliquot carries none."""
    a = aliquot(7, vial="C", plate="B123")
    assert gdc.barcode(a.sample, a.case) == {"vial": "C", "portion": "11", "plate": "B123",
                                             "site": "A1"}  # fmt: skip
    assert gdc.barcode(a.sample, "FAKE-A1-0008") == dict.fromkeys(("vial", "portion", "plate",
                                                                   "site"), "")  # fmt: skip
    assert set(gdc.barcode("CPT0001230001", "C3X-00001").values()) == {""}


def test_gather_writes_the_sheet_counts_and_clinical_tables(tmp_path):
    """The CLI path end to end on a made-up folder: one tumour per case through chosen(), the
    counts' columns in the sheet's order, a page. Catches: a counts column out of the sheet's
    order, and a table written without its README."""
    items = [aliquot(1), aliquot(1, vial="B"), aliquot(2, code="11"), aliquot(3)]
    tests = [{"gene_symbol": "ESR1", "molecular_analysis_method": "IHC", "test_result": "Positive"}]
    build(tmp_path / "gdc", items, cases={aliquot(1).case: {"tests": tests}})
    spec = spec_file(tmp_path / "fake_brca.toml")
    done = subprocess.run([str(CTX), "--root", str(tmp_path / "root"), "gather", "--dataset",
                           "fake_brca", "--gdc-dir", str(tmp_path / "gdc"), "--dataset-file",
                           str(spec)], capture_output=True, text=True, check=False)  # fmt: skip
    assert done.returncode == 0, done.stderr
    out = tmp_path / "root" / "outputs" / "01_gather"
    sheet = pl.read_parquet(next(out.glob("*_fake_brca_samplesheet_v1.parquet")))
    counts = pl.read_parquet(next(out.glob("*_fake_brca_counts_v1.parquet")))
    assert counts.columns[1:] == sheet["sample"].to_list()
    assert gdc.chosen(sheet, "Primary Tumor")["case"].to_list() == [aliquot(1).case,
                                                                    aliquot(3).case]  # fmt: skip
    tested = pl.read_parquet(next(out.glob("*_fake_brca_molecular_tests_v1.parquet")))
    assert tested.rows() == [(aliquot(1).case, "ESR1", "IHC", "Positive")]
    for table in out.glob("*.parquet"):
        assert table.with_name(table.name.replace(".parquet", "_README.md")).is_file()
    assert len(list(out.glob("*_overview_v1.html"))) == 1


def _page_words(tmp_path: Path, cdr: str | None) -> str:
    """The gather page of a made-up folder, with a CDR workbook named ``cdr`` or none."""
    from gdc_fixture import write_cdr  # noqa: PLC0415

    build(tmp_path / "gdc", [aliquot(1), aliquot(2)])
    if cdr is not None:
        columns = ["bcr_patient_barcode", "type", "OS", "OS.time"]
        write_cdr(tmp_path / "gdc" / cdr, [dict(zip(columns, (aliquot(1).case, "BRCA", 1, 365),
                                                     strict=True))], columns)  # fmt: skip
    spec = spec_file(tmp_path / "fake_brca.toml", cdr=cdr)
    done = subprocess.run([str(CTX), "--root", str(tmp_path / "root"), "gather", "--dataset",
                           "fake_brca", "--gdc-dir", str(tmp_path / "gdc"), "--dataset-file",
                           str(spec)], capture_output=True, text=True, check=False)  # fmt: skip
    assert done.returncode == 0, done.stderr
    (page,) = (tmp_path / "root" / "outputs" / "01_gather").glob("*_overview_v1.html")
    return unpacked(page.read_text())


@pytest.mark.skipif(not XLSX, reason=NO_XLSX)
def test_the_page_names_a_cdr_only_for_a_folder_that_holds_one(tmp_path):
    """H.4: a folder with no TCGA-CDR, as CPTAC's, is described without one, in Input and in
    Reproduce's inputs; a folder with one names it. Catches the input text chosen once for
    every folder, which named a workbook a CPTAC folder never holds."""
    without = _page_words(tmp_path / "cptac", None)
    assert "TCGA-CDR" not in without
    assert "counts, cases and annotations" in without
    held = _page_words(tmp_path / "tcga", "cdr.xlsx")
    assert "the TCGA-CDR workbook" in held and "annotations and TCGA-CDR" in held


@pytest.mark.skipif(not XLSX, reason=NO_XLSX)
def test_the_cdr_is_copied_verbatim_for_its_type_alone(tmp_path):
    """Catches: rows of another cancer type kept, a '#N/A' or a whole number reinterpreted, and
    the unnamed row-number column carried."""
    from gdc_fixture import write_cdr  # noqa: PLC0415

    path = tmp_path / "cdr.xlsx"
    columns = ["bcr_patient_barcode", "type", "OS", "OS.time"]
    rows = [dict(zip(columns, r, strict=True)) for r in (
        ("FAKE-A1-0001", "BRCA", 1, 365), ("FAKE-A1-0002", "BRCA", "#N/A", 12.0),
        ("FAKE-B1-0001", "OV", 0, 1))]  # fmt: skip
    write_cdr(path, rows, columns)
    frame = gdc.cdr(path, "BRCA")
    assert frame.columns == ["bcr_patient_barcode", "type", "OS", "OS.time"]
    assert frame.rows() == [("FAKE-A1-0001", "BRCA", "1", "365"),
                            ("FAKE-A1-0002", "BRCA", "#N/A", "12")]  # fmt: skip
    # openpyxl reads a whole number written without a point as an int; one stored as a float
    # (12.0) prints as the workbook shows it.
    assert gdc._cell(12.0) == "12" and gdc._cell(12.5) == "12.5" and gdc._cell(None) == ""


# The plan's table, section 1.1, measured 2026-10-07 on the folders on this Mac: per project,
# tumour files, cases, FFPE dropped, other protocol dropped, cases kept, cases lost, ties
# decided by the first aliquot ID, and tumour-normal pairs kept of those with both.
FOUR = {
    "tcga_brca": (1111, 1095, 13, 5, 1093, 2, 0, (113, 113)),
    "tcga_kirc": (541, 533, 4, 4, 533, 0, 0, (72, 72)),
    "tcga_luad": (540, 517, 12, 13, 515, 2, 0, (58, 58)),
    "tcga_lusc": (511, 501, 0, 9, 501, 0, 1, (51, 51)),
}
FOLDERS = {"tcga_brca": "ctx_breast", "tcga_kirc": "ctx_kidney", "tcga_luad": "ctx_lung",
           "tcga_lusc": "ctx_lung"}  # fmt: skip


@pytest.mark.realdata
@pytest.mark.skipif(not XLSX, reason=NO_XLSX)
@pytest.mark.parametrize("slug", sorted(FOUR))
def test_the_four_tcga_projects_give_the_plans_numbers(slug):
    """Opt-in (-m realdata): each project's folder on this Mac, read as gather reads it, gives
    the numbers the plan measured with its own prototype. Catches rows dropped or duplicated by
    the reader or the rule (lessons.md, a test at real size)."""
    folder = Path(os.environ.get("CTX_GDC_HOME", Path.home())) / FOLDERS[slug] / "data/raw" / slug
    if not folder.is_dir():
        pytest.skip(f"{folder} is not on this machine")
    found = gdc_run.read_all(folder, load(ROOT / "conf" / "datasets" / f"{slug}.toml"))
    sheet = found.sheet
    tumours = sheet.filter(pl.col("condition") == "tumour")
    reasons = tumours["drop_reason"].value_counts()
    count = dict(reasons.rows())
    taken = gdc.chosen(tumours)
    candidates = tumours.filter(pl.col("drop_reason").is_in(["", "not_first"]))
    ties = candidates.group_by("case").len().filter(pl.col("len") > 1).height
    pairs = set(gdc.chosen(sheet).filter(pl.col("condition") == "normal")["case"])
    both = set(sheet.filter(pl.col("condition") == "normal")["case"]) & set(tumours["case"])
    measured = (tumours.height, tumours["case"].n_unique(), count.get("ffpe", 0),
                count.get("protocol", 0), taken.height, tumours["case"].n_unique() - taken.height,
                ties, (len(pairs & set(taken["case"])), len(both)))  # fmt: skip
    assert measured == FOUR[slug]


def test_the_acknowledgement_names_the_project_and_quotes_its_own_policy(tmp_path):
    """GDC's FAQ asks for the GDC's citation and defers to each project's own policy (ui.toml
    [gdc], read 2026-10-09). Catches the project accession or release left out, TCGA's
    sentence printed for a project that does not declare it, and an unknown policy accepted."""
    plain = load(spec_file(tmp_path / "plain.toml")).gdc
    assert plain is not None
    (only,) = gdc.acknowledgement(plain)
    assert f"Data of {PROJECT} were retrieved" in only and plain.release in only
    assert "10.1038/s41588-021-00791-5" in only and "TCGA" not in only
    text = (tmp_path / "plain.toml").read_text()
    tcga = tmp_path / "tcga.toml"
    tcga.write_text(text.replace("[guard]", 'attribution = "tcga"\n[guard]', 1))
    found = load(tcga).gdc
    assert found is not None
    first, second = gdc.acknowledgement(found)
    assert first.startswith("The results shown here are in whole or part based upon data")
    assert first.endswith("TCGA Research Network: https://www.cancer.gov/tcga.")
    assert second == only
    tcga.write_text(text.replace("[guard]", 'attribution = "target"\n[guard]', 1))
    with pytest.raises(DatasetError, match="gdc.attribution"):
        load(tcga)
