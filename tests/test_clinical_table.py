"""ctx.clinical.harmonise on a made-up GDC folder with a made-up CDR; and, opt-in with -m
realdata, the breast cohort against Liu et al. 2018 and the plan's numbers, typed here apart
from the code (lessons.md: a test that reads the constant it checks passes at any value of it)."""

import importlib.util
import os
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from gdc_fixture import aliquot, build, spec_file

from ctx.clinical import harmonise, survival
from ctx.dataset import ENDPOINTS, load
from ctx.gather import gdc, gdc_run

ROOT = Path(__file__).resolve().parents[1]
XLSX = importlib.util.find_spec("openpyxl") is not None
TUMOUR = "Primary Tumor"
THERAPY = "Hormone Therapy"
CDR_TEXT = {  # one CDR row as the workbook shows it, every value text
    "type": "BRCA", "age_at_initial_pathologic_diagnosis": "61", "gender": "FEMALE",
    "ajcc_pathologic_tumor_stage": "Stage IIA",
    "histological_type": "Infiltrating Ductal Carcinoma", "histological_grade": "[Not Available]",
    "menopause_status": ("Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior "
                         "hysterectomy)"),
    "Redaction": "", **{e: "0" for e in ENDPOINTS}, **{f"{e}.time": "730" for e in ENDPOINTS},
}  # fmt: skip


def cdr(rows: dict[str, dict[str, str]]) -> pl.DataFrame:
    """A CDR copy: one row per patient, the defaults above changed by each row's values."""
    return pl.DataFrame([{"bcr_patient_barcode": case, **CDR_TEXT, **changes}
                         for case, changes in rows.items()])  # fmt: skip


def tables(tmp_path: Path, **kw) -> gdc_run.Read:
    """Gather's tables from a made-up folder of five patients; patient 1 has two tumours."""
    items = [aliquot(1), aliquot(1, vial="B"), aliquot(2), aliquot(3), aliquot(4), aliquot(5),
             aliquot(5, code="11")]  # fmt: skip
    build(tmp_path / "gdc", items, **kw)
    return gdc_run.read_all(tmp_path / "gdc", load(spec_file(tmp_path / "x.toml")))


def made(found: gdc_run.Read, rows: dict[str, dict[str, str]], sheet: pl.DataFrame | None = None
         ) -> harmonise.Clinical:  # fmt: skip
    """The clinical table of the made-up folder with this CDR."""
    given = sheet if sheet is not None else found.sheet
    return harmonise.build(given, cdr(rows), found.cases, found.tests, found.treatments,
                           found.annotations, tumour=TUMOUR, therapy=THERAPY)  # fmt: skip


def case(n: int) -> str:
    """Made-up patient n's barcode."""
    return aliquot(n).case


def test_one_row_per_patient_through_chosen_whatever_the_sheet_holds(tmp_path):
    """The full sheet and the sheet cut to its chosen rows give the same table. Catches a
    consumer bypassing gdc.chosen (patient 1 twice, its second tumour counted) and a normal
    file taken as a tumour."""
    found = tables(tmp_path)
    rows = {case(n): {} for n in range(1, 6)}
    full, cut = made(found, rows), made(found, rows, gdc.chosen(found.sheet))
    assert full.patients.equals(cut.patients)
    assert full.patients["patient"].to_list() == [case(n) for n in range(1, 6)]
    assert full.patients["sample"][0] == aliquot(1).sample
    assert full.flow == (("tumour", 5), ("rule", 5), ("cdr", 5), ("redactions", 5))


def test_days_become_months_and_missing_tokens_become_reasons(tmp_path):
    """Catches days over 30 (or 30.5) for months, Stage X read as a stage, #N/A read as 0,
    and an unknown menopause kept as a value."""
    rows = {case(1): {"PFI.time": "1461"}, case(2): {"ajcc_pathologic_tumor_stage": "Stage X"},
            case(3): {"age_at_initial_pathologic_diagnosis": "#N/A",
                      "menopause_status": "[Unknown]"},
            case(4): {"ajcc_pathologic_tumor_stage": "Stage IIIC"},
            case(5): {"ajcc_pathologic_tumor_stage": "[Discrepancy]"}}  # fmt: skip
    p = made(tables(tmp_path), rows).patients
    assert p["PFI_months"][0] == 48.0  # four Julian years of 365.25 days
    assert p["stage"].to_list() == [2, None, 2, 3, None]
    assert harmonise.reasons(p, "stage") == [("Stage X", 1), ("discrepancy", 1)]
    assert p["age"][2] is None and harmonise.reasons(p, "age") == [("#N/A", 1)]
    assert p["menopause"].to_list() == ["post", "post", None, "post", "post"]


def test_a_flag_without_a_time_leaves_that_analysis_only(tmp_path):
    """Catches a flag without a time kept with a time of 0, and the patient dropped from every
    endpoint for one."""
    rows = {case(n): {} for n in range(1, 6)}
    rows[case(1)] = {"OS": "1", "OS.time": "#N/A"}
    rows[case(2)] = {"DFI": "#N/A", "DFI.time": "#N/A"}
    p = made(tables(tmp_path), rows).patients
    assert harmonise.analysed(p, "OS").height == 4 and harmonise.analysed(p, "PFI").height == 5
    assert harmonise.reasons(p, "OS") == [("flag without time", 1)]
    assert harmonise.reasons(p, "DFI") == [("no flag: #N/A", 1)]


def test_a_redaction_removes_the_patient_and_a_notification_does_not(tmp_path):
    """Catches the prototypes' 1,093: redacted cases kept (GS plan); and a notification taken
    for an exclusion."""
    note = {"case_submitter_id": case(3), "entity_type": "case",
            "entity_submitter_id": case(3), "classification": "Notification",
            "category": "Neoadjuvant therapy"}  # fmt: skip
    rows = {case(n): {} for n in range(1, 6)}
    rows[case(2)] = {"Redaction": "Redacted"}
    result = made(tables(tmp_path, annotations=[note]), rows)
    p = result.patients
    assert case(2) not in p["patient"].to_list() and p.height == 4
    assert p.filter(pl.col("neoadjuvant"))["patient"].to_list() == [case(3)]
    assert result.flow[-1] == ("redactions", 4)


def test_receptor_status_reads_every_record(tmp_path):
    """Patient 1 has ER positive then negative; patient 2 HER2 IHC equivocal and FISH
    positive; patient 3 hormone therapy recorded. Catches the first test winning, HER2 left
    equivocal where FISH decides, and an unrecorded therapy written as no."""
    ihc = {"molecular_analysis_method": "IHC"}
    extra = {
        case(1): {"tests": [{"gene_symbol": "ESR1", "test_result": r, **ihc}
                            for r in ("Positive", "Negative")]},
        case(2): {"tests": [{"gene_symbol": "ERBB2", "test_result": "Equivocal", **ihc},
                            {"gene_symbol": "ERBB2", "molecular_analysis_method": "FISH",
                             "test_result": "Positive"},
                            {"gene_symbol": "ESR1", "test_result": "Positive", **ihc}]},
        case(3): {"treatments": [{"treatment_type": "Hormone Therapy",
                                  "treatment_or_therapy": "yes"}]},
    }  # fmt: skip
    p = made(tables(tmp_path, cases=extra), {case(n): {} for n in range(1, 6)}).patients
    assert p["er"].to_list()[:2] == [None, "positive"]
    assert p["er_missing"].to_list()[:3] == ["conflicting", None, "not tested"]
    assert p["her2"][1] == "positive"
    assert p["therapy"].to_list() == [None, None, "yes", None, None]
    assert harmonise.reasons(p, "therapy") == [("not recorded", 4)]


def test_table_3_counts_every_flag_of_the_type():
    """Table 3 counts event flags over the CDR's rows, the patients without counts too."""
    rows = cdr({"P1": {"OS": "1"}, "P2": {"OS": "#N/A"}, "P3": {}})
    assert harmonise.flags(rows)["OS"] == (1, 1)


# Liu et al. 2018, Table 3, BRCA: events and censored per endpoint; Table 2: median follow-up;
# Figure 2A-D: ER-positive against ER-negative, log-rank p. The plan's own measurements on this
# Mac (section 1.2 and 2): the flow, the analysis rows, and our ER p-values.
TABLE_3 = {"OS": (151, 946), "PFI": (145, 952), "DFI": (84, 869), "DSS": (83, 995)}
ROWS = {"OS": (1087, 151), "PFI": (1087, 145), "DFI": (945, 84), "DSS": (1068, 83)}
FLOW = (("tumour", 1095), ("rule", 1093), ("cdr", 1092), ("redactions", 1088))
LIU_ER_P = {"OS": 0.097, "PFI": 0.005, "DFI": 0.001, "DSS": 0.009}
OUR_ER_P = {"OS": 0.086, "PFI": 0.0018, "DFI": 0.0004, "DSS": 0.0095}


@pytest.mark.realdata
@pytest.mark.skipif(not XLSX, reason="openpyxl ==3.1.5 is pinned in pixi.toml; run pixi install")
def test_tcga_brca_gives_liu_table_3_and_the_er_pattern():
    """Opt-in (-m realdata): the CDR's flags equal Liu's Table 3 in all 8 counts; the analysis
    rows, the flow and the median follow-up are the plan's; ER-negative is worse under all
    four endpoints, significant for PFI, DFI and DSS and not for OS, as in Liu's Figure 2.
    Catches days over 30, a redaction kept, a flag without time kept, and a second tumour."""
    folder = Path(os.environ.get("CTX_GDC_HOME", Path.home())) / "ctx_breast/data/raw/tcga_brca"
    if not folder.is_dir():
        pytest.skip(f"{folder} is not on this machine")
    found = gdc_run.read_all(folder, load(ROOT / "conf/datasets/tcga_brca.toml"))
    assert found.cdr is not None
    result = harmonise.build(found.sheet, found.cdr, found.cases, found.tests, found.treatments,
                             found.annotations, tumour=TUMOUR, therapy=THERAPY)  # fmt: skip
    assert dict(result.flags) == TABLE_3
    assert result.flow == FLOW
    p = result.patients
    for endpoint, (n, events) in ROWS.items():
        rows = harmonise.analysed(p, endpoint)
        assert (rows.height, int(rows[endpoint].sum())) == (n, events)
    days = [harmonise.whole(v)[0] for v in found.cdr["OS.time"]]
    months = np.median([d for d in days if d is not None]) / (365.25 / 12)
    assert round(months, 1) == 27.7
    for endpoint in ENDPOINTS:
        rows = harmonise.analysed(p, endpoint).filter(pl.col("er").is_not_null())
        t, e = rows[f"{endpoint}_days"].to_numpy(), rows[endpoint].to_numpy()
        _, _, p_value = survival.logrank(t, e, rows["er"].to_numpy())
        negative = (rows["er"] == "negative").to_numpy().astype(float)
        hr = survival.wald(survival.cox(t, e, negative[:, None])).hr[0]
        assert hr > 1
        assert float(f"{p_value:.2g}") == OUR_ER_P[endpoint]
        assert (p_value < 0.05) == (LIU_ER_P[endpoint] < 0.05)


# Liu et al. 2018, Table 3, read through PubMed Central's E-utilities (PMC6066282) on
# 2026-10-09: events and censored per endpoint for the kidney and lung cohorts of M10.10.
TABLE_3_M10_10 = {
    "tcga_kirc": {"OS": (177, 360), "PFI": (162, 375), "DFI": (15, 102), "DSS": (110, 415)},
    "tcga_luad": {"OS": (188, 334), "PFI": (213, 309), "DFI": (92, 217), "DSS": (116, 370)},
    "tcga_lusc": {"OS": (219, 285), "PFI": (149, 355), "DFI": (63, 241), "DSS": (91, 361)},
}
FOLDERS = {"tcga_kirc": "ctx_kidney", "tcga_luad": "ctx_lung", "tcga_lusc": "ctx_lung"}


@pytest.mark.realdata
@pytest.mark.skipif(not XLSX, reason="openpyxl ==3.1.5 is pinned in pixi.toml; run pixi install")
@pytest.mark.parametrize("slug", sorted(TABLE_3_M10_10))
def test_kidney_and_lung_give_liu_table_3(slug):
    """Opt-in (-m realdata), M10.10: the CDR's flags equal Liu's Table 3 in all 8 counts of
    each cohort, and the table holds one row per patient. Catches another cohort's CDR rows,
    a flag read as text, and a patient twice."""
    home = Path(os.environ.get("CTX_GDC_HOME", Path.home()))
    folder = home / FOLDERS[slug] / "data/raw" / slug
    if not folder.is_dir():
        pytest.skip(f"{folder} is not on this machine")
    spec = load(ROOT / f"conf/datasets/{slug}.toml")
    found = gdc_run.read_all(folder, spec)
    assert found.cdr is not None and spec.clinical is not None
    result = harmonise.build(found.sheet, found.cdr, found.cases, found.tests, found.treatments,
                             found.annotations, tumour=TUMOUR,
                             therapy=spec.clinical.therapy)  # fmt: skip
    assert dict(result.flags) == TABLE_3_M10_10[slug]
    assert result.patients["patient"].n_unique() == result.patients.height
