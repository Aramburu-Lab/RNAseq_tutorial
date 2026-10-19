"""The clinical table: one row per patient, from gather's GDC tables (plan section 1.2).

The patient's tumour comes through ``gdc.chosen()``; endpoints, age, stage, grade, histology
and menopause from the TCGA-CDR; receptor status, sex and treatment from GDC's own records; the
notifications from gather's annotation log. Every missing value is null with a reason in a
``<column>_missing`` column, so "stage missing 24: 13 Stage X, 6 discrepancy, 5 not available"
is counted, not guessed. A CDR redaction removes the patient; a notification removes no one
(CONVENTIONS.md, Data and scope). Months are days over ``DAYS_PER_MONTH`` (the CDR's
notes sheet, row 23, gives days).
"""

import re
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import polars as pl

from ctx.clinical import ClinicalError
from ctx.dataset import ENDPOINTS
from ctx.gather import gdc

# Mean days per month, Julian year over twelve: the one conversion, mirrored by ctxSurv.
DAYS_PER_MONTH = 365.25 / 12
# A value that stands for no value, and the reason the table records for it.
MISSING = {
    "[Not Available]": "not available", "[Discrepancy]": "discrepancy",
    "[Unknown]": "unknown", "[Not Evaluated]": "not evaluated",
    "[Not Applicable]": "not applicable", "#N/A": "#N/A", "Stage X": "Stage X",
    "not reported": "not reported", "Unknown": "unknown", "": "empty",
}  # fmt: skip
STAGE = re.compile(r"Stage (IV|I{1,3})[A-C]?")
ROMAN = {"I": 1, "II": 2, "III": 3, "IV": 4}
# Histological grade G1 to G4; GX, grade not assessable, is a reason, never a grade.
GRADE = re.compile(r"G([1-4])")
MENOPAUSE = {"Pre": "pre", "Peri": "peri", "Post": "post", "Indeterminate": "indeterminate"}
RECEPTORS = {"er": "ESR1", "pr": "PGR"}
CALLED = ("Positive", "Negative")
EQUIVOCAL = "Equivocal"
# The notifications a patient's row flags, by GDC category.
NOTIFICATIONS = {
    "neoadjuvant": "Neoadjuvant therapy",
    "prior_malignancy": "Prior malignancy",
    "prior_treatment": ("History of unacceptable prior treatment related to a prior/other "
                        "malignancy"),
}  # fmt: skip
REDACTED = "Redacted"
# Liu et al. 2018 Table 3's use marks as a dataset file's [clinical.published.table3_use] copies
# them: a tick, a tick with an asterisk, a cross.
USES = ("recommended", "caution", "not_recommended")
CDR_COLUMNS = ("bcr_patient_barcode", "age_at_initial_pathologic_diagnosis",
               "ajcc_pathologic_tumor_stage", "histological_type", "histological_grade",
               "menopause_status",
               "Redaction", *ENDPOINTS, *(f"{e}.time" for e in ENDPOINTS))  # fmt: skip


@dataclass(frozen=True)
class Clinical:
    """The patient table and how it was reached.

    Attributes:
        patients: One row per analysed patient.
        flow: (step, patients) in order: tumour patients, after the rule, in the CDR, after
            the CDR's redactions.
        flags: Per endpoint, the CDR's event flags over every row of its type: (events,
            censored), what Liu et al. 2018 Table 3 counts.
    """

    patients: pl.DataFrame
    flow: tuple[tuple[str, int], ...]
    flags: Mapping[str, tuple[int, int]]


def missing(raw: str | None) -> str | None:
    """The reason a raw value is no value, or None for a real value."""
    return MISSING.get((raw or "").strip()) if raw is not None else "empty"


def stage(raw: str) -> tuple[int | None, str | None]:
    """AJCC stage group I to IV as 1 to 4; ``Stage X`` and the missing tokens are no stage."""
    if (reason := missing(raw)) is not None:
        return None, reason
    found = STAGE.fullmatch(raw.strip())
    return (ROMAN[found[1]], None) if found else (None, f"unread: {raw}")


def grade(raw: str) -> tuple[int | None, str | None]:
    """Histological grade G1 to G4 as 1 to 4; GX and the missing tokens are no grade."""
    if (reason := missing(raw)) is not None:
        return None, reason
    found = GRADE.fullmatch(raw.strip())
    return (int(found[1]), None) if found else (None, raw.strip())


def whole(raw: str) -> tuple[int | None, str | None]:
    """A whole number, or the reason it is none; ``#N/A`` is never 0."""
    if (reason := missing(raw)) is not None:
        return None, reason
    try:
        return int(raw), None
    except ValueError:
        return None, f"unread: {raw}"


def menopause(raw: str) -> tuple[str | None, str | None]:
    """pre, peri, post or indeterminate, from the CDR's long labels."""
    if (reason := missing(raw)) is not None:
        return None, reason
    first = raw.split(" ", maxsplit=1)[0]
    return (MENOPAUSE[first], None) if first in MENOPAUSE else (None, f"unread: {raw}")


def receptor(results: set[str]) -> tuple[str | None, str | None]:
    """One IHC call from every result a patient's records hold.

    One distinct call is the status; Equivocal alone, or two calls that disagree, are no
    status, with that reason. Every record is read: the first one does not win.
    """
    called = results & {*CALLED, EQUIVOCAL}
    if not called:
        return None, "not tested"
    if len(called) > 1:
        return None, "conflicting"
    value = called.pop()
    return (value.lower(), None) if value != EQUIVOCAL else (None, "equivocal")


def her2(ihc: set[str], fish: set[str]) -> tuple[str | None, str | None]:
    """HER2 by IHC, and by FISH where IHC is equivocal or absent (the rule set at build).

    IHC positive or negative decides; IHC that disagrees with itself is conflicting; then a
    single FISH call decides; otherwise the IHC reason stands (equivocal, or not tested).
    """
    status, reason = receptor(ihc)
    if status is not None or reason == "conflicting":
        return status, reason
    by_fish, fish_reason = receptor(fish & set(CALLED))
    if by_fish is not None:
        return by_fish, None
    return None, "conflicting" if fish_reason == "conflicting" else reason


def _tests(tests: pl.DataFrame) -> dict[str, dict[str, Any]]:
    """ER, PR and HER2 per patient, each with its reason when it is no call."""
    held: dict[tuple[str, str, str], set[str]] = {}
    for case, gene, method, result in tests.select("case", "gene_symbol", "method",
                                                   "result").rows():  # fmt: skip
        held.setdefault((case, gene, method), set()).add(result)
    out: dict[str, dict[str, Any]] = {}
    for case in {c for c, _, _ in held}:
        row: dict[str, Any] = {}
        for name, gene in RECEPTORS.items():
            row[name], row[f"{name}_missing"] = receptor(held.get((case, gene, "IHC"), set()))
        row["her2"], row["her2_missing"] = her2(
            held.get((case, "ERBB2", "IHC"), set()), held.get((case, "ERBB2", "FISH"), set())
        )
        out[case] = row
    return out


def _endpoint(row: Mapping[str, str], name: str) -> dict[str, Any]:
    """One endpoint's event, days and months, or null with the reason it leaves the analysis."""
    flag, flag_reason = whole(row[name])
    days, days_reason = whole(row[f"{name}.time"])
    reason = flag_reason and f"no flag: {flag_reason}"
    if flag is not None and (days is None or days < 0):
        reason = "flag without time"
    if reason or flag not in (0, 1):
        reason = reason or f"unread flag: {row[name]}"
        flag = days = None
    months = days / DAYS_PER_MONTH if days is not None else None
    return {name: flag, f"{name}_days": days, f"{name}_months": months, f"{name}_missing": reason}


def _cdr_row(row: Mapping[str, str]) -> dict[str, Any]:
    """Age, stage, histology, menopause and the four endpoints of one CDR row."""
    age, age_reason = whole(row["age_at_initial_pathologic_diagnosis"])
    group, stage_reason = stage(row["ajcc_pathologic_tumor_stage"])
    histology = row["histological_type"]
    meno, meno_reason = menopause(row["menopause_status"])
    graded, grade_reason = grade(row["histological_grade"])
    out = {
        "age": age, "age_missing": age_reason, "stage": group,
        "stage_label": None if stage_reason else row["ajcc_pathologic_tumor_stage"],
        "stage_missing": stage_reason, "grade": graded, "grade_missing": grade_reason,
        "histology": None if missing(histology) else histology,
        "histology_missing": missing(histology), "menopause": meno,
        "menopause_missing": meno_reason,
    }  # fmt: skip
    for name in ENDPOINTS:
        out |= _endpoint(row, name)
    return out


def flags(cdr: pl.DataFrame) -> dict[str, tuple[int, int]]:
    """(events, censored) of each endpoint's flag over every CDR row given: Table 3's count."""
    out = {}
    for name in ENDPOINTS:
        values = [whole(v)[0] for v in cdr[name]]
        out[name] = (values.count(1), values.count(0))
    return out


def build(  # noqa: PLR0913 - gather's six tables, each read once
    sheet: pl.DataFrame, cdr: pl.DataFrame, cases: pl.DataFrame, tests: pl.DataFrame,
    treatments: pl.DataFrame, annotations: pl.DataFrame, *, tumour: str, therapy: str,
) -> Clinical:  # fmt: skip
    """The patient table from gather's GDC tables.

    Args:
        sheet: Gather's sample sheet; the tumour of each patient is taken by ``gdc.chosen``.
        cdr: Gather's copy of the TCGA-CDR rows, as text.
        cases: Gather's GDC cases.
        tests: Gather's molecular tests.
        treatments: Gather's treatments.
        annotations: Gather's annotation log.
        tumour: The dataset file's tumour sample type.
        therapy: The GDC treatment type the ``therapy`` column records, e.g. Hormone Therapy.

    Returns:
        The table, the flow to it and the CDR's flags per endpoint.

    Raises:
        ClinicalError: If the CDR lacks a column the table needs, or holds a patient twice.
    """
    if absent := [c for c in CDR_COLUMNS if c not in cdr.columns]:
        raise ClinicalError(f"the TCGA-CDR copy has no {absent}")
    if cdr["bcr_patient_barcode"].n_unique() != cdr.height:
        raise ClinicalError("the TCGA-CDR copy holds a patient twice")
    taken = gdc.chosen(sheet, tumour)
    by_cdr = {r["bcr_patient_barcode"]: r for r in cdr.iter_rows(named=True)}
    in_cdr = taken.filter(pl.col("case").is_in(list(by_cdr)))
    redacted = sorted(c for c, r in by_cdr.items() if r["Redaction"] == REDACTED)
    kept = in_cdr.filter(~pl.col("case").is_in(redacted))
    receptors, sex = _tests(tests), dict(cases.select("case", "sex_at_birth").rows())
    given = (pl.col("treatment_type") == therapy) & (pl.col("treatment_or_therapy") == "yes")
    treated = set(treatments.filter(given)["case"])
    noted = {k: set(annotations.filter(pl.col("category") == v)["case"])
             for k, v in NOTIFICATIONS.items()}  # fmt: skip
    rows = []
    for case, sample, plate, site in kept.select("case", "sample", "plate", "site").rows():
        recorded = sex.get(case) or ""
        status = receptors.get(
            case,
            {
                f"{n}{s}": v
                for n in ("er", "pr", "her2")
                for s, v in (("", None), ("_missing", "not tested"))
            },
        )
        rows.append({
            "patient": case, "sample": sample, "plate": plate, "site": site,
            "sex": None if missing(recorded) else recorded, "sex_missing": missing(recorded),
            **_cdr_row(by_cdr[case]), **status,
            "therapy": "yes" if case in treated else None,
            "therapy_missing": None if case in treated else "not recorded",
            **{k: case in v for k, v in noted.items()},
        })  # fmt: skip
    # A column no row fills is typed Null, which compares with nothing: text, or whole numbers.
    patients = (pl.DataFrame(rows, infer_schema_length=None).sort("patient")
                .with_columns(pl.col(pl.Null).cast(pl.String))
                .with_columns(pl.col("age", "stage", "grade").cast(pl.Int64)))  # fmt: skip
    tumours = sheet.filter(pl.col("sample_type") == tumour)["case"].n_unique()
    flow = (("tumour", tumours), ("rule", taken.height), ("cdr", in_cdr.height),
            ("redactions", kept.height))  # fmt: skip
    return Clinical(patients, flow, flags(cdr))


def analysed(patients: pl.DataFrame, endpoint: str) -> pl.DataFrame:
    """The patients with an event flag and a time for this endpoint: its analysis rows."""
    return patients.filter(pl.col(endpoint).is_not_null())


def reasons(patients: pl.DataFrame, column: str) -> list[tuple[str, int]]:
    """(reason, patients) a column is missing for, most first, then by reason."""
    counted = patients[f"{column}_missing"].drop_nulls().value_counts()
    return sorted(counted.rows(), key=lambda r: (-r[1], r[0]))


def unread_uses(marks: Mapping[str, Any]) -> tuple[tuple[str, ...], str]:
    """The endpoints a copy of Table 3's Use marks leaves blank, then why M10 cannot read it.

    One reading for M10's page and for ``ctx run --check``: while any endpoint is blank the page
    says the dataset file holds no copy; a mark not in ``USES``, or a caution mark without the
    row's Explanation/Caution, stops the page.
    """
    blank = tuple(e for e in ENDPOINTS if not marks.get(e))
    if blank:
        return blank, ""
    if bad := {e: marks[e] for e in ENDPOINTS if marks[e] not in USES}:
        return (), f"{bad}, each one of {USES}"
    caution = [e for e in ENDPOINTS if marks[e] == "caution"]
    if caution and not marks.get("caution"):
        return (), f"copy Table 3's Explanation/Caution for {caution} as caution"
    return (), ""
