# ruff: noqa: SIM115, T201 - a proof record, run once and kept as run
"""Package S: the clinical workbook's non-clear-cell cases, and what each layer holds of them.

The 7 cases go into cptac_ccrcc.toml's [gdc] exclude_cases; this prints them from the workbook
and counts, for the 110 tumour cases, the protein columns and GDC files by sample type.
"""

import csv
import io
import json
import tarfile
from collections import Counter

import openpyxl

RAW = "/Users/rytis/ctx_kidney/data/raw/"
book = openpyxl.load_workbook(
    RAW + "cptac_ccrcc/S050_S044_CPTAC_ccRCC_Discovery_Cohort_Clinical_Data_r4_Sept2019.xlsx",
    read_only=True,
    data_only=True,
)
rows = list(book["Patient_Clinical_Attributes"].iter_rows(values_only=True))
head = rows[0]
hist = {r[head.index("case_id")]: r[head.index("histologic_type")] for r in rows[1:]}
print("workbook cases", len(hist), dict(Counter(hist.values())))
other = sorted(c for c, h in hist.items() if h != "Clear cell renal cell carcinoma")
print("non-clear-cell", len(other), other)

bio = json.load(open(RAW + "cptac_ccrcc_star/pdc_biospecimen_PDC000127.json"))["data"][
    "biospecimenPerStudy"
]
kind = {b["aliquot_submitter_id"]: (b["case_submitter_id"], b["sample_type"]) for b in bio}
with tarfile.open(RAW + "cptac_ccrcc/Supplementary_Data_Proteome_TMT.tar.gz") as tar:
    member = next(m for m in tar.getmembers() if m.name.endswith(".tsv"))
    first = io.TextIOWrapper(tar.extractfile(member)).readline()
cols = next(csv.reader([first], delimiter="\t"))[4:]
for label in ("Primary Tumor", "Solid Tissue Normal"):
    cases = {kind[c][0] for c in cols if c in kind and kind[c][1] == label}
    print(label, "protein cases", len(cases), "of them clear cell", len(cases - set(other)))
files = list(
    csv.DictReader(open(RAW + "cptac_ccrcc_star/gdc_files_star_counts.tsv"), delimiter="\t")
)
gdc: dict[str, set[str]] = {}
for f in files:
    gdc.setdefault(f["cases.0.samples.0.sample_type"], set()).add(f["cases.0.submitter_id"])
for label, cases in sorted(gdc.items()):
    print(label, "GDC cases", len(cases), "of them clear cell", len(cases - set(other)))
both = gdc["Primary Tumor"] & gdc["Solid Tissue Normal"]
print("GDC cases with both", len(both), "clear cell", len(both - set(other)))
nat = {kind[c][0] for c in cols if c in kind and kind[c][1] == "Solid Tissue Normal"}
print(
    "clear-cell cases with GDC normal and protein normal",
    len(gdc["Solid Tissue Normal"] & nat - set(other)),
)
