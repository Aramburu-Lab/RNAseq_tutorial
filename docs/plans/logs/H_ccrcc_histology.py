# ruff: noqa: SIM115, PLR2004 - a proof record, run once and kept as run
"""Package H: the clinical workbook's histologic types.

And which cases have GDC tumour and normal files.
"""

import csv
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
print("cases", len(hist), dict(Counter(hist.values())))
files = list(
    csv.DictReader(open(RAW + "cptac_ccrcc_star/gdc_files_star_counts.tsv"), delimiter="\t")
)
kinds: dict[str, set[str]] = {}
for f in files:
    for k in f:
        if k.endswith("sample_type") and f[k]:
            kinds.setdefault(f["cases.0.submitter_id"], set()).add(f[k])
pairs = {c for c, k in kinds.items() if len(k) == 2}
print("GDC cases with both", len(pairs))
print("their histology", dict(Counter(hist.get(c, "not in workbook") for c in pairs)))
print(
    "non-clear-cell cases with both:",
    sorted(c for c in pairs if hist.get(c) != "Clear cell renal cell carcinoma"),
)
