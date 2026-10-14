"""Package H: the ccRCC external ceilings, counted from the files on disk.

PDC000127's biospecimens, the TMT10 mapping workbook, the paper's protein table columns and
GDC's STAR-Counts files table, each counted on its own and then against each other.
"""

import csv
import io
import json
import tarfile
from collections import Counter

import openpyxl

RAW = "/Users/rytis/ctx_kidney/data/raw/"
PDC = RAW + "cptac_ccrcc/"
GDC = RAW + "cptac_ccrcc_star/"

bio = json.load(open(GDC + "pdc_biospecimen_PDC000127.json"))["data"]["biospecimenPerStudy"]
print("PDC biospecimens:", len(bio), "rows")
print("  aliquots", len({b["aliquot_submitter_id"] for b in bio}),
      "cases", len({b["case_submitter_id"] for b in bio}))
by_type = Counter(b["sample_type"] for b in bio)
print("  aliquots by sample_type", dict(by_type))
cases_type = {t: len({b["case_submitter_id"] for b in bio if b["sample_type"] == t}) for t in by_type}
print("  cases by sample_type", cases_type)
aliquot_type = {b["aliquot_submitter_id"]: b["sample_type"] for b in bio}
aliquot_case = {b["aliquot_submitter_id"]: b["case_submitter_id"] for b in bio}

book = openpyxl.load_workbook(PDC + "S044_TMT10_Label_to_Sample_Mapping_File_CCRCC_r2_Jan2019.xlsx",
                              read_only=True, data_only=True)
rows = list(book["TMT10_mapping_File"].iter_rows(values_only=True))
head = rows[0]
labels = [i for i, h in enumerate(head) if h and str(h).endswith("Specimen Label")]
plex_of = {}
proteome = [r for r in rows[1:] if r[1] and "_Proteome_" in r[1]]
for r in proteome:
    for i in labels:
        if r[i] and r[i] != "pooled sample":
            plex_of.setdefault(r[i], []).append(r[1])
print("Mapping workbook: proteome plexes", len(proteome), "labelled channels",
      sum(len(v) for v in plex_of.values()), "distinct labels", len(plex_of))
print("  labels in two plexes", sum(len(v) > 1 for v in plex_of.values()))

with tarfile.open(PDC + "Supplementary_Data_Proteome_TMT.tar.gz") as tar:
    member = next(m for m in tar.getmembers() if m.name.endswith(".tsv"))
    first = io.TextIOWrapper(tar.extractfile(member)).readline()
cols = next(csv.reader([first], delimiter="\t"))[4:]
print("Protein table:", member.name.split("/")[-1], len(cols), "sample columns,",
      len(set(cols)), "distinct")
print("  by PDC sample_type", dict(Counter(aliquot_type.get(c, "not in PDC") for c in cols)))
print("  not in the mapping workbook", [c for c in cols if c not in plex_of])
print("  mapping labels not in the table", [c for c in plex_of if c not in cols])
tum = {aliquot_case[c] for c in cols if aliquot_type.get(c) == "Primary Tumor"}
nat = {aliquot_case[c] for c in cols if aliquot_type.get(c) == "Solid Tissue Normal"}
print("  cases with a tumour column", len(tum), "with a normal column", len(nat),
      "with both", len(tum & nat))

files = list(csv.DictReader(open(GDC + "gdc_files_star_counts.tsv"), delimiter="\t"))
print("GDC files:", len(files))
types = Counter()
gcase = {"Primary Tumor": set(), "Solid Tissue Normal": set()}
for f in files:
    kinds = {f[k] for k in f if k.endswith("sample_type") and f[k]}
    types[" + ".join(sorted(kinds))] += 1
    for k in kinds:
        gcase.setdefault(k, set()).add(f["cases.0.submitter_id"])
print("  files by sample_type", dict(types))
print("  cases by sample_type", {k: len(v) for k, v in gcase.items()})
manifest = list(csv.DictReader(open(GDC + "gdc_manifest_star_counts.txt"), delimiter="\t"))
print("GDC manifest:", len(manifest), "files")
both = gcase["Primary Tumor"] & gcase["Solid Tissue Normal"]
print("  GDC cases with tumour and normal files", len(both))
print("  of them with protein tumour and normal", len(both & tum & nat))
