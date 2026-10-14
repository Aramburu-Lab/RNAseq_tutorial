# CPTAC ccRCC: M7's rules, registered before any M7 run (package H, 2026-10-14)

Registered by this file's first commit, before M7 has run on CPTAC ccRCC: no protein number of this cohort has been computed in the course. They wait in a plan, not in `conf/datasets/cptac_ccrcc.toml`, because that file cannot declare `[proteo]` yet: the paper's protein table is a TSV of log2 abundances named by PDC aliquot, with no column annotation, and M7's reader takes only a GCT (`ctx.proteo.tables.read_gct`) of tumours (`ctx.proteo.run.across`). The block below goes under `[proteo.sanity]` unchanged once the reader and the schema take this table; lesson 17 applies: compare the rules in force with this commit's.

Not blind in one respect, stated: the sources' numbers were read on 2026-10-14 from Clark et al. 2019 Cell (doi:10.1016/j.cell.2019.10.007, PMC7331093) for this registration. Directions come from the sources; thresholds are M5's ALPHA 0.05 and LFC_MIN 1, Clark's own BH 0.01, or 0.5 for "a majority", none tuned on this cohort.

Clark analysed 103 clear cell tumours and 72 NATs for the mRNA-protein correlations; the 7 non-clear-cell tumours of the 110 were left out. A rule below that compares with those numbers must be run on the same 103, or say that it was not.

The metric names are proposals for the reader that will produce them; `class.<name>.p` in today's cross design tests "above", so OXPHOS gets its own "below" name.

```toml
[proteo.sanity]
registered_at = "2026-10-14 15:45"

[[proteo.sanity.rules]]
id = "C0"
text = "Recorded, not gating: the share of genes with a positive Spearman rho at BH 0.01 across tumours and across NATs, beside Clark 2019's 74 % and 52 %."
source = "Clark et al. 2019 Cell, doi:10.1016/j.cell.2019.10.007: 74 % and 52 % of 7,027 mRNA-protein pairs significant and positive, BH 0.01, over 103 tumours and 72 NATs"
gate = false
conditions = [{ metric = "pergene.tumour.positive_share", op = ">", value = 0 }]

[[proteo.sanity.rules]]
id = "C1"
text = "Most genes have a positive Spearman rho between mRNA and protein across tumours, at BH 0.01."
source = "Clark et al. 2019 Cell: 74 % of mRNA-protein pairs over 103 tumours; 0.5 is a course convention for a majority"
conditions = [{ metric = "pergene.tumour.positive_share", op = ">", value = 0.5 }]

[[proteo.sanity.rules]]
id = "C2"
text = "More genes correlate positively across tumours than across normal-adjacent tissue."
source = "Clark et al. 2019 Cell: 74 % over 103 tumours against 52 % over 72 NATs"
conditions = [{ metric = "pergene.tumour.positive_share", op = ">", other = "pergene.normal.positive_share" }]

[[proteo.sanity.rules]]
id = "C3"
text = "OXPHOS proteins fall in tumours against their paired normals, more than their mRNA does."
source = "Clark et al. 2019 Cell: OXPHOS downregulated in tumours at the protein level, its mRNA and protein uncoupled; MSigDB v2026.1.Hs HALLMARK_OXIDATIVE_PHOSPHORYLATION"
conditions = [
    { metric = "class.oxphos.protein_lfc.median", op = "<", value = 0 },
    { metric = "class.oxphos_below_other.p", op = "<", value = 0.05 },
]

[[proteo.sanity.rules]]
id = "C4"
text = "Glycolysis proteins rise in tumours against their paired normals."
source = "Clark et al. 2019 Cell: proteins of the glycolysis pathway and their mRNAs upregulated in tumours; MSigDB v2026.1.Hs HALLMARK_GLYCOLYSIS"
conditions = [{ metric = "class.glycolysis.protein_lfc.median", op = ">", value = 0 }]

[[proteo.sanity.rules]]
id = "C5"
text = "CA9 protein rises in clear cell tumours against their paired normals."
source = "Tostain et al. 2010 Eur J Cancer, doi:10.1016/j.ejca.2010.07.020: CA9 absent from healthy renal tissue, expressed in most clear cell carcinomas; thresholds M5's"
conditions = [
    { metric = "protein.tumour_normal.CA9.lfc", op = ">=", value = 1.0 },
    { metric = "protein.tumour_normal.CA9.padj", op = "<", value = 0.05 },
]

[[proteo.sanity.rules]]
id = "C6"
text = "The missingness figures come before every test on the page."
source = "TODO.md, Verification: M7 sanity, missingness plotted before any test"
conditions = [{ metric = "page.missingness_first", op = "==", value = 1 }]
```

## What the reader needs (package E owns the schema)

| Need | Where today |
|---|---|
| A protein table that is a TSV: gene symbol column `Index`, four descriptor columns, log2 abundance, samples named by aliquot | `src/ctx/proteo/tables.py:126` reads GCT 1.3 only |
| Each aliquot's sample type and case, from PDC's biospecimens; its plex, from the TMT10 mapping workbook | `src/ctx/dataset.py:264` ProteoSpec has `plex` as a GCT annotation name only |
| A design that compares tumour with paired normal at both layers | `src/ctx/dataset.py:27` PROTEO_DESIGNS and `:576` their keys; `src/ctx/proteo/run.py:112` pairs tumours only |
| Per-gene correlation across normals as well as tumours | `src/ctx/proteo/run.py:112` |
| The 103 clear cell tumours, from the clinical workbook | no field for a case exclusion list in a GDC dataset |
