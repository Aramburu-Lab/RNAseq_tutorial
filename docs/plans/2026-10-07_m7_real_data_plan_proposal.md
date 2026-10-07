# M7 RNA versus protein on real data: plan, 2026-10-07

**Outcome.** M7 can run on real data: CPTAC prospective BRCA, protein from PDC000120 and RNA from GDC's open STAR-Counts through M10's ingest; for Step-1, Rieckmann 2017 copy numbers against a new paired NDN versus monocyte contrast in M5. 121 of the 125 tumour cases at PDC have open GDC RNA. No data file was downloaded. Krug 2020's median RNA-protein r is in an unread supplement, so that sanity rule waits.

## 1. Data

Why: CPTAC BRCA measured the same tumours at both layers and has a published reference; Step-1 has no proteome, and Rieckmann measured its two blood cell types in healthy donors, D11, [TODO.md:782](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md).

Metadata queries, 2026-10-07: [PDC GraphQL](https://pdc.cancer.gov/graphql), [GDC API](https://api.gdc.cancer.gov/status) at Data Release 46.0, [PRIDE](https://www.ebi.ac.uk/pride/ws/archive/v2/projects/PXD004352), [Zenodo](https://zenodo.org/api/records/8394329), [PyPI](https://pypi.org/pypi/cptac/json).

| Fact | Value | Source |
|---|---|---|
| Study | PDC000120 "Prospective Breast BI Proteome", TMT10, CPTAC2 Confirmatory | PDC `study` |
| Plex design | 17 TMT10 plexes, 25 fractions each; channel 131 a pooled reference in all 17, of 40 tumours: 9 triple negative, 12 HER2+, 19 ER+; 15 tumour plexes, 2 NAT-only plexes | PDC `studyExperimentalDesign`; Krug 2020 STAR Methods |
| Tumour channels | 133 channels, 125 aliquots of 125 case IDs; 7 cases in more than one plex; the paper reports 122 tumours plus 3 process replicates | PDC `biospecimenPerStudy`; [Krug 2020](https://doi.org/10.1016/j.cell.2020.10.036) |
| Normalisation applied | protein ratio = median PSM ratio to the reference, PSMs below 50 % precursor purity dropped, reporter ions at MS2 with no MS3 in the methods, so ratios compress; then a two-component Gaussian mixture per sample | Krug 2020, PMC8077737 |
| NATs | 18, excluded on QC | Krug 2020 |
| RNA at GDC | CPTAC-2 Breast: 133 open STAR-Counts, all Primary Tumor, 558,934,156 B; 121 match PDC tumour cases, 508,400,341 B; 4 PDC-only, 12 GDC-only; raw reads controlled, dbGaP phs000892 | GDC `files` |
| TCGA overlap | 0 of 139 case IDs at PDC and GDC are TCGA barcodes; the retrospective set, PDC000173, iTRAQ4, 109 cases, is TCGA tumours | PDC; [Mertins 2016](https://doi.org/10.1038/nature18003) |
| Subtypes | PAM50: 29 Basal, 57 LumA, 17 LumB, 14 HER2-E, 5 Normal-like; ER and/or PR+ 83 | Krug 2020 |
| Rieckmann donors | steady-state single-shot runs for Neutrophil 01-04, MO.classical 01-04, MO.intermediate and MO.nonclassical 01-04 | PRIDE file list |
| Copy numbers | Supplementary Table 6, "Protein copy numbers of immune cell types", by the proteomic ruler | [nature.com SI](https://doi.org/10.1038/ni.3693); [Wiśniewski 2014](https://doi.org/10.1074/mcp.M113.037309) |

The `cptac` package is not used: 1.5.14 is from 2024-06-06, it needs pandas and pyranges, against [CONVENTIONS.md:29](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/CONVENTIONS.md), and its downloader hard-codes a Zenodo token, [download_tools.py:16](https://github.com/PayneLab/cptac/blob/master/cptac/tools/download_tools.py). Its BRCA files are a pan-cancer reanalysis, not Krug 2020's. [cptac 0.9.7 brca.py:53-60](https://github.com/PayneLab/cptac/blob/v0.9.7/cptac/brca.py) read the paper's v5.4-public freeze as `proteome-ratio-norm-NArm.gct`, `rnaseq-fpkm-log2.gct` and a sample annotation with TMT.Plex, TMT.Channel, PAM50, ER.Updated.Clinical.Status and ERBB2.Gene.Amplified, brca.py:210-230; that the PDC tarball holds these is inferred from its version name, UNVERIFIED.

Downloads, each needing Rytis's permission:

| File | Source | Bytes | Terms | Kit |
|---|---|---|---|---|
| S060_Breast_Cancer_data_freeze_GCTfiles_v5.4-public.tar.gz, md5 1f108a42725a0b107e6f34e3d4d60c43 | PDC000120, CloudFront signed URL per request | 34,878,978 | [CPTAC DUA](https://cptac-data-portal.georgetown.edu/cptac/s?id=73217): acknowledgement; redistribution not addressed | derived tables only |
| 121 STAR-Counts | GDC, through M10's ingest | 508,400,341 | GDC open access, acknowledge GDC | as M10 ships TCGA |
| c5.go.cc.v2026.1.Hs.symbols.gmt | MSigDB 2026.1 | 778,106 | [CC BY 4.0](https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp) | yes |
| c2.cp.reactome.v2026.1.Hs.symbols.gmt | MSigDB 2026.1 | 898,680 | CC BY 4.0 | yes |
| Rieckmann Supplementary Table 6, XLSX | nature.com SI | 10,930 kB as displayed; URL UNVERIFIED | not open access in PubMed's record; reuse terms UNVERIFIED | no |
| Rieckmann Supplementary Text and Figures PDF | nature.com SI | 2,591 kB as displayed | as above | no |
| Krug 2020 supplement PDF, NIHMS1687926-supplement-1.pdf | PMC | UNVERIFIED | all rights reserved | no |

Needed only if a decision changes: the CDAP gene report, `CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv`, 61,431,484 B; PRIDE `total.proteomes.single.shot.zip`, 680,721,621 B, EBI terms of use; Rieckmann Table 8, RNA-seq, 2,557 kB. [immprot.org](http://www.immprot.org/) is personal and non-commercial use only.

## 2. Analyses and figures

Each reuses a preview contract, [preview_proteo.py:232-395](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/scripts/python/preview_proteo.py), and the controls in [the controls plan:148-155](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/docs/plans/2026-10-01_figure_controls_plan_approved.md).

| Figure | Contract | Statistics |
|---|---|---|
| F1 Protein groups quantified per tumour, by plex; control colours by plex or PAM50 | BarData | counts; PAM50 by plex: chi-square permutation p, Cramér's V |
| F2 Quantified share per protein, one curve per RNA-abundance tertile, the 70 % line | CurveData | descriptive, drawn before any test |
| F3 Per-gene RNA-protein correlation by class | CurveData | Spearman, Pearson as a control; BH; Kruskal-Wallis, Dunn with BH, ε² |
| F4 RNA against protein fold change: Basal vs LumA, ER- vs ER+ | ScatterData, 50/90 % contours | M5's DESeq2 calls; protein Welch t with BH, called by [stats.py:33](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/de/stats.py) `classify`; Spearman with bootstrap CI |
| F5 Paired RNA and protein, tumour by tumour | HeatmapData, tracks PAM50 and plex | per-gene ρ; no further test |
| F6 ERBB2, GRB7, STARD3: RNA against protein per tumour, coloured by ERBB2 amplification | ScatterData | per-gene Spearman |
| F7 Step-1: NDN over monocyte, RNA fold change against copy-number fold change, granule genes marked | ScatterData, contours | Spearman with bootstrap CI; granule against other genes on protein minus RNA fold change: Mann-Whitney, rank-biserial r |

Rules: missing is never zero; tests use proteins quantified in at least 70 % of samples, decided in exact arithmetic, [CONVENTIONS.md:86](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/CONVENTIONS.md); the page states the plex design; F7 is labelled cross-study, [CONVENTIONS.md:88](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/CONVENTIONS.md). One row per gene before the filter, then correlation, then BH, with a test of that order. Classes: GO:CC cytosolic ribosome 121 genes, spliceosomal complex 219, RNA polymerase complex 117, proteasome 66; Reactome amino acid metabolism 358, fatty acid metabolism 175; granules: azurophil 157, specific 160, tertiary 164, counted on 2026-10-07 from MSigDB's gene-set JSON.

Sanity rules, registered before any protein value is read, [TODO.md:724](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md):

| Rule | Reference |
|---|---|
| Median per-gene Pearson r within an agreed band of Krug 2020's, over genes in both layers, our count beside their 8,362 and the 10,107 proteins | Krug 2020 STAR Methods; median UNVERIFIED |
| ESR1, PGR, ERBB2, GATA3 each r above 0.5 | Mertins 2016, 77 TCGA tumours, iTRAQ: 0.74, 0.74, 0.84, 0.83; median Pearson 0.39; 6,135 of 9,302 pairs at FDR < 0.05 |
| Complex subunits below metabolic enzymes, one-sided Mann-Whitney | Mertins 2016: metabolism enriched among positive, ribosome, RNA polymerase and splicing among negative |
| F1 and F2 precede every test on the page | TODO.md:724 |
| Step-1: ρ above 0 with its CI; granule genes above the others | [Rørvig 2013](https://doi.org/10.1189/jlb.1212619); [Cowland and Borregaard 2016](https://doi.org/10.1111/imr.12440) |

Rieckmann's own RNA-protein Pearson values are in Supplementary Figure 6, not yet read.

## 3. The build

- `src/ctx/proteo/`: `tables.py` reads GCT members from the tarball with `tarfile` and the XLSX; `measure.py` holds the statistics; `pdc.py` resolves a fresh signed URL for `ctx setup`; `figures.py`, `page.py`, `run.py`. No `cptac.py` and no second `concordance.py`, unlike [TODO.md:406](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md).
- `ctx.concordance._spearman` becomes public for signed values; `agreement()` calls it, D-040's tests unchanged, [concordance.py:174](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/concordance.py). log2(CPM + 1) gets one definition, shared with M10; today CPM is inline at [filters.py:41](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/de/filters.py).
- Inputs: the dataset's gather counts and sheet, M5's results for its contrasts, and the protein table and GMTs as `ctx setup` resources with SHA-256 measured on arrival, [dataset.py:27](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/dataset.py).
- CLI: `ctx proteo --dataset <slug> --method spearman --min-quantified 0.7 --lfc-min 1 --gene-class all`. The dataset file's `[proteo]` block chooses across-samples or cross-dataset; PDC and PXD identifiers stay in `conf/datasets/`, [test_repo.py:26](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/tests/test_repo.py).
- Content: `proteo.toml` rewritten for real data, "Made-up data." moved into the preview builder; cards; glossary gains label-free quantification and proteomic ruler, and `tmt` loses its wrong LFQ aliases, [glossary.yaml:1743](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/report/glossary.yaml).
- Fixture tests on planted data, each naming its mutation, under `PYTHONDONTWRITEBYTECODE=1`: zero imputation, the 70 % edge both ways, pairwise-complete ranks, BH over tested genes only, the order selftest, duplicate isoforms, a planted complex class, compressed ratios.
- Budgets, Track B under 10 min and 4 GB, page under 15 MB. Measured on synthetic data at the paper's sizes: per-gene Spearman with BH over 3,534 kept genes of 10,107 took 0.134 s, equal to scipy within 1e-12 on 200 genes; a 1,000-resample bootstrap over 8,362 genes 1.284 s; the process 2.39 s, peak 864,026,624 B. The made-up page is 6,115,268 B; the real page is measured at real size.

## 4. Dependencies and order

- CPTAC needs M10's GDC ingest. `~/ctx_fetch_gdc.sh`, [TODO.md:2594](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md), fetches TCGA projects only; CPTAC-2 needs a primary-site filter, no CDR step, and a strand column chosen from the data, since the libraries are TruSeq Stranded.
- PAM50 and ER reach M5's sheet through a declared annotation join, as Step-1's SDRF join, [tier1.toml:60-63](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/conf/datasets/tier1.toml).
- Step-1 needs a tier1 contrast: steady-state blood has 9 donors with both NDN and monocytes, each donor in one read-length batch, 7 short and 2 long.
- Parallel now: the concordance refactor, readers, statistics, fixtures, content, glossary, Step-1.
- For the M8 planner: GDC holds 104 open CPTAC-2 breast miRNA files, the tri-layer's own miRNA.

## 5. Decisions for Rytis

1. **CPTAC sources: RNA from GDC STAR-Counts, protein and annotation from the paper's freeze.** Cost: 543 MB and waiting for M10. Freeze alone: 35 MB and the paper's numbers, but its RNA is hg19 FPKM, not counts, so M5 cannot test it. Zenodo reanalysis: fixed URLs and a CC BY record, but not the paper's numbers, and PAM50 in it UNVERIFIED.
2. **Approve the downloads in section 1.** Without the supplements the published-median rule cannot be registered.
3. **Kit ships derived tables and figure data only.** Shipping the freeze rests on a DUA silent on redistribution; Rieckmann's table is publisher-owned.
4. **Protein calls use D-012's rule, padj < 0.05 and |log2 FC| at least 1.** One meaning of "changed"; compression then hides protein calls, which F4's control shows. 0.5 for protein keeps more calls, two definitions.
5. **Step-1: contrast `ndn_vs_mono_steady` in tier1.toml**, variable cell_type, `~donor + group`, paired, steady-state blood, n [9, 9] by the one-sample-per-individual rule, [tier1.toml:129](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/conf/datasets/tier1.toml), with Monocyte matched to MO.classical. Cost: one more M5 contrast, and the 92 recorded metrics checked again, D-046. The alternative is a second DE engine in M7.
6. **Pin fastexcel 0.21.0**, MIT, osx-arm64 and linux-64, for the XLSX. Cost: a shared-environment change. Alternative: a stdlib XLSX reader.

Also: drop `--evidence`, which changes no test.

## 6. Steps

- [ ] P7.0 S: downloads approved; Krug's supplement and Rieckmann's Figure 6 read; values recorded with sources.
- [ ] P7.1 S: sanity rules registered in the dataset files and committed before any protein value is read.
- [ ] P7.2 M: public Spearman interval in `ctx.concordance`; one CPM definition with M10; mutation tests.
- [ ] P7.3 M: `proteo.tables`, one row per gene with a one-to-many report.
- [ ] P7.4 M: `proteo.measure` with the order selftest and mutation sweep.
- [ ] P7.5 S: `proteo.pdc` signed-URL resolver into `setup._fetch_one`, md5 then SHA-256.
- [ ] P7.6 S: tier1 contrast; `ctx de` again; 92 metrics unchanged.
- [ ] P7.7 L: `ctx proteo` page, content, cards, glossary, end-to-end fixture test.
- [ ] P7.8 M: Step-1 on real data; sanity; browser at 1280 and 375 px, both themes.
- [ ] P7.9 M, after M10's ingest: CPTAC BRCA through gather, M5 and M7; sanity, and the same median on the freeze's own RNA with the GDC copy made unreachable; budget measured.
- [ ] P7.10 M: D-051 controls F1 to F4 after re-measuring their teaching value.
- [ ] P7.11 S: real-size budget test, kit rebuilt and checked, lessons, review.

## Not done or not verified

The freeze's contents and what "NArm" removed; which 3 of 125 tumours are replicates; Krug's median r; Rieckmann's Figure 6 values, table layout, file URLs and whether donors 01-04 are the same people across cell types; CPTAC library strandedness in the counts; PDC000173's files, whose listing returned empty fields.
