# M8 Non-coding RNA on real data: plan, 2026-10-07

**The Step-1 half of M8 can be built now from files on this Mac. The TCGA half needs M10's GDC ingest, M5 run on TCGA dataset files, and five downloads of about 66 MB. Measured while planning: miR-210-3p rises in tumours of all four cohorts, by 2.78 to 4.82 log2, and anticorrelates with ISCU within 515 KIRC tumours, rho -0.236.**

Corrections to the brief, and one answer:

- Step-1's genes table holds `gene_id` and `gene_name` only, per its [README](/Users/rytis/ctx_tier1/outputs/01_gather/2026-10-07_gather_tier1_genes_v1_README.md). Biotypes come from the GENCODE 50 GTF through `gene_types`, [checkpoint.py:56](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/qc/checkpoint.py).
- TCGA is not all poly-A: GDC read groups list 10 rRNA-depleted RNA-seq libraries, 1 LUAD and 9 LUSC, and those 9 LUSC tumours also have a poly-T library. GDC API, Data Release 46.
- Step-1 is poly-A and unstranded: `library_selection` Oligo-dT on 210 of 210 runs, Salmon `expected_format` U on 210 of 210, in `~/ctx_tier1/data/raw/tier1_rnaseq/salmon/*/lib_format_counts.json`.

Every number below was measured by the scripts in `scratchpad/plan_m8/`. Line numbers are the working tree's at hand-back; other writers had uncommitted edits in `TODO.md` and `decisions.md`.

## 1. Data

Step-1 carries the lncRNA half because it is the course's own dataset; it has no small-RNA library, so no miRNA quantity comes from it, [CONVENTIONS.md:87](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/CONVENTIONS.md). TCGA carries the miRNA half, BRCA first, [TODO.md:692](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md).

| Input on this Mac | Size | Measured |
|---|---|---|
| `~/ctx_tier1/outputs/01_gather/2026-10-07_gather_tier1_counts_v1.parquet` | 10,852,713 B | 78,470 genes, all typed by the 78,941 gene rows of the GENCODE 50 GTF, 124,650,284 B |
| `~/ctx_tier1/outputs/05_de/2026-10-07_de_tier1_*_results_v1.parquet` | 10 contrasts | featured `ndn_gcsf`: 12,796 genes tested |
| `~/ctx_breast/data/raw/tcga_brca_mirna`, `tcga_brca_isoform` | 62 MB, 401 MB | 1,207 files, 1,079 cases: 1,096 tumour, 104 normal, 7 metastatic; 103 pairs |
| `~/ctx_kidney/data/raw/tcga_kirc_mirna`, `tcga_kirc_isoform` | 32 MB, 193 MB | 616 files: 544 tumour, 71 normal; 71 pairs |
| `~/ctx_lung/data/raw/tcga_luad_mirna`, `tcga_luad_isoform` | 29 MB, 215 MB | 567 files: 519 tumour, 46 normal; 46 pairs |
| `~/ctx_lung/data/raw/tcga_lusc_mirna`, `tcga_lusc_isoform` | 27 MB, 204 MB | 533 files: 488 tumour, 45 normal; 45 pairs |
| STAR counts in `tcga_brca`, `tcga_kirc`, `tcga_luad`, `tcga_lusc` | 4.9, 2.4, 2.4, 2.2 GB | GENCODE v36, 60,660 genes, 1,881 of them miRNA |

Every miRNA file lists the 1,881 hairpins of miRBase 21, the external ceiling; 1,502 to 1,604 carry a read in some file. Mature accessions from isoform rows: 2,102 to 2,241 per cohort. In BRCA 99.95 % of isoform reads are mature, 0.88 % of mature reads are cross-mapped, and a mature miRNA has a median 7 isoforms. Samples with both layers on one barcode: BRCA 1,198, KIRC 591, LUAD 535, LUSC 523; tumour-normal pairs with both: 103, 71, 19, 38.

| Download, needs Rytis's approval | Source | Bytes | Terms |
|---|---|---|---|
| `Predicted_Targets_Context_Scores.default_predictions.txt.zip` | [TargetScan 8.0](https://www.targetscan.org/vert_80/vert_80_data_download/) | 17,673,917 | none on the home or download page; cite McGeary 2019 |
| `miR_Family_Info.txt.zip`, also MIMAT to name | same | 162,636 | same |
| `Gene_info.txt.zip`, the prediction universe | same | 622,317 | same |
| `miRTarBase_SE_WR.csv`, reporter assay or Western blot | [miRTarBase](https://awi.cuhk.edu.cn/miRTarBase/downloads/files/10.0/miRTarBase_SE_WR.csv) | 3,348,778 | a 705 B LICENSE, a warranty disclaimer calling it "provided at no cost in the public domain"; cite Cui 2025, doi:10.1093/nar/gkae1072 |
| `gencode.v36.annotation.gtf.gz`, the STAR counts' annotation | [GDC](https://api.gdc.cancer.gov/data/be002a2c-3b27-43f3-9e0f-fd47db92a6b5) | 44.5 MB on GDC's page, bytes UNVERIFIED | open; md5 c03931958d4572148650d62eb6dec41a |
| Later, with M7: CPTAC-2 breast miRNA and isoform files, 104 tumour files each | GDC, open | 5,251,801 and 47,005,217 | GDC's |

Sizes are HTTP headers or GDC metadata, read 2026-10-07. Not recommended: `hsa_MTI.csv`, 337,103,345 B, which adds weak evidence, CLIP-seq among it. miRBase did not answer; RNAcentral's API named MIMAT0000267 hsa-miR-210-3p.

## 2. Analyses and figures

Contracts are the preview's `BarData`, `CurveData` and `ScatterData` with `Calls`, [preview_ncrna.py:169-356](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/scripts/python/preview_ncrna.py).

| Figure | Statistic | Control | Measured |
|---|---|---|---|
| F1 reads by biotype per sample; TCGA by sample and library type | descriptive, no test | subset by biotype, cell type, condition | Step-1 medians: protein-coding 88.3 %, lncRNA 8.93 %, Mt_rRNA 1.64 % with range 0.06 to 17.7 %, miRNA genes 0.0018 %. TCGA tumour medians by cohort: protein-coding 96.1 to 97.3 %, lncRNA 0.95 to 1.11 % |
| F2 expression by biotype against the filter | Mann-Whitney, r | min-count presets 1 to 50 | `ndn_gcsf` tests 2,687 of 34,948 lncRNAs |
| F3 volcano by biotype, M5's featured contrast | M5's p, BH within each biotype | Rytis's: M5's cut-offs, biotype colours, BH per biotype, [TODO.md:659](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md) | lncRNA calls pooled against within: 1,057 and 1,061 in `ndn_gcsf`, 66 and 38 in `mono_gcsf` |
| F4 lncRNA against its nearest protein-coding gene | Spearman on values and on residuals from cell type by condition means, BH over pairs; distant pairs as background | Rytis's drop-down, pairs named by both genes, [TODO.md:660](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md) | 2,682 pairs: 1,349 antisense overlaps, 354 same-strand, 411 within 10 kb, 568 farther; on the 26 `ndn_gcsf` samples, q < 0.05 on residuals for 20, 22, 12 and 8 % |
| F5 library type within one tumour, LUSC page | Wilcoxon signed-rank, n = 9 | none, two fixed groups | rRNA-depleted against poly-T: histone mRNA 0.526 and 0.041 %, small RNA 3.48 and 0.013 %, Mt_rRNA 0.017 and 0.630 %; 9 of 9 each, p = 0.0039 |
| F6 miRNA volcano, tumour against normal | `ctx.de.deseq.fit` unchanged, `~donor + group`, BH | M5's calls kind | tested, up, down: BRCA 350, 56, 69; KIRC 337, 38, 87; LUAD 354, 56, 47; LUSC 420, 75, 76 |
| F7 target shift: validated, predicted-only and non-target fold changes | KS and Mann-Whitney, median shift, BH over tests | target source; context++ presets precomputed | after download |
| F8 miRNA against target across samples | Spearman within tumours, BH over pairs | pair select over a precomputed table | miR-210-3p and ISCU: -0.450 over 586 KIRC samples, -0.236 within 515 tumours |
| Target enrichment among down genes | hypergeometric `ora`, [lite.py:66](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/function/lite.py); background the tested genes in TargetScan's universe | none | after download |

Sanity rules, registered in each dataset file before its M8 run, [TODO.md:740](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md):

| Rule | Direction from | Measured 2026-10-07 |
|---|---|---|
| miRNA genes hold under 0.1 % of counts in every poly-A library | course convention | Step-1 max 0.0142 %; TCGA max 0.049 % |
| Mature miRNAs absent from mRNA-seq: same-sample agreement low | size selection | KIRC, 200 samples, 1,246 hairpins placed by coordinates: median rho 0.045 |
| KIRC miR-210-3p up; ISCU anticorrelated within tumours | McCormick 2013, doi:10.1038/bjc.2013.56 | log2 FC 3.79, padj 6e-138; rho -0.236, p 6e-08 |
| KIRC miR-200c-3p and miR-141-3p down | Nakada 2008, doi:10.1002/path.2437 | -5.84, -6.02 |
| BRCA miR-21-5p and miR-145-5p | Iorio 2005, doi:10.1158/0008-5472.CAN-05-1783, whose abstract names them without a direction | +2.12, -2.25 |
| rRNA-depleted libraries carry more histone mRNA than poly-T ones | library chemistry | 9 of 9 |

These numbers were measured before any rule was registered. Directions come from the papers and thresholds from M5's `ALPHA` and `LFC_MIN`, [stats.py](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/de/stats.py); the record must say so.

Teaching notes from the data: Step-1 and TCGA lncRNA shares, 8.93 % and about 1 %, differ in three things at once: 34,948 GENCODE 50 lncRNAs against 16,901 in v36; Salmon splits overlapping reads where STAR drops them, a median 6.37 million per BRCA file; the libraries. HIF co-targets rise with miR-210-3p, CA9 at rho +0.363 within tumours, so correlation is not targeting.

## 3. The build

- Package `src/ctx/ncrna/`, [TODO.md:407](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md): `biotypes.py` with one grouping of at most six colours, importing `MT_RRNA` and `RRNA` from checkpoint; `neighbours.py`; `targets.py`; `run.py`; `page.py`. `adjust_within(p, strata, method)` sits beside `adjust`, [stats.py:65](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/de/stats.py); `de/` is free since the record, [decisions.md:241](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/knowledge/decisions.md). `gene_types` gains start, end and strand.
- CLI: `ctx ncrna --dataset --dataset-file` with the five planned flags, [ncrna.toml:46](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/report/content/ncrna.toml); outputs in `outputs/08_ncrna/`, [naming.py:21](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/io/naming.py).
- Inputs: gather outputs, M5 results, `annotation_gtf`, the ingest's miRNA tables, target resources through `ctx setup` with SHA-256 on arrival. The dataset file gains `[ncrna]`: featured contrast, miRNA contrast, featured miRNAs, pair count, rules.
- Content: `content/ncrna.toml` filled in its nine sections, `[preview]` removed at the release. Glossary entries: context++ score, ceRNA, cross-mapped read, RPM of miRNA-mapped reads, replication-dependent histone mRNA. The GDC acknowledgement on TCGA pages, [TODO.md:702](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md).
- Tests: a committed fixture of real GDC files cut to rows, with let-7a's three hairpins and miR-210, and synthetic counts with planted lncRNA pairs. Mutations to catch: one hairpin for the sum; precursor rows counted; pooled BH; values for residuals; aliquot as join key; validated and predicted merged; a flipped sign; the JS Spearman twin against Python under node.
- Budgets, Track B 10 min and 4 GB, [runtime.py:37](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/runtime.py): the GTF's gene rows 6.4 s at 170 MB by a scratch parser like `gene_types`; neighbours 0.9 s; isoform files aggregated per file 1.2 to 2.2 s at 276 MB at most, against 4.6 GB for BRCA when every row is held; miRNA fits 4.4 to 7.0 s at 405 MB plus a 257 MB child. The budget test measures the page against 15 MB at real size.

## 4. Dependencies and order

- From M10's ingest, without a second reader, per TCGA slug: the gather sheet, counts and genes with `gene_type`; the four `N_` rows per file; each STAR file's `library_selection` from its BAM's read groups, both libraries of a sample kept; the `mirna` and `isoform` kinds aggregated per file into mature and hairpin matrices with region and cross-mapped totals, the two totals checked equal per file, as they are over KIRC, 1,852,722,023 reads; joins on the sample barcode; a sortable `run_accession` for `rank_replicates`, [curation.py:41](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/gather/curation.py).
- M5 on Step-3, yes: M8 stratifies M5's TCGA results. M5's fit took 66.7 s at 2.32 GB for KIRC and 128.8 s at 3.57 GB for BRCA, in scratch processes that also held every count column; BRCA sits near the 4 GB cap.
- In parallel: the Step-1 half alongside M10; target readers on fixtures; downloads once approved. F3's control waits for GS.3's column table, [TODO.md:666](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md).

## 5. Decisions for Rytis

1. **BH within biotype over M5's tested genes, no refit.** Recommended; costs nothing new, and the filter stays M5's. A refit with per-biotype minimum counts costs one fit per contrast, 128.8 s and 3.57 GB for BRCA.
2. **Validated means miRTarBase strong evidence.** Recommended; 3,348,778 B. All human interactions take 337,103,345 B and hold CLIP-only pairs a filter must remove.
3. **The kit ships derived target rows, the raw files linked.** Recommended; students cannot re-derive targets offline. Shipping the raw files adds 21,807,648 B under terms TargetScan does not publish.
4. **Anticorrelation within tumours, normals apart.** Recommended; normals then form small strata, 20 in LUAD. Pooled, the tumour-normal gap drives rho, -0.450 against -0.236 for ISCU.
5. **F3's control after GS.3.** Recommended; F3 ships first without it. Building it now means a second column table, merged later.
6. **The library-type lesson from the 9 LUSC tumours.** Recommended; it costs one GDC metadata query in the ingest. Without it no sample shows both library types.
7. **Download GDC's v36 GTF.** Recommended; 44.5 MB. GENCODE 50 coordinates by stable ID lose 1,128 of 16,901 v36 lncRNAs.

## 6. Steps

1. [ ] Decisions and downloads approved; this plan copied into `docs/plans/`. S
2. [ ] `gene_types` returns start, end and strand; every caller's columns unchanged. S
3. [ ] `biotypes.py` and `adjust_within`, with mutations. M
4. [ ] `neighbours.py`, background pairs, the Spearman twin under node. M
5. [ ] `ctx ncrna`, `[ncrna]`, content, cards; the Step-1 page with F1, F2, F4, and F3 without its control. M
6. [ ] Step-1 run measured; browser at 1280 and 375 px in both themes; lessons. S
7. [ ] F3's calls control with biotype colours and per-biotype BH, after GS.3. M
8. [ ] With M10's owner: section 4's ingest requirements and the cut GDC fixture. M
9. [ ] Downloads through `ctx setup` with checksums. S
10. [ ] `targets.py`: TargetScan and miRTarBase readers, MIMAT names, target classes, with mutations. M
11. [ ] TCGA-BRCA page: F1, F3, F4 on v36, F6 to F8, enrichment, the sanity record. L
12. [ ] KIRC, LUAD and LUSC pages; F5 on LUSC. M
13. [ ] The miRNA, mRNA and protein view with M7, on CPTAC miRNA. L
14. [ ] Budget test at real size, kit rebuilt, TODO and lessons. S

## Not done

No download, no `pixi`, no page built; target analyses wait for the downloads. TargetScan publishes no terms I could find, and miRBase's were not read: one source sat behind a CAPTCHA, which was not passed. M10's ingest design was not seen. The scratch fits kept the first file per sample by aliquot name, and the mRNA fits dropped FFPE files: rules the ingest has not set.
