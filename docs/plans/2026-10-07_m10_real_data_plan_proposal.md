# M10 Clinical and survival, breast first, with the GDC count ingest

**Outcome.** Breast needs no download: counts, GDC clinical, annotations and the TCGA-CDR are on this Mac and pass every check below, measured 2026-10-07 by scripts in the scratchpad's `plan_m10/`. 1,088 patients enter M10, PFI 1,087 with 145 events. The CDR's event flags equal Liu et al. 2018 Table 3 in all 16 counts of BRCA, KIRC, LUAD and LUSC, and the ER exercise gives the paper's pattern. The one GDC reader is designed for TCGA and CPTAC alike; three decisions wait in section 5.

**Corrections to recorded text, measured**

| Finding | Measured | Corrects |
|---|---|---|
| 52 of 3,008 TCGA counts files come from stranded libraries | first-strand share 0.0444 to 0.1243, the rest 0.4948 to 0.5332 | the first-aliquot rule [gene_selection_plan.md:143](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/docs/plans/2026-09-30_gene_selection_plan.md) avoids them in BRCA by sort order only |
| 2 BRCA patients have only FFPE tumour files | both censored | D-049's 1,089 PFI patients becomes 1,087; 145 events stand [decisions.md:112](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/knowledge/decisions.md) |
| Data Release 46.0 has no `demographic.gender` | the API warns "unrecognized values"; `sex_at_birth` matches the CDR, 1,085 female, 12 male | [TODO.md:284](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md) |
| CDR times are days | the file's own TCGA-CDR_Notes, row 23 | [TODO.md:280](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md), [decisions.md:1183](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/knowledge/decisions.md) |

## 1. Data

### 1.1 The GDC ingest, one reader

`ctx.gather.gdc`, run as `ctx gather --dataset tcga_brca --gdc-dir data/raw/tcga_brca --dataset-file conf/datasets/tcga_brca.toml`, is the only code that opens a GDC folder, as gather alone reads fetchngs output [CONVENTIONS.md:30](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/CONVENTIONS.md); M7, M8, M10 and Step-3 DE read `outputs/01_gather/` [TODO.md:171](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md). It works from GDC's metadata fields, never from barcode shape, so TCGA and CPTAC pass alike; a TCGA barcode adds vial, portion and plate. The dataset file names project, release, workflow, primary-site filter, library protocol and clinical source [CONVENTIONS.md:96](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/CONVENTIONS.md). CPTAC-2 breast, queried today and not downloaded: 133 open STAR-Counts files, 133 cases, all primary tumour and frozen, no CDR.

Checks, before any write:
- manifest; metadata TSV by column name, CRLF stripped; release in `gdc_status.json`; retrieval log; primary site where declared, a field CPTAC's fetch must request. The miRNA and isoform folders share the 11 metadata columns, so M8 adds parsers only.
- md5: 1,231 of 1,231 BRCA files match, in 7.7 s.
- all 3,008 TCGA files: `# gene-model: GENCODE v36`, one gene order, four `N_` rows, 60,660 rows. The 44 `_PAR_Y` rows are zero everywhere; without them the 60,616 unversioned IDs, D-049's gene identity, are unique.
- strand share, `stranded_first / (stranded_first + stranded_second)` over genes, classed by nf-core's 0.8 and 0.1 thresholds in a registry entry beside `gene_strand_bias` [qc_metrics.toml:396](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/conf/qc_metrics.toml); no file falls between classes. Files of a protocol other than the declared one are dropped, and counts come from the protocol's column: `unstranded` for TCGA, the matching stranded column for CPTAC-2. Library size leaves out the `N_` rows.

One tumour per case, in this order, each drop logged:

| Step | BRCA | KIRC | LUAD | LUSC |
|---|---|---|---|---|
| Tumour files; cases | 1,111; 1,095 | 541; 533 | 540; 517 | 511; 501 |
| FFPE dropped | 13 | 4 | 12 | 0 |
| Other protocol dropped | 5 | 4 | 13 | 9 |
| Cases kept; lost | 1,093; 2 | 533; 0 | 515; 2 | 501; 0 |
| Ties to the first aliquot ID | 0 | 0 | 0 | 1 |
| Tumour-normal pairs kept | 113 of 113 | 72 of 72 | 58 of 58 | 51 of 51 |

In BRCA, `N_noFeature` takes a median 61.8 % of FFPE files' counted plus no-feature reads, against 4.9 % in frozen unstranded files.

Outputs, in gather's pattern [harmonise.py:374](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/gather/harmonise.py): a sheet of every file with `sample` the aliquot, `case`, `sample_type`, `condition`, `batch`, preservation, strand class and share, library size, `chosen` and `drop_reason`; int32 counts of every file, for BRCA 93,695,475 bytes of Parquet and 204,748,235 of CSV; genes; the annotation log; GDC cases, treatments and molecular tests; the CDR rows. Consumers select only through `gdc.chosen()`.

Annotations: none of BRCA's 6,934 redacts a counts file, and the 52 "Item flagged DNU" are DNA aliquots. Logged per patient: neoadjuvant therapy 7, prior malignancy 16, unacceptable prior treatment 11, a DCIS case the CDR also redacts, and 5 barcode observations. Notifications exclude no one [CONVENTIONS.md:98](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/CONVENTIONS.md); 28 of the 1,088 carry one of the first three.

### 1.2 The clinical table

One row per patient [TODO.md:284](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md), built by `ctx.clinical.harmonise` from gather's tables:
- From the CDR: endpoints, age, stage, histology, menopause; months are days over 365.25/12. GDC agrees on stage for 1,072, both missing 24, and on age, in whole years, for 1,066.
- From GDC's IHC tests: ER positive 801, negative 237, equivocal 2 and conflicting 8 of 1,048 cases, the last two null with a reason; PR alike; HER2's IHC and FISH rule is set at build.
- `[Not Available]`, `[Discrepancy]`, `[Unknown]`, `[Not Evaluated]`, `#N/A`, `Stage X`, `not reported` and `Unknown` become null plus a `<column>_missing` reason: "stage missing 24: 13 Stage X, 6 discrepancy, 5 not available".
- In the cohort, hormone therapy is "yes" for 517 and unrecorded for 571, never "no", so an absent record stays missing.
- An event flag without a time leaves that analysis, logged: BRCA OS 1 censored row, KIRC PFI 2 events, LUAD OS 4.
- CDR redactions remove 5 BRCA patients, 4 with counts.

Flow: 1,095 tumour cases, 1,093 after the rule, 1,092 in the CDR, 1,088 after redactions.

| Endpoint | Patients | Events | Table 3, events and censored |
|---|---|---|---|
| OS | 1,087 | 151 | 151, 946 |
| PFI | 1,087 | 145 | 145, 952 |
| DFI | 945 | 84 | 84, 869 |
| DSS | 1,068 | 83 | 83, 995 |

Table 3 counts flags, the analysis rows with a time; the test asserts both. [Liu et al. 2018](https://doi.org/10.1016/j.cell.2018.02.052) was read through PubMed Central's [E-utilities](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=6066282).

### 1.3 Files

On this Mac, per the fetch logs [TODO.md:3343](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md): `tcga_brca` 4.9 GB, `tcga_kirc` and `tcga_luad` 2.4 GB, `tcga_lusc` 2.2 GB, each with cases, annotations and the CDR, plus miRNA and isoform folders. Optional downloads, each needing Rytis's permission, sized by GDC's PanCan manifest and GEO's FTP listing:

| File | Source | Bytes | Licence | For |
|---|---|---|---|---|
| clinical_PANCAN_patient_with_followup.tsv | GDC `/data/0fc78496-818b-4896-bd83-52db1f533c5c` | 18,633,685 | GDC open access | its ER field, UNVERIFIED, against Liu's p |
| TCGA_mastercalls.abs_tables_JSedit.fixed.txt | GDC `/data/4f277128-f793-4354-a13d-30cc7fe9f6b5` | 901,812 | GDC open access | purity |
| GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz | NCBI GEO FTP | 591,676,211 | UNVERIFIED [TODO.md:103](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md) | SCAN-B |
| GSE96058-GPL11154 and -GPL18573 series matrices | NCBI GEO FTP | 218,719 and 42,331 | UNVERIFIED | SCAN-B outcomes |
| PAM50 labels | no verified source | | | subtype [TODO.md:670](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md) |

Pages drawing GDC data acknowledge the TCGA project and the GDC [TODO.md:2594](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md); no wording is in the repo, so it is copied from GDC into `ui.toml` at build.

## 2. Analyses and figures

The preview's contracts [preview_clinical.py](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/scripts/python/preview_clinical.py) on real data, with CONVENTIONS' statistics [CONVENTIONS.md:56](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/CONVENTIONS.md):

| Figure | Statistics | Measured on BRCA |
|---|---|---|
| ER exercise, KM under all four endpoints | log-rank, Cox HR, Liu's p beside ours | below |
| Expression by stage, ER, histology, menopause, age | Kruskal-Wallis with Dunn, Mann-Whitney, Spearman trend, BH, effect sizes, missing counts | ESR1 by IHC ER: AUC 0.960 |
| KM per build-list gene with risk table and censor marks; Cox table | log-rank, RMST at 60 months; HR per SD, alone and adjusted for age and stage; PH test on ranked time as R's `cox.zph`; Harrell's C | below |
| Cut-point scan, 20th to 80th percentile | minimum p beside the permutation-corrected p | |
| Power map | Hsieh-Lavori, alpha stated | 145 PFI events: 59.3 % power at HR 1.2 per SD, 88.5 % at 1.3 |
| Cox screen; confounding map | BH over the panel; Cramér's V over stage, ER, histology, age, plate, 40 sites | |
| Treatment selection | descriptive | hormone therapy recorded for 473 ER-positive, 13 ER-negative |

Sanity checks, measured while planning, so not pre-registered:

| Check | Published | Ours |
|---|---|---|
| Event counts | Liu Table 3 | 16 of 16 equal |
| Median follow-up | Liu Table 2: 27.7 months | 27.7 |
| ER-negative worse; p for OS, PFI, DFI, DSS | Liu Figure 2: 0.097, 0.005, 0.001, 0.009 | 0.086, 0.0018, 0.0004, 0.0095; 797 ER-positive, 235 ER-negative, by GDC's field, not the paper's |
| Proliferation adverse, ER genes protective, disease extent adds | [Wirapati 2008](https://doi.org/10.1186/bcr2124), [Györffy 2010](https://doi.org/10.1007/s10549-009-0674-9), via PubMed | adjusted PFI HR per SD: MKI67 1.19, CI 1.01 to 1.41; ESR1 0.77, 0.65 to 0.91; PGR 0.79, 0.67 to 0.94; stage 2.39 per step, 1.90 to 3.01; 1,063 patients, 136 events |

## 3. The build

- `ctx clinical --dataset tcga_brca`, options `--endpoint PFI --split median --horizon 60 --covariates age,stage --strata --genes-from --min-log2cpm 1`, into `outputs/12_clinical/` [TODO.md:296](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md), reading gather's tables only. Build list in the dataset file: ESR1, PGR, ERBB2, MKI67 proposed.
- `ctx.clinical`: `harmonise`; `survival` and `score`, GS.7's one definitions from the preview with the PH test on R's ranks [gene_selection_plan.md:434](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/docs/plans/2026-09-30_gene_selection_plan.md); `cutpoint`; `page`; the JS twin `_surv.js.j2`. Text in `clinical.toml` and `gather.toml`, a card per option [CONVENTIONS.md:68](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/CONVENTIONS.md).
- Engine numpy, as GS.7 and D-049 decided; oracles statsmodels 0.15.0, installed, and a committed R 4.4.2 survival 3.7.0 fixture [gene_selection_plan.md:379](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/docs/plans/2026-09-30_gene_selection_plan.md). lifelines waits for notebook 10 [TODO.md:3651](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md); conda-forge's 0.30.0 accepts pandas 3, 0.30.1 and 0.30.3 do not. The lifelines boundary [CONVENTIONS.md:29](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/CONVENTIONS.md) and "fitted in `lifelines`" [clinical.toml:10](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/report/content/clinical.toml) change with a D-entry.
- openpyxl ==3.1.5, MIT, with Python 3.12 builds on conda-forge for linux-64, osx-arm64 and win-64, reads the CDR's xlsx and writes its fixture.
- Curve and surface contracts gain the controls and export hook that scatter, linked and bar have [save.py:146](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/viz/save.py).

Tests on GDC-shaped fixtures, barcodes prefixed for the cohort guard [decisions.md:151](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/knowledge/decisions.md), with the mutation each catches:

| Test | Catches |
|---|---|
| shuffled metadata columns, CRLF | positional reads, `\r` in values |
| a v38 gene-model line; two rows swapped; a nonzero PAR_Y row | another GENCODE; misaligned genes; versions stripped before the PAR_Y drop |
| library size; a stranded cohort; strand edges in exact arithmetic [lessons.md:62](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/knowledge/lessons.md) | `N_` rows summed; `unstranded` read for stranded libraries; float misclassification |
| FFPE and other-protocol aliquots that sort first | the tie-break run before either drop |
| a redaction, a notification, an md5 mismatch | notifications excluding; redactions ignored; unchecked files |
| filtered and unfiltered sheets before M10 | a consumer bypassing `gdc.chosen()` |
| CDR days, missing tokens, flag without time, redaction; conflicting ER | days/30; `Stage X` as a stage; `#N/A` as 0; the prototypes' 1,093 [gene_selection_plan.md:144](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/docs/plans/2026-09-30_gene_selection_plan.md); first test wins |
| GS plan tests 3 to 7 | Breslow for Efron; event-time PH ranks; scan bounds; another generator |
| one non-zero value, censored, shaped like ENSG00000248769 | β = -256.8 reported as converged [gene_selection_plan.md:170](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/docs/plans/2026-09-30_gene_selection_plan.md) |
| real data, opt-in, section 1's numbers typed in [lessons.md:71](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/knowledge/lessons.md) | rows dropped or duplicated |
| the page at real size: bytes, offline, contrast, lints, acknowledgement | a 15 MB overrun [lessons.md:23](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/knowledge/lessons.md) |

Budgets, Track B under 10 min and 4 GB [CONVENTIONS.md:28](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/CONVENTIONS.md), measured on BRCA:

| Step | Wall | Peak RSS |
|---|---|---|
| md5 of 1,231 files | 7.7 s | 54 MB |
| read to a matrix | 15.7 s | 479 MB |
| frame, Parquet, CSV | 1.2 s | 1.17 GB |
| log2 CPM, z, Cox fits, float64 | 15.1 s | 2.26 GB |

KIRC, LUAD and LUSC read in 8.3 to 9.3 s at 317 to 321 MB. The screen fitted 3,150 genes in 0.7 s [gene_selection_plan.md:181](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/docs/plans/2026-09-30_gene_selection_plan.md), so BRCA's 14,998 expressed genes take about 3.3 s. Page: the GS plan's 12,495,267 bytes [gene_selection_plan.md:307](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/docs/plans/2026-09-30_gene_selection_plan.md) plus the export's 1,607,890 [figure_controls_plan_approved.md:27](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/docs/plans/2026-10-01_figure_controls_plan_approved.md) leave 896,843 under 15 MB. Both are arithmetic.

## 4. Dependencies and order

- Now: ingest, clinical table, survival definitions; GS.1 is done [TODO.md:649](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md).
- Before links from M5: GS.2 and GS.4 [TODO.md:650](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md); meanwhile M10 takes pasted lists.
- In parallel once `ctx.gather.gdc` merges: M8's parsers, M7's CPTAC-2 dataset file, Step-3 DE through the same door.
- GS.8 after cohort bundles [TODO.md:635](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md); kidney and lung after breast; SCAN-B after its download.
- Breast DE sanity rules [TODO.md:721](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md) stay unseen: no tumour-normal comparison was computed.

## 5. Decisions for Rytis

1. **One tumour per case.** Recommended: drop FFPE, then other protocols, then take the first aliquot ID. First-barcode-only picks the same file wherever a frozen unstranded one exists; the rules differ on BRCA's 2 FFPE-only and LUAD's 2 stranded-only patients. Cost: those 4 leave; D-049's 1,089 becomes 1,087. Alternative: first barcode only keeps the 4, mixing protocols, by sort order alone.
2. **The breast panel until Step-3 DE.** Recommended: the most variable of BRCA's 13,457 expressed protein-coding genes up to the byte cap, about 3,600 by arithmetic. Cost: other genes show "no data on this page". Alternatives: waiting for Step-3 DE delays live selection; uint8 for every gene moved 4 to 12 patients across the median [gene_selection_plan.md:149](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/docs/plans/2026-09-30_gene_selection_plan.md).
3. **CDR data in pages.** The CDR article is CC BY-NC-ND 4.0 by its PubMed Central record, the course's materials CC BY 4.0 [CONVENTIONS.md:99](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/CONVENTIONS.md). Recommended: embed CDR columns verbatim, which ND permits for non-commercial use, compute months in the page, and print the licence and Liu's citation. Cost: those data stay non-commercial in the kit. Alternative: precomputed curves only, ending live KM and Cox. Whether GDC's terms supersede the article's licence is UNVERIFIED.

Assumptions, each overridable: gather is the only GDC reader, extending D-006; no exclusion beyond CDR redactions; 12 men kept and flagged [TODO.md:670](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md); openpyxl added; lifelines deferred; SCAN-B later, GEO's 3,273 [TODO.md:281](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md); PAM50 out until sourced [gene_selection_plan.md:162](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/docs/plans/2026-09-30_gene_selection_plan.md).

## 6. Steps

- [ ] M10.0 [S] Decisions; D-entries: GDC reader, tumour rule, `sex_at_birth`, CDR days, engine, licence notice.
- [ ] M10.1 [M] `conf/datasets/tcga_brca.toml`: `[gdc]`, `[guard]`, build list; schema in `ctx.dataset`.
- [ ] M10.2 [L] `ctx.gather.gdc` as section 1.1, with `chosen()`, its order test and the strand registry entry.
- [ ] M10.3 [M] `ctx gather --gdc-dir`: page, content, cards; opt-in four-project test; budget.
- [ ] M10.4 [M] `ctx.clinical.harmonise` as section 1.2; the Table 3 test.
- [ ] M10.5 [L] `survival`, `score`: Efron, `cox.zph` ranks, monotone refusal; `_surv.js.j2`; GS tests 3 to 7, R fixture.
- [ ] M10.6 [M] KM contract with risk rows and censor marks; controls and export on curves and surfaces.
- [ ] M10.7 [L] `ctx clinical` page: section 2, acknowledgement, lints.
- [ ] M10.8 [M] Live selection, GS.7: panel block, paste, tries counter, cross-study label, `--genes-from`.
- [ ] M10.9 [S] Real-size budget; both themes at 1280 and 375 px [lessons.md:24](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/knowledge/lessons.md); kit; review; lessons.
- [ ] M10.10 [M] Kidney and lung.
- [ ] M10.11 [M] Swimmer plot on the linked view; cumulative incidence.
- [ ] M10.12 [M] SCAN-B OS validation, after its download.
- [ ] GS.8 [M] Step-1 root, after cohort bundles.
