# Real data for M4, M10, M8, M9 and M7, and the custom-data guide: decisions for Rytis

Written 2026-10-07 from the six plans in this folder, each a proposal until answered. Every
number is the planner's measurement, cited in its plan. The recommended answer comes first.

## What each plan found

| Plan | Real data | Ready now | Waits for |
|---|---|---|---|
| [M4](2026-10-07_m4_real_data_plan_proposal.md) | Step-1's 210 runs, Salmon outputs nf-core already wrote | library type, decoys, mapping, sample correlation from files on this Mac | Q0, a read-only listing on Kebnekaise; then Q1, a copy of the per-sample `quant.sf` files; the cross-method figure needs STAR on the same reads (M3.11) |
| [M10](2026-10-07_m10_real_data_plan_proposal.md) | TCGA-BRCA, 1,088 patients; GDC clinical and the TCGA-CDR | everything: no download; the CDR's event counts equal Liu et al. 2018 Table 3 in 16 of 16 | nothing; it builds the one GDC reader M8 and M7 use |
| [M8](2026-10-07_m8_real_data_plan_proposal.md) | Step-1 lncRNA by biotype; TCGA miRNA, BRCA first | the Step-1 half | M10's GDC reader, M5 on TCGA, about 66 MB of target databases |
| [M9](2026-10-07_m9_real_data_plan_proposal.md) | GTEx, HPA, Ensembl per gene, cached | everything: no bulk download; GTEx's licence read today allows our own figures | M10 first, by the work order |
| [M7](2026-10-07_m7_real_data_plan_proposal.md) | CPTAC prospective BRCA, 125 tumours; Rieckmann 2017 for Step-1 | the Step-1 part, after a new M5 contrast | M10's GDC reader; downloads |
| [Guide](2026-10-07_custom_data_guide_plan_proposal.md) | a worked example of 4 Step-1 runs | the tooling: site files, a `Code` block, the trace reader | the laptop proof (2.2 GB of downloads) and a Kebnekaise rehearsal |

## Decisions

### M4 Quantification

| # | Recommended | Cost | Alternative |
|---|---|---|---|
| 1 | Featured dataset Step-1 | single-end effective lengths rest on Salmon's prior | a breast dataset, which has no dataset file yet |
| 2 | Q1 copies every Salmon file per sample, after Q0 measures the size | disk on the Mac, size unknown until Q0 | merged tables only, which lose per-sample ambiguity |
| 3 | Gene summaries by ctx's own tximport port, checked against R tximport | about 40 lines and a D-entry, since CONVENTIONS names pytximport | pytximport: GPL-3.0, 5 packages and pandas, which CONVENTIONS bars there |
| 4 | Cross-method comparison: M3.11 adds featureCounts; Q2 runs Salmon on the same 8 trimmed reads | one cluster job, 4 cores, 16 GB, 1 h | the cohort's own Salmon, which mixes in trimmer differences |
| 5 | Track A, plus a Track B test on the cohort path | none | |
| 6 | Salmon defaults are nf-core's, bias correction off and offered as a knob | the preview's cards, which say it is on, are corrected | |
| 7 | F4, reads against TPM by length, carries no test, as M3's R6 exception | none | |

### M10 Clinical and survival, with the GDC reader

| # | Recommended | Cost | Alternative |
|---|---|---|---|
| 1 | One tumour per case: drop FFPE, then other library protocols, then take the first aliquot | 4 patients leave (BRCA 2, LUAD 2); D-049's 1,089 PFI patients become 1,087 | first barcode only, which keeps them and mixes protocols by sort order |
| 2 | Until Step-3 DE: the most variable of BRCA's 13,457 expressed protein-coding genes up to the page cap, about 3,600 | other genes show "no data on this page" | wait for Step-3 DE, delaying live selection |
| 3 | Embed the CDR's columns verbatim, compute months in the page, print the CC BY-NC-ND 4.0 licence and Liu's citation | those data stay non-commercial in the kit; whether GDC's terms supersede the article's is UNVERIFIED | precomputed curves only, ending live KM and Cox |

Assumed unless you object: gather is the only GDC reader; no exclusion beyond the CDR's redactions; the 12 men kept and flagged; survival fitted in numpy with statsmodels and an R fixture as oracles, lifelines deferred to notebook 10; SCAN-B and PAM50 later.

### M8 Non-coding RNA

| # | Recommended | Cost | Alternative |
|---|---|---|---|
| 1 | BH within biotype over M5's tested genes, no refit | none | a refit per biotype, 128.8 s and 3.57 GB per BRCA contrast |
| 2 | Validated targets are miRTarBase's strong evidence | 3,348,778 B | all human interactions, 337,103,345 B, with CLIP-only pairs to filter out |
| 3 | The kit ships derived target rows; the raw files are linked | students cannot re-derive targets offline | shipping 21,807,648 B under terms TargetScan does not publish |
| 4 | miRNA against target within tumours, normals apart | small normal strata, 20 in LUAD | pooled, where the tumour-normal gap drives rho (-0.450 against -0.236 for ISCU) |
| 5 | The biotype volcano's control after GS.3 | the figure ships first without it | a second column table now, merged later |
| 6 | The library-type lesson from 9 LUSC tumours with both libraries | one GDC metadata query in the reader | no sample shows both library types |
| 7 | Download GDC's GENCODE v36 GTF, the TCGA counts' own annotation | 44.5 MB | GENCODE 50 by stable ID, losing 1,128 of 16,901 v36 lncRNAs |

### M9 Gene cards

| # | Recommended | Cost | Alternative |
|---|---|---|---|
| 1 | Draw GTEx data as our own figures, with the dated acknowledgement and a line that newer releases exist; CONVENTIONS' UNVERIFIED clause amended | one footer line | links only: tissue, exon and RNA-against-protein panels become links |
| 2 | GTEx v10 through its API | GENCODE v39 against our 50; genes newer than v39 get no GTEx panel | V11 as files only: 2.4 GB of TPM and a 20.3 GB exon table, no API |
| 3 | Every GTEx sample in the page, packed, so the donor control recomputes in the page | 87 cards fit under 15 MB | summaries: 254 cards, the control becomes three fixed options |

Assumed: MANE exons from the GTF, single-cell as text, interactions linked, notebook 09 ships precomputed cards and fetches one gene live when online.

### M7 RNA versus protein

| # | Recommended | Cost | Alternative |
|---|---|---|---|
| 1 | RNA from GDC's STAR counts through M10's reader; protein and annotation from the paper's data freeze | 543 MB and waiting for M10 | the freeze alone: 35 MB, but hg19 FPKM that M5 cannot test |
| 2 | Approve the downloads below | see the list | without the supplements, the published-median sanity rule cannot be registered |
| 3 | The kit ships derived tables and figure data only | none | shipping the freeze rests on an agreement silent on redistribution |
| 4 | Protein calls use D-012's rule, padj < 0.05 and abs(log2 FC) of at least 1 | ratio compression hides protein calls, which a control shows | 0.5 for protein, a second definition of "changed" |
| 5 | A new Step-1 contrast, neutrophils against monocytes in steady-state blood, paired, 9 against 9 | one more M5 contrast; the 92 recorded metrics checked again | a second DE engine inside M7 |
| 6 | One Excel reader for M7 and M10: openpyxl 3.1.5, MIT, which also writes the test fixtures | one environment change | M7's planner proposed fastexcel and M10's openpyxl: two readers, one definition broken |

### The custom-data guide

| # | Recommended | Cost | Alternative |
|---|---|---|---|
| 1 | One generated page, "Your own data", from `ctx guide`, every command tested; READMEs point to it | medium | a docs/ file, a second copy checked at best by pattern |
| 2 | Counting on a laptop above a floor measured on this Mac; the cluster below it | you raise this Mac's Docker memory, 11.7 of 24 GiB now, under the 16.8 GiB index peak seen once | laptops only fetch, or a transcript-only index whose counts are not comparable with Step-1's |
| 3 | One launcher per pipeline, with a site file for the laptop, Kebnekaise or another SLURM site | G.1's tests rerun | separate laptop scripts, second copies of the checks |
| 4 | Worked example: the 4 smallest Step-1 runs, 2 neutrophils against 2 monocytes, confounding stated | 1.27 GB of reads | 2 paired donors, 3.79 GB |
| 5 | A newcomer reproduces it from the page before the release | one colleague's afternoon | after the release |
| 6 | Head-job memory 12G in rnaseq.sbatch (P4.5) | none on a site with the memory | 8G, peaks at 6.16 to 7.49 GiB already |

## Downloads that need your permission

None is fetched until you approve it.

| For | File | Source | Bytes | Terms |
|---|---|---|---|---|
| M8 | `Predicted_Targets_Context_Scores.default_predictions.txt.zip` | TargetScan 8.0 | 17,673,917 | none published; cite McGeary 2019 |
| M8 | `miR_Family_Info.txt.zip` | TargetScan 8.0 | 162,636 | same |
| M8 | `Gene_info.txt.zip` | TargetScan 8.0 | 622,317 | same |
| M8 | `miRTarBase_SE_WR.csv` | miRTarBase 10.0 | 3,348,778 | public-domain disclaimer; cite Cui 2025 |
| M8 | `gencode.v36.annotation.gtf.gz` | GDC | about 44.5 MB | GDC open access |
| M7 | `S060_Breast_Cancer_data_freeze_GCTfiles_v5.4-public.tar.gz` | PDC, PDC000120 | 34,878,978 | CPTAC data use agreement: acknowledgement |
| M7 | 121 CPTAC-2 breast STAR counts files | GDC | 508,400,341 | GDC open access |
| M7 | `c5.go.cc.v2026.1.Hs.symbols.gmt`, `c2.cp.reactome.v2026.1.Hs.symbols.gmt` | MSigDB 2026.1 | 778,106 and 898,680 | CC BY 4.0 |
| M7 | Rieckmann 2017 Supplementary Table 6 and its supplementary PDF | nature.com | about 10.9 MB and 2.6 MB | reuse terms UNVERIFIED; not shipped |
| M7 | Krug 2020 supplement PDF | PubMed Central | UNVERIFIED | all rights reserved; read, not shipped |
| M4 | Q1, the Salmon files from Kebnekaise | your cluster | Q0 measures it | the study's IDF names no licence |
| M4 | nf-core yeast test data at 626c8fab, to rebuild the test fixture | GitHub | 24,675,021 | MIT |
| Guide | 4 Step-1 runs, GENCODE 50 genome and GTF, nf-core images, for the laptop proof | ENA, EBI, quay.io | 1,268,992,130 + 845,635,028 + 124,650,284, images unmeasured | open |
| M10, later | `clinical_PANCAN_patient_with_followup.tsv`; ABSOLUTE purity calls | GDC | 18,633,685; 901,812 | GDC open access |
| M10, later | SCAN-B, GSE96058 expression and series matrices | NCBI GEO | 591,676,211 and 261,050 | UNVERIFIED |

## Cluster commands

Only M4's Q0 is ready now; it reads, writes nothing, and its output decides Q1's size. It is in
the M4 plan, section 1. The guide's Kebnekaise rehearsal comes after its tooling is built.

## Order

Your priority order holds: M4, M10, M8, M9, M7, the guide alongside. Planned in parallel
where files allow: M10's GDC reader and M4's tximport port first; M8's Step-1 half, M9's
network layer and the guide's site files beside them; M8's TCGA half and M7's tumour part
after the reader merges.
