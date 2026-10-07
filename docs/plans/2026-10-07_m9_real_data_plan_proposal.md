# M9 Gene cards: external-data research and build plan

This plan is also at [PLAN.md](/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/plan_m9/PLAN.md).

**Outcome.** GTEx, HPA and Ensembl all answer a browser from any origin, so notebook 09 can join JupyterLite. Every panel fills from per-gene requests; the recommended plan downloads no bulk file. GTEx's licence, read today, allows our own figures with an acknowledgement, closing the UNVERIFIED clause of [CONVENTIONS.md:99](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/CONVENTIONS.md). Three decisions are open. The network layer can start beside M10.

Measured on NPTX2, ENSG00000106236, 2026-10-07, 13:04 to 13:10 UTC, unless marked arithmetic or UNVERIFIED. Responses: [plan_m9/](/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/plan_m9/).

## 1. Provider research

CORS: every response carried `access-control-allow-origin: *`, HPA's 302 included, and Chromium 152 fetched all three from another origin.

### GTEx Portal API v2, [openapi.json](https://gtexportal.org/api/v2/openapi.json)

| Item | Finding |
|---|---|
| Endpoints | `reference/gene`; `expression/geneExpression` with `attributeSubset=sex`; `expression/medianExonExpression`, `medianJunctionExpression`, `medianTranscriptExpression`; `dataset/collapsedGeneModelExon` |
| Bytes | per-sample TPM by sex 175,140; exons 61,223; junctions 48,219; transcripts 23,050; gene model 1,656; gene 468 |
| Values | TPM arrays: 19,616 samples, 54 tissues, 11 to 818 each; 13,164 male, 6,452 female, 8 sex-specific groups empty; no sample or donor ID, so the preview's donor hover goes |
| Release | `datasetId` in every row; `metadata/dataset`: gtex_v10 on GENCODE v39, 19,788 samples, 946 donors |
| Traps | expression endpoints default to gtex_v10, `tissueSiteDetail` to gtex_v8, `reference/gene` to GENCODE v26; pages hold 250 rows by default, at most 100,000; NPTX2's exon table has 270 |
| Rate limit | none in the spec or headers |
| Terms | [licence](https://gtexportal.org/home/license): acknowledge the Portal with a dated phrase; imply no endorsement; a redistributor keeps data current or states they do "not reflect the most current/accurate data available from the GTEx Portal" |
| Citation | [FAQ](https://gtexportal.org/home/faq#citePortal): NIH Common Fund support, and the Portal with the retrieval date |

GTEx V11, GENCODE 47 on V10's samples, exists as files only, [news](https://gtexportal.org/home/home/news?id=630); the API's `DatasetId` has no v11. NPTX2 is version 4 in v39, version 6 in our GENCODE 50.

### Human Protein Atlas, [data access](https://www.proteinatlas.org/about/help/dataaccess)

| Item | Finding |
|---|---|
| Endpoints | `/<ENSG>.xml`; `/<ENSG>.json`; `api/search_download.php` |
| Bytes | XML 1,064,622, 55,883 gzipped, 774,163 of it antibody sections; JSON 9,381; one search of 12 columns 4,032 |
| XML | 51 consensus and 40 HPA tissues in nTPM; 18 immune cell types plus total PBMC, per-sample values with sex; 154 single-cell types; IHC per tissue with reliability; ICC locations |
| JSON | categories; blood immunoassay 2490000 pg/L, an integer, mass spectrometry "220000", a string; taus "0.66" and "Not detected", strings |
| Search | cell-line nTPM by group column: `cell_RNA_leukemia` gave 93 lines, HL-60 0.0, THP-1 1.9, U-937 0.4; the XML holds no per-line values |
| Release | XML `entry version="25"`, Ensembl 109, GENCODE 43; [release history](https://www.proteinatlas.org/about/releases): 25.1 on 2026-05-25; `v25.proteinatlas.org` redirects to www today, `v24` answers 200; no documented rate limit, a search with too many results gets 400 |
| Terms | [licence](https://www.proteinatlas.org/about/licence): CC BY 4.0 for copyrightable parts; Ensembl, GTEx, TCGA, Tabula Sapiens, IntAct, BioPlex, OpenCell, TFClass, DepMap under user terms; BioGRID MIT; AlphaFold 3 CC BY-NC-SA 4.0 |
| Citation | a primary publication per section, proteinatlas.org, image credit, and `v25.proteinatlas.org/<ENSG>`; content never shown without citation |

Two licence facts reach the card: the consensus tissue RNA includes GTEx samples, [v25.0 notes](https://www.proteinatlas.org/about/releases), and single-cell data come from 35 sources while the XML names none per cell type.

### Ensembl REST

| Item | Finding |
|---|---|
| Endpoint | `lookup/id/<ENSG>?expand=1;mane=1;content-type=application/json`, content type in the query: 8,877 B, 6 transcripts, MANE Select ENST00000265634.4, NM_002523.3, 5 exons; 7.9 s; release from `info/data`, 116 |
| Rate limit | 55,000 per 3,600 s in headers; 429 with Retry-After, [wiki](https://github.com/Ensembl/ensembl-rest/wiki/Rate-Limits) |
| Terms | www.ensembl.org legal pages returned 403 to curl, so were not read. [EMBL-EBI terms](https://www.ebi.ac.uk/about/terms-of-use), revised 2024-02-05: no restriction beyond the data owners', attribution expected, overuse blocked |
| Citation | according to PubMed, Ensembl 2026, [doi:10.1093/nar/gkaf1239](https://doi.org/10.1093/nar/gkaf1239), which says release 116 is the last on the current platform |

The GENCODE 50 GTF registered at [tier1.toml:38](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/conf/datasets/tier1.toml) is Ensembl 116 by its header and tags 19,299 MANE Select transcripts, NPTX2's included. NPTX2's 5 MANE exons equal GTEx's 5 collapsed exons and its 4 introns equal GTEx's 4 junctions, base for base. REST adds nothing a card needs.

## 2. Data

| Item | Source | Bytes | Licence | Status |
|---|---|---|---|---|
| GTEx, 6 calls per gene | API v2, gtex_v10 | 309,756 raw | GTEx licence | fetched by `ctx genecard`, cached gzipped |
| HPA, 3 calls per gene | XML, JSON, search | 1,078,035 raw | CC BY 4.0 | same |
| MANE exons | GENCODE 50 GTF | 124,650,284 | unread | on disk, SHA-256 in tier1.toml |
| Own evidence | M5 results, gather genes table | on disk | ours | read |

57 genes are about 79 MB raw by arithmetic. No download needs permission. Bulk alternatives, bytes by HEAD: GTEx v10 median TPM 8,846,936; HPA `rna_tissue_consensus` 5,293,680, `rna_immune_cell_sample` 34,398,257, `rna_celline` 205,861,141, `proteinatlas.json.gz` 12,831,449.

Identity is the unversioned Ensembl ID, [decisions.md:98](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/knowledge/decisions.md). Tier1's genes table, 78,470 rows, names 35,370 genes by ID and has 446 symbols shared by 1,975 IDs. GENCODE 50 renamed NDUFA4L2 to COXFA4L2, ENSG00000185633, which [TODO.md:721-722](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md) still name.

## 3. Panels on real data

The preview's contracts carry over, [preview_genecard.py:106-213](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/scripts/python/preview_genecard.py). Reference panels are descriptive: n, median, IQR, tau.

| Panel | Data | On the card |
|---|---|---|
| Own evidence | M5, 10 contrasts; box of the featured contrast | BoxData; Wald, BH, D-012 call; Mann-Whitney, BH, rank-biserial r |
| Normal tissues | GTEx per-sample TPM by sex | StripData, log axis; donor control all, female, male |
| Exons, junctions, isoforms | GTEx medians and model; GTF MANE | HeatmapData with CV and junction medians; isoform shares; Primer-BLAST link |
| Immune cells | HPA XML | BarData by lineage; our tau beside HPA's |
| RNA against protein | HPA consensus nTPM, IHC level and reliability | drawn if decision 1 is yes |
| Single cell, subcellular, blood, cell lines | HPA | text with citation; lines named in the dataset file |
| Prognosis, TCGA RNA, structure, interactions | HPA, AlphaFold, STRING | linked; the JSON holds only an interaction count, 9 |
| Clinical row | M10 | after GS.7 |

Step-1's default card list, [gene_selection_plan.md:130](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/docs/plans/2026-09-30_gene_selection_plan.md), is the 7 registered box genes and the top 25 up and 25 down calls of ndn_gcsf: 57 genes, no overlap, 12 named only by ID.

Sanity checks, registered before the first real card:

| Check | Expected |
|---|---|
| NPTX2, HPA 25.1 | immunoassay 2.49 µg/L, the 2.5 of [TODO.md:716](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md); top of 154 single-cell types Corticotrophs, 335.9 nCPM; pituitary 200.0 nTPM |
| NPTX2, GTEx v10 | highest median Pituitary 218.6 TPM; lowest EBV lymphocytes 0.052 |
| NPTX2, tau | ours 0.6710 on the 51 consensus values; HPA prints 0.66; both shown, HPA's definition UNVERIFIED |
| NPTX2, structure | GTEx exons and junctions equal MANE exons and introns |
| ACTB | HPA "Detected in all" and "Low tissue specificity"; our tau below 0.671; UNVERIFIED until fetched |

NPTX2 and CA9 are in none of the 10 Step-1 results tables; their own-evidence panel says so. GAPDH is called up in 4 of 10 contrasts, ndn_gcsf log2 FC 1.32 at padj 1.6e-7: one pitfall line.

## 4. Build

Package `src/ctx/genecard/`: `gtex.py`, `hpa.py`, `annotation.py`, `cache.py`, `cards.py`, `page.py`, `run.py`, each under the 600-line limit of [.labtemplate.yml:33](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.labtemplate.yml). Parsers use only the standard library and numpy, like `ctx.function.lite`, so notebook 09 runs them in Pyodide. HTTP is `urllib` with an injectable opener, as at [setup.py:201](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/setup.py); no new dependency.

CLI: `ctx genecard [GENE...] --dataset --genes-from --contrast --gtex-dataset gtex_v10 --hpa-release v25 --offline --refresh`. The HPA default drops ".1" because a URL pins only the major release. `--genes-from` takes a GS.1 list or an M5 table by `gene_name`.

Cache: `resources/cache/<provider>/<release>/<ENSG>/<endpoint>.json.gz`, refining [TODO.md:274](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md) since one gene needs several responses. `resources/cache/manifest.json` holds URL, UTC time, status, bytes and the SHA-256 of each body as received. A refetch whose body differs under the same key stops the run; `--refresh` moves the old folder to `resources/cache/superseded/<date>/`. A failed call leaves an "unavailable" panel naming the error.

[genecard.toml](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/report/content/genecard.toml) is rewritten for the real panels and options.

Fixtures are the committed NPTX2 responses; network tests stay deselected, [pyproject.toml:43](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/pyproject.toml).

| Test | Mutation caught |
|---|---|
| Two-page exon response | default 250 rows, losing 20 of NPTX2's 270 |
| Stub refuses GTEx URLs without `datasetId` or `gencodeVersion` | drift to v8 or v26 |
| Empty `data` | read as zero expression, not absent |
| Join by unversioned ID | v39's .4 never meeting our .6 |
| HPA types and stamp | integer against string concentrations; "Not detected" tau; XML version 26 under key v25 |
| Cache drift, offline, writers | silent overwrite; a call under `--offline`; `resources/manifest.json` touched |
| Degraded panel | whole card dropped when HPA fails |
| Tau order | tau on linear nTPM, 0.9575 for NPTX2; total PBMC counted |
| Exon and junction match | off-by-one intron coordinates |
| Symbols | first ID of an ambiguous symbol taken; NDUFA4L2 accepted |
| Card cap, budget-marked, 57 cards | a count-based cap; a page over 15,000,000 B |
| Footer | missing release, date or versioned URL; unfilled placeholder |

Budgets, Track B: under 10 min and 4 GB. Calls took 0.21 to 0.49 s at GTEx and 0.06 to 0.13 s at HPA, so 57 genes take about 190 s at the slowest, by arithmetic. Render time and peak RAM are UNVERIFIED.

| GTEx samples in the page as | Bytes per card, NPTX2 | Cards that fit |
|---|---|---|
| int16 byte planes, gzip, base64 | 65,250 | 87 |
| JSON, 3 significant digits | 120,463 | 47 |
| summaries per tissue and sex | 22,351 | 254 |

Fit is arithmetic: 15,000,000 less the kit v6 preview's 5,231,627, the evidence parts of [gene_selection_plan.md:287-296](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/docs/plans/2026-09-30_gene_selection_plan.md), and export assets of 507,629 plus 195,608, [figure_controls_plan.md:27](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/docs/plans/2026-10-01_figure_controls_plan_approved.md) and [TODO.md:439](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md), leaves 5,687,673. IHC and text bytes are not counted.

## 5. Dependencies and order

- Work order item 12 follows M10, [TODO.md:449-450](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md); GS.9 follows the spike and M10 and uses GS.2 to GS.4, [TODO.md:685](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md).
- In parallel now: steps 1 to 7, new files only, no page script, no pixi change.
- After M10 merges: steps 8 to 11, sharing the page script and the controls kinds of [figure_controls_plan.md:171-181](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/docs/plans/2026-10-01_figure_controls_plan_approved.md).

## 6. Decisions for Rytis

1. **Draw GTEx data, from its API or inside HPA's consensus, as our own figures?** Recommended: yes, with the dated acknowledgement and a line that newer GTEx releases exist; amend CONVENTIONS.md:99 in a D-entry. Cost: one footer line. If no: tissue, exon and RNA-against-protein panels become links, and the qPCR brief loses its data.
2. **GTEx release.** Recommended: gtex_v10 through the API. Cost: GENCODE v39 against our 50; genes newer than v39 get no GTEx panel, and v39's coverage of the 12 ID-named top calls is UNVERIFIED. Alternative: V11 files on GENCODE 47, 2.4 GB of TPM and a 20.3 GB exon parquet, with no API.
3. **GTEx samples in the page.** Recommended: every sample, int16 gzip, 65,250 B per card, so the approved donor control recomputes quantiles in the page. Cost: 87 cards fit; DecompressionStream, with M10's fallback. Alternative: summaries, 22,351 B and 254 cards; the control becomes three precomputed options and the points go.

Assumed unless you object: MANE from the GTF, not REST; single-cell as text; interactions linked; notebook 09 ships precomputed cards and fetches one gene live when online.

## 7. Steps

1. [ ] M9.0 On approval: plan into `docs/plans/`; D-entry; CONVENTIONS.md:99; TODO M9 and lines 721-722. S
2. [ ] M9.1 `cache.py`: layout, manifest, drift stop, `--refresh`, `--offline`. M
3. [ ] M9.2 `gtex.py`: client, parser, paging, version check; NPTX2 fixtures with provenance README. M
4. [ ] M9.3 `hpa.py`: XML, JSON, search; stamp, coercion, one row per ID. M
5. [ ] M9.4 `annotation.py`: MANE exons and introns from the GTF; GTEx match; the preview's qPCR pair rule. S
6. [ ] M9.5 Resolution and card list: GS.1 parser; ambiguous and unknown symbols refused; order by padj, absolute log2 FC, then ID. S
7. [ ] M9.6 Per-gene tables with README; `ctx genecard` with cards and degraded panels. M
8. [ ] M9.7 After M10: card page, footer, offline test, NPTX2 and ACTB checks. L
9. [ ] M9.8 GS.9: evidence view, box, cards packed to the cap, missing-card command, entry page, controls F1 to F4; budget test at 57 cards. L
10. [ ] M9.9 Notebook 09 and its lite bundle. M
11. [ ] M9.10 GS.10's M9 part: four browsers, 1280 and 375 px, both themes; kit; sizes, wall time, RAM; review, lessons. M

## Not done

No Pyodide run, no `file://` origin test, no second gene: ACTB, GTEx resolution by unversioned ID and v39 coverage of new genes are UNVERIFIED. HPA's tau definition, BioGRID access and IHC bytes were not researched. Requests: 8 GTEx, 4 HPA, 2 Ensembl, 4 browser fetches, 13 HEAD; no data file downloaded. A web search found no GTEx rate limit, [gtexr docs](https://docs.ropensci.org/gtexr/).
