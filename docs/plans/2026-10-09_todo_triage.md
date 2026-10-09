# TODO.md triage for work-order item 13 (main ecf03fc, 2026-10-09)

Read-only. Every unchecked `- [ ]` and `N. [ ]` in TODO.md: 123 rows (grep `^\s*(- |[0-9]+\. )\[ \]`).
Evidence is a commit, a file:line or a test. Line numbers are TODO.md's at ecf03fc.

## Counts

| Status | Rows |
|---|---|
| DONE | 28 |
| SUPERSEDED | 10 |
| BUILD | 52 |
| RYTIS | 20 |
| BLOCKED | 13 |
| Total | 123 |

## Finding that changes several rows: CI runs again, and main is red

GitHub Actions now runs jobs (run 37807649888, 2026-10-08, 21 min), so the "CI billing" block of K.5 and D-047 no longer holds. No run in `gh run list` succeeded. Run 37807649888 on main: 5 failed, 1,312 passed, 9 skipped, 72 errors, 1,171 s. Causes read from its log (copy in this folder, `ci_37807649888.log`):

| Tests | Cause in the log |
|---|---|
| 72 errors in `test_export_browser.py`, `test_heatmap_browser.py`; `test_qc_fastp_grid.py::test_the_page_draws_the_floor_grey_and_every_band` | Chromium: "No usable sandbox!" on ubuntu-latest, whose AppArmor blocks unprivileged user namespaces (`tests/browser.mjs`) |
| `test_io.py::test_provenance_not_from_cwd` | `git commit` exit 128 in a temp repo: no git identity on the runner |
| `test_runtime.py::test_threads_applied_before_import`, `test_cli.py::test_threads_set_before_import` | BLAS pool 2 where 3 expected on the 4-core runner; cause not established |
| `test_runtime.py::test_child_peak_attributed` | parent peak 1,110,167,552 B against the 256 MiB bound on Linux; cause not established |

The containers workflow has never run; `containers/analysis.Dockerfile` has `build` and `production` stages and no `figures` target.

## Table

Pkg is the BUILD package below. Size S/M/L/XL.

### Priority (Rytis, 2026-10-07)

| Line | Id | Item | Status | Evidence, or what remains |
|---|---|---|---|---|
| 439 | P1 | M4 on real data | BUILD | Open parts buildable: thinned heatmap's tall triangle (S), strip colours repeating with 7 cell types (S), M4.4 `salmon.py` and probe + M4.11 laptop exercise on the 8 teaching runs already on the Mac (M), M4.15 preview through `page.py` (S), browser 1280/375 (S); Pkg C1. "A Track B budget test" is DONE: `tests/test_budget.py:225`. F5 and F8 BLOCKED on M3.11 and Q2 (RYTIS Q1). M4.14 controls in Pkg E |
| 440 | P2 | M10 on real data | BUILD | Pkg C2: Dunn's test (reuse `ctx.proteo.measure`, src/ctx/proteo/measure.py:115, one definition) S; retire the made-up preview's own Cox and KM S; box labels at 375 px S; kit rebuild with TCGA roots S; M10.10 kidney and lung (all GDC folders on the Mac: ~/ctx_kidney, ~/ctx_lung) M; M10.11 swimmer and cumulative incidence M; GDC acknowledgement wording, row 705. Cox screen on the live selection and live captions in Pkg E. M10.12 SCAN-B is RYTIS (download) |
| 441 | P3 | M8 on real data | BUILD | Pkg C2: F1 by biotype S, LUSC histone rule measured at the kit rebuild S, BRCA peak RSS 2.30 against 3.77 GB (row 746) M, named-gene box order differing between runs (a determinism defect) S. F2 preset and F7 controls in Pkg E; F3's control BLOCKED on GS.3. Plan step 13 (miRNA, mRNA, protein on CPTAC miRNA) needs a download, RYTIS Q6 |
| 442 | P4 | M9 on real data | BUILD | Pkg C1: HPA cell-line search S; NDUFA4L2 in TODO is COXFA4L2 in GENCODE 50 (resolve by previous symbol, state it) S; browser check S. M9.8 (GS.9) L in Pkg E; M9.9 BLOCKED on Phase 5 lite-build; M9.10 BLOCKED on M9.8 and Safari. Junction panel: no V11 junction file on disk, RYTIS Q6. HPA's own tau absent from v25: RYTIS Q7 |
| 443 | P5 | M7 on real data | BUILD | Pkg C1: F4 (M5 on a GDC sheet now exists, D-058) M; F5 heatmap M; the two freeze columns named by aliquot S; contours from the existing `DensityData` (src/ctx/viz/density.py:59) S; P7.11 budget test and kit S. P7.10 controls in Pkg E |
| 444 | P6 | Custom-data guide | RYTIS | U.1 to U.9 DONE (cd0d38d, a75a4fb); U.10's fixes merged (f735001, 274947b) but Part 2 cannot run: Docker has 12,525,977,600 B (11.67 GiB) against the 16.8 GiB index peak (src/ctx/report/guide.py:106), Q3. U.11 Kebnekaise rehearsal Q2, U.14 newcomer Q10. U.12 and U.13 are BUILD S in Pkg G |

### Work order (2026-10-01)

| Line | Id | Item | Status | Evidence, or what remains |
|---|---|---|---|---|
| 456 | WO3 | P2Q.9 QC bands | BUILD | L, Pkg B; same work as row 639 |
| 457 | WO4 | P2S.4, P2S.5 | BUILD | L, Pkg D; rows 657, 660 |
| 458 | WO5 | M1, M2 controls | BUILD | M, Pkg E; `ctx.qc.identity` first (Pkg B); item 2's leftovers (M1 filter selects labelled CELL_TYPE and CONDITION; depth legend cut at 375 px) S, Pkg E |
| 459 | WO6 | M5 controls | BUILD | M, Pkg E; per-contrast volcano pages after GS.2 to GS.4 |
| 460 | WO7 | M6 controls | BUILD | M, Pkg E; after GS.5 |
| 461 | WO8 | M3 controls | BUILD | M, Pkg E: sample picker, gene-type filter on points, live strand thresholds (D-051 consequence) |
| 462 | WO9 | GS.2 to GS.4 | BUILD | L; GS.3 in Pkg D, GS.2 and GS.4 in Pkg E |
| 463 | WO10 | M8 controls | BLOCKED | Neighbour drop-down DONE (2b2246f; `templates/_pairs.js.j2`); the biotype volcano's cut-offs wait for GS.3 (Pkg D), then Pkg E |
| 464 | WO11 | M10, breast first | SUPERSEDED | by Priority 2: built and merged bc7426b, D-055; remainders listed at line 440 |
| 465 | WO12 | External research, then M9 | SUPERSEDED | by Priority 4: research in docs/plans/2026-10-07_m9_real_data_plan_proposal.md section 1 (CORS at line 11), M9 merged 7e6e1af, D-056; GS.9 is M9.8 |
| 466 | WO13 | This triage, then build | BLOCKED | closes when every other row closes |

### Phase 2

| Line | Id | Item | Status | Evidence, or what remains |
|---|---|---|---|---|
| 516 | M1c | Step-1 metadata, curation, gather, counts | DONE | `conf/datasets/tier1.toml`; D-046 evidence: gather 210 samples, 78,470 genes, job 47485533, rnaseq 3.27.0 a1fcddd. "Released with the teaching subset" is the Zenodo release, RYTIS Q10 |
| 517 | M5 | pydeseq2 on Step-1 | DONE | D-046 (R1 to R8, PDAC naive and within-batch recorded); ~/ctx_tier1/outputs/05_de holds ndn_pdac_naive, _adjusted, _within, pca by read length. The "precomputed knob grid" is D-005's JupyterLite grid, unbuilt, carried by row 711 |
| 518 | M6 | ORA/GSEA, lite against gseapy | DONE | tests/test_enrich.py:40 (ORA), :128 and :153 (prerank, NES against gseapy); D-046 R4, R5 |
| 519 | Reports | M1, M5, M6, hub; sanity | DONE | D-046: pages 4,963,184, 12,773,948, 5,365,560, 67,666 B; record ac4a0b5; 7 of 8 rules pass, R4 fails as registered and is taught |

### Phase 2M

| Line | Id | Item | Status | Evidence, or what remains |
|---|---|---|---|---|
| 562 | M3.10 | Fake-STAR refusal test in CI on linux-64 and the Mac | BUILD | S, Pkg A. CI runs again but is red; the test (tests/test_process_align.py:131) is not among the failures of run 37807649888, and `-q` output does not show that it ran rather than skipped: add `-rA` or a junit artifact, make main green |
| 563 | M3.11 | Cluster full-genome STAR pilot, 8 runs | RYTIS | Q1. An agent writes the cluster block first (S) |
| 564 | M3.12 | M2b judges M3's alignment table | BUILD | M, Pkg B. `conf/qc_metrics.toml:377,397` have `stage = "align"` metrics; `src/ctx/qc/run.py` reads no align table |
| 565 | M3.13 | Full-cohort M3 on Kebnekaise | RYTIS | Q1, after M3.11 measures RAM, wall time and disk |
| 566 | M3.14 | Suite back under 5 min | BUILD | M, Pkg A. CI 1,171 s (run 37807649888); 1,445 s locally on 2026-10-07 |

### Phase 2P, 2V, 2G, 2K, 2L

| Line | Id | Item | Status | Evidence, or what remains |
|---|---|---|---|---|
| 581 | P2P.3 | M5, M6, hub to content; `input` required; fields renamed | BUILD | M, Pkg D. `ModulePage` still has question, concept, practice, think, mistakes and `input: str = ""` (src/ctx/report/builder.py:299-313); de/page.py:17 still holds prose |
| 588 | P2V.2 | Made-up previews, build script, hub 11 of 11 | DONE | f924be9, 891ab62, scripts/python/build_preview.py and preview_*.py; "published" is superseded by D-050 (made-up pages leave at the release) |
| 595 | G.2 | `Code` result kind | DONE | src/ctx/report/builder.py:156; cd0d38d (U.5) |
| 596 | G.3 | dataset_template.toml | DONE | conf/dataset_template.toml; tests/test_dataset_template.py:56 |
| 597 | G.4 | guide and trace modules | SUPERSEDED | by U.3 and U.9: src/ctx/gather/trace.py, src/ctx/report/guide.py, content/guide.toml (not cluster.toml), tests/test_trace.py, tests/test_guide.py; cd0d38d |
| 598 | G.5 | Kebnekaise rehearsal | RYTIS | Q2 (= U.11) |
| 599 | G.6 | CONVENTIONS student route; items closed; D-entries | BUILD | S, Pkg G (= U.12). CONVENTIONS.md:31 names only the lab launcher; the "D-046 and D-047" G.6 names were taken by other decisions, so the guide's D-entry takes the next free number |
| 600 | G.7 | Guide browser check, every Copy | BUILD | S, Pkg G (= U.13); Track B measured in U.6 |
| 601 | G.8 | Newcomer reproduces from the page | RYTIS | Q10 (= U.14) |
| 610 | K.4 | Kit proof on Windows and macOS | RYTIS | Q10: a person at a Windows screen, and a browser-downloaded zip, which waits for the release |
| 611 | K.5 | Containers in the kit | RYTIS | Q9. Billing no longer blocks; the `figures` target does not exist (Pkg F builds it) |
| 614 | K.6 | Release kit with no made-up page | BUILD | S, Pkg C2: a release mode in `scripts/python/build_kit.py` that refuses `--made-up` and asserts 0 Made-up links, with a test; the release itself is Q10 |
| 625 | L.6 | Short all-capital aliases matched by case | BUILD | S, Pkg B. `_names()` matches every alias without `_EXACT` case-insensitively (src/ctx/report/glossary_lint.py:88-99) |

### Phase 2Q and P4.5

| Line | Id | Item | Status | Evidence, or what remains |
|---|---|---|---|---|
| 633 | P2Q.3 | Read fates and composition | BUILD | S, Pkg B. Built: bubble (src/ctx/qc/figures.py:43), Sankey (:87), status table (`*_qc_tier1_status_*`). Not built: the MT % against rRNA % composition scatter; after P2Q.10 corrects rrna_fraction |
| 635 | P2Q.5 | `ctx qc` page | DONE | M2.6a to M2.6d, M2.7, M2.8 ticked: Step-1 210 samples, 50.2 s, 1.71 GB, page 6.76 MB; planted page 5,156,086 B checked at 1280 and 375 px |
| 637 | P2Q.7 | VST heatmap, PCA brushing in M5; set_scores and Kruskal-Wallis in M6 | BUILD | L, Pkg D. No heatmap in src/ctx/de, no `set_scores` in src/ctx/function |
| 638 | P2Q.8 | M2 read stage on 210 runs | RYTIS | Q1, with M3.13 |
| 639 | P2Q.9 | QC background bands | BUILD | L, Pkg B. `plot_bands` appears in registry.py (4) and qc_metrics.toml (2) only, 0 in src/ctx/viz/strip.py; the 210 raw FastQC zips are on the Mac (210 files in ~/ctx_tier1/data/raw/tier1_rnaseq/fastqc/raw) |
| 640 | P2Q.10 | Correct top_gene_share and rrna_fraction | BUILD | M, Pkg B. Definitions still as before (conf/qc_metrics.toml:541-547, 598-604). Step-1 is libType U, so the approved route counts 18S through the antisense lncRNAs (docs/plans/2026-10-01_qc_bands_proposal_approved.md:69-71). The bands then go to Rytis, Q5 |
| 644 | P4.5 | Head-job memory 12G | DONE | scripts/bash/rnaseq.sbatch:5 `--mem=12G`, by the recommended answer (TODO.md:444, cd0d38d). The recorded blocks of past runs (TODO.md:2303, 2572, 3284) keep the 8G they ran with |

### Phase 2D

| Line | Id | Item | Status | Evidence, or what remains |
|---|---|---|---|---|
| 650 | 2D.run | `ctx run` | BUILD | S, Pkg F. Built: src/ctx/cli.py:520, tests/test_pipeline.py:105 (stops at first failure), :186 (tables byte-equal). The proof runs gather, qc and report only; de and function on a dataset with counts are not compared |
| 651 | 2D.new | New dataset in one file | SUPERSEDED | by the guide (plan decision 1: one generated page, not `docs/new_dataset.md`): conf/dataset_template.toml, `ctx guide`, `ctx run --check` (tests/test_pipeline.py:125, :144), loader problems listed together (src/ctx/dataset.py:345) |
| 652 | 2D.cohort | Cohort bundles, `--cohort`, availability grid | BUILD | L, Pkg F. No conf/cohorts, no `--cohort` in cli.py |
| 653 | 2D.img | `figures` image on GHCR | RYTIS | Q9; the Dockerfile target is BUILD S in Pkg F |
| 654 | 2D.spec | FigureSpec, `ctx figure`, `source`, `verify-figure`, Reproduce panel | BUILD | L, Pkg F; its CI Docker step waits for row 653 |
| 655 | 2D.contour | Two-knob and power-map templates with click-to-set sliders | BUILD | M, Pkg F. Density grid and HDR levels exist (src/ctx/viz/contour.py:46, 66) |
| 656 | 2D.ts | Figure files without a timestamp | BUILD | S, Pkg F. The Step-1 volcano PDF holds `/CreationDate (D:20260930122634+02'00')` and its SVG `<dc:date>2026-09-30T12:26:34`; `fig.savefig` passes no metadata (src/ctx/viz/save.py:237). Liberation Sans half DONE: containers/analysis.Dockerfile:52 installs fonts-liberation |
| 657 | P2S.4 | M5, M6, hub text into content; plain contrast labels | BUILD | L, Pkg D. content/de.toml holds 18 lines and says the rest is still in de/page.py and de/run.py |
| 658 | P2S.7 | ss_ndn VST trend fallback stated | BUILD | M, Pkg D. `deseq.pca` records no trend (src/ctx/de/deseq.py:280) |
| 659 | P2S.8 | A room per contract for vertical axis titles | RYTIS | Q4 |
| 660 | P2S.5 | Selectable FDR on M6 | BUILD | M, Pkg D; builds the `threshold` kind |
| 661 | P2S.6 | A control on every figure | SUPERSEDED | Audit approved (D-051 amendment; docs/plans/2026-10-01_figure_controls_plan_approved.md); the builds are WO5 to WO10, M4.14, P7.10, M9.8 |
| 665 | 2D.knobs | QC exclusion knob in M5; volcano raw or shrunk x | BUILD | M, Pkg D. `ctx de` options are min_count and lfc_shrink only (src/ctx/cli.py:332-340). Default off, so no registered number moves |
| 668 | GS.2 | Entry pages | BUILD | M, Pkg E. No `index.html` written anywhere in src/ctx, build_kit.py or check_kit.py |
| 669 | GS.3 | Unthinned volcano | BUILD | M, Pkg D. Columnar rows exist (src/ctx/report/figures.py:49, e731c57); the volcano still keeps only raw p < 0.1 (src/ctx/de/run.py:245) |
| 670 | GS.4 | M5 picker | BUILD | L, Pkg E; `ui.toml [select]` exists (line 153) |
| 671 | 2D.M5contour | M5 contour figures | BLOCKED | on row 655 |
| 672 | 2D.review | Phase 2D review | BLOCKED | on every Phase 2D row |

### Phase 3

| Line | Id | Item | Status | Evidence, or what remains |
|---|---|---|---|---|
| 676 | 3.subset | Instructor-side subset on Zenodo | RYTIS | Q10. Done: P3.0 ran, 8 runs in ~/ctx_tier1/data/raw/tier1_teach, chr22 indices built by M3 (M3.9 ticked, D-041). Zenodo hosting remains |
| 677 | M1a/M1b | Accession routes, ID hierarchy, budget; exercises | BUILD | M, Pkg G. No docs/accessions.md; written as page text in content (guide plan decision 1), not a docs file |
| 678 | 3.registry | QC registry | DONE | M2.2 ticked; tests/test_qc_registry.py:23, :70, :155 |
| 679 | M2a | M2a on the laptop subset | DONE | TODO.md:682 ticks the same work; planted adapter and quality caught (tests/test_qc_reads.py:96); a planted sex swap flagged (tests/test_qc_checkpoint.py:222-233) |
| 680 | M2b | QC checkpoint | SUPERSEDED | by M2.3 and M2.8 (Step-1 count checkpoint, flags only, ticked) and M3.12 (the align stage) |
| 681 | 3.M3 | M3 on the Track A subset | DONE | M3.6, M3.9 ticked; D-041 (750.55 s, 1,599,750,144 B); the STAR-on-Linux half is M3.10 |
| 683 | 3.fullref | M3 full-reference comparison | SUPERSEDED | by M3.11 |
| 684 | 3.M4 | Salmon + pytximport; Venn | SUPERSEDED | by Priority 1: D-054 (ctx.quant.tximport, not pytximport); the cross-method comparison is M4's F5 after M3.11 and Q2 |

### Phase 4

| Line | Id | Item | Status | Evidence, or what remains |
|---|---|---|---|---|
| 688 | 4.step3 | Step-3 cohorts, TOIL tumour against GTEx | RYTIS | Q6. Done: GDC acquisition for four projects (P4.2, P4.3 commits fb30193, d7a0998, bdcf207), paired M5 (conf/datasets/tcga_brca.toml:111-119, e731c57), 12 men flagged (D-055). Open: TOIL and GTEx download, a PAM50 source, adipose composition |
| 689 | 4.cassetta | Breast Step-1 Cassetta; lung scout | RYTIS | Q8. Counts done (P4.1 56087ea, P4.4 b347cfb); no dataset file for breast_mac or breast_mono in conf/datasets. Lung scout is research, BUILD S after Q8 |
| 690 | P4.0 | Tissue downloads | DONE | breast_mac and breast_tissue counted and checked (56087ea, bdcf207); TCGA counts through GDC (P4.2, P4.3) |
| 691 | 4.step2 | TCGA-LAML or DLBC | RYTIS | Q8 |
| 692 | 4.upset | Pan-cancer Venn and UpSet | BUILD | M, Pkg C2: tumour_normal DE exists for BRCA, KIRC, LUAD, LUSC (e731c57) |
| 693 | 4.M7s3 | M7 Step-3 | RYTIS | Q6 for ccRCC, LUAD, LSCC tumour-against-NAT downloads. Breast subtype concordance and ERBB2 are F4 and F5 of Priority 5 (BUILD, Pkg C1); ERBB2 genes in conf/datasets/cptac_brca.toml:74-75 |
| 694 | 4.M7s1 | M7 Step-1 against Rieckmann | DONE | P7.6, P7.8; 849a8f7; D-057 |
| 695 | 4.M8 | M8 biotypes and TCGA miRNA targets | RYTIS | Q6 for the protein footprint (CPTAC-2 breast miRNA and isoform files). Rest DONE: 2b2246f, e202a43, D-058 |
| 696 | 4.spike | External-data spike | DONE | m9 plan section 1 (CORS, line 11), D-056; the plan's own "no Pyodide run" (line 163) moves to M9.9 |
| 697 | 4.M9 | M9 clients, cards, worked examples | BUILD | S, Pkg C1. Done: 7e6e1af, D-056, tests/test_genecard.py:207, :291, :319; qPCR pair src/ctx/genecard/annotation.py. Not done: breast and kidney worked-example cards (only ~/ctx_tier1/outputs/11_genecard exists) |
| 698 | 4.M10 | M10 breast, then kidney and lung | SUPERSEDED | by Priority 2: D-055; M10.10 to M10.12 |
| 699 | GS.5 | M6 live ORA | BUILD | L, Pkg E. No `_ora.js.j2`; `--genes-from` only on clinical and genecard (cli.py:421, 478) |
| 700 | GS.6 | Open beside, `postMessage` | BUILD | M, Pkg E. No `postMessage` in templates |
| 701 | GS.7 | M10 selection | DONE | bc7426b, e415097; `templates/_surv.js.j2`; `ctx clinical --genes-from` (cli.py:421); tests/test_survival.py:252, :312, :324; tests/test_clinical_page.py:32. Remainders in Priority 2 |
| 702 | GS.8 | M10 for the Step-1 root | BLOCKED | on row 652 (the pancreas cohort bundle) |
| 703 | GS.9 | M9 selection | SUPERSEDED | by M9.8 in Priority 4; built in Pkg E |
| 704 | GS.10 | Proof across four browsers | BLOCKED | on GS.2 to GS.9; Safari by Rytis |
| 705 | 4.ack | GDC acknowledgement on every GDC page | BUILD | S, Pkg C2 for the wording and M7, Pkg D for M6. Carried by de/run.py, ncrna/page.py, clinical/page.py, genecard/page.py; absent from function and proteo; the sentence names no project accession and is "UNVERIFIED ... written from memory" (src/ctx/report/content/ui.toml:145-149) |
| 706 | 4.ER | ER+ against ER- under four endpoints | DONE | D-055; tests/test_clinical_table.py:155; TODO.md:440 (p 0.086, 0.0095, 0.0004, 0.0018 against Liu's) |
| 707 | D28 | Treatment dataset | RYTIS | Q8 |

### Phase 5 and 6

| Line | Id | Item | Status | Evidence, or what remains |
|---|---|---|---|---|
| 711 | 5.lite | `ctx lite-build`, notebooks, Pages | BUILD | L, Pkg F; no notebooks folder, no `lite-build` command. Deploying to GitHub Pages is Q10 |
| 712 | 5.cards | Cards, prompts, quizzes, glossary, mistakes | DONE | `ModulePage` requires every section non-empty (builder.py:299-300); tests/test_cli.py:197, :203; glossary D-048 |
| 713 | 5.repro | Reproduce panels everywhere; Zenodo per cohort | BLOCKED | on row 654, and Q10 for Zenodo |
| 714 | 5.ci | Notebooks headless in CI | BLOCKED | on row 711 |
| 715 | 5.teacher | Teacher guide with timing | BUILD | M, Pkg G, from the measured Reproduce times |
| 719 | 6.img | `full` image on GHCR, Apptainer, smoke tests | RYTIS | Q9 |
| 720 | 6.bench | Runtime and RAM in reports | DONE | every page's `Reproduction` (builder.py:203-213) |
| 721 | 6.pilot | Pilot, v1.0, Zenodo DOI | RYTIS | Q10 |

### Verification

| Line | Id | Item | Status | Evidence, or what remains |
|---|---|---|---|---|
| 727 | V1 | Every module end to end in CI under 5 min | BUILD | Pkg A; CI red, 1,171 s |
| 728 | V2 | Offline, no http script, under 15 MB | DONE | tests/test_report.py:160, :179; the kit's link rule (D-047) |
| 729 | V3 | Figure variants; captions carry test or "no test" | BUILD | S, Pkg A. Files: viz/save.py:233-243. No test reads every caption for its test or "no test" |
| 730 | V4 | README and naming lint for every output | BUILD | S, Pkg A. Per module only (e.g. tests/test_process_run.py:83-97, test_proteo_page.py:44); one test over every end-to-end root fixture |
| 731 | V5 | gather fails loudly; no-md5 reported | BUILD | S, Pkg A. Failures tested (tests/test_gather.py:59-165, test_curation.py:229); `no_md5` exists (src/ctx/gather/fetchngs.py:135) but no test reaches it |
| 732 | V6 | QC registry and MultiQC parity | DONE | tests/test_qc_registry.py:23, :52, :70, :155; tests/test_multiqc_io.py:55 |
| 733 | V7 | Planted problems caught | DONE | tests/test_qc_reads.py:96; tests/test_qc_checkpoint.py:222-233 |
| 734 | V8 | Gene cards offline, degrade, NPTX2 | DONE | tests/test_genecard.py:291, :319, :345; D-056 (2.49 µg/L, 335.9 nCPM) |
| 735 | V9 | KM, log-rank, Table 3, risk table, corrected p | DONE | tests/test_survival.py:41, :79, :252, :312; tests/test_clinical_table.py:137, :155; tests/test_clinical_page.py:136. Oracle is R survival 3.7.0 and statsmodels, not lifelines (D-055) |
| 736 | V10 | Contour figure checks | BLOCKED | on rows 655, 671 |
| 737 | V11 | Reproduce panels run in CI | BLOCKED | on rows 653, 654 |
| 738 | V12 | `ctx run` byte for byte on a second file | BUILD | S, Pkg F, with row 650 |
| 739 | V13 | Cohort bundles; Step-3 sanity rules registered | RYTIS | Q5. tcga_*.toml carry only `[ncrna.sanity]`; the mRNA rules were never registered and M5 has already run on these counts |
| 740 | V14 | Biology sanity, Step-1 and Step-3 | BLOCKED | Step-1 DONE (D-046 R1, R2, PDAC recorded); KIRC waits for Q5 |
| 741 | V15 | Notebooks reproduce the CLI | BLOCKED | on row 711 |
| 742 | V16 | M7 sanity | DONE | conf/datasets/cptac_brca.toml K0 to K3 (lines 110-140); D-057; K0 records, not gates, because Krug prints no median |
| 743 | V17 | M8 sanity | DONE | N1 PASS (TODO.md:441); KIRC miR-210-3p +3.79 (D-058) |
| 744 | V18 | Paired against unpaired | DONE | D-046 R7: 6,290 paired, 4,570 unpaired |
| 745 | V19 | Notebook headless; cards match defaults | BLOCKED | notebook half on row 711; card half DONE (tests/test_cli.py:203) |
| 746 | V20 | Peak RAM and wall time within budget | BUILD | M, Pkg C2: BRCA `ctx ncrna` read 2.30 then 3.77 GB on the same inputs (TODO.md:441), cause not established |

### Phase 2 implementation-order sections

| Line | Id | Item | Status | Evidence, or what remains |
|---|---|---|---|---|
| 1625 | P2A.3 | Cluster download and rnaseq on 210 runs | DONE | D-046: job 47485533's counts on the Mac, 210 samples joined, 634 files hash-stable |
| 2003 | P2C.4 | Track B bundle in an empty root | BUILD | M, Pkg F; measure against the 100 MB cap first, over it goes to Rytis |
| 2019 | P2C.5 | Close | BUILD | M, Pkg A. D-007 to D-010 do not exist in knowledge/decisions.md; the fresh-clone under-5-min run waits for M3.14 |

## BUILD packages

Waves run in order; packages inside a wave run in parallel on disjoint files. A hot file has one owner per wave and passes on at that owner's merge.

| Hot file | Wave 1 | Wave 2 | Wave 3 | Wave 4 |
|---|---|---|---|---|
| src/ctx/cli.py | C1 | D | E | F |
| src/ctx/report/content/cards.toml | C1 | D | E | F |
| src/ctx/report/glossary.yaml | C2 | D | E | F |
| src/ctx/dataset.py | C2 | D | E | F |
| tests/test_cli.py | A | D | E | F |
| src/ctx/report/templates/_page.js.j2 | none | D | E | F |
| src/ctx/report/content/ui.toml | C2 | D | E | F |

Wave 1 packages add no CLI option, card or glossary entry unless they own that file.

**A. CI green and the suite under 5 minutes** (wave 1; Phase 2M, Verification, P2C.5). Size L.
- Items: M3.10, M3.14, V1 (727), V3 (729), V4 (730), V5 (731), P2C.5.
- Files: .github/workflows/ci.yml, tests/browser.mjs, tests/conftest.py, src/ctx/runtime.py, tests/test_runtime.py, tests/test_io.py, tests/test_cli.py (threads test only), tests/test_process_align.py, tests/test_process_count.py, tests/test_gather.py, a new tests/test_outputs_contract.py, knowledge/decisions.md (D-007 to D-010).
- Root causes to reach: the runner's AppArmor userns rule (fix in the runner step, not `--no-sandbox` in the page tests); git identity inside the test's own repo; the BLAS pool count and the Linux parent-peak reading, both not yet explained.
- Proof: a green CI run with `-rA` or junit showing M3.10's test passed on ubuntu; the same test on the Mac; `/usr/bin/time -l pixi run test` under 300 s; mutations: drop the runner fix (72 errors return), plant a caption with no test (V3 fails), a table with no README (V4 fails), a run with no md5 file failing instead of reported (V5 fails).

**B. QC bands and definitions** (wave 1; work order 3, Phase 2Q). Size L.
- Items: P2Q.9, P2Q.10 (definitions), P2Q.3 composition scatter, M3.12, L.6, `ctx.qc.identity` (the sex check moved from de/run.py:607, CONVENTIONS' after-the-record rule) so WO5's sex view can follow.
- Files: conf/qc_metrics.toml, src/ctx/qc/** (new identity.py), src/ctx/viz/strip.py, src/ctx/gather/rnaseq.py (`COPIED_FILES` with the FastQC zips), src/ctx/report/content/qc.toml, src/ctx/report/glossary_lint.py, tests/test_qc_*.py, tests/test_glossary.py, tests/test_multiqc_io.py. Ordered exception: the two-line import change in src/ctx/de/run.py, merged before D branches.
- Proof: the registry refuses a band out of order and an unknown key; strips draw plot_bands; FastQC per-position curves; `_show` never prints "100 is below 100"; full-cohort FastQC flags counted and recorded; top_gene_share median 5.34 % falls to 3.22 % as the proposal measured; Step-1's 92 recorded metrics reproduce; "OS" no longer suggests survival_endpoints. Mutations per rule.

**C1. M4, M7, M9 remainders on real data** (wave 1; Priority 1, 4, 5). Size L.
- Items: P1 open parts (heatmap triangle, strip colours, M4.4, M4.11, M4.15, browser), P5 (F4, F5, aliquot column names, contours, P7.11), P4 (cell-line search, COXFA4L2 resolution, browser), row 697 (breast and kidney worked-example cards).
- Files: src/ctx/quant/**, src/ctx/proteo/**, src/ctx/genecard/**, scripts/python/preview_quant.py, content/quant*.toml, proteo*.toml, genecard*.toml, src/ctx/cli.py, content/cards.toml, their tests.
- Proof: M4.11 inside Track A's budget with `/usr/bin/time -l`; 7 cell types with 7 distinct colours; F4 calls by D-012's rule; K1 to K3 still pass; budget tests at real size; browser at 1280 and 375 px both themes; mutations per new function.

**C2. M10, M8 and TCGA remainders** (wave 1; Priority 2, 3, Phase 4). Size L.
- Items: P2 (Dunn's through proteo.measure's one definition, preview Cox and KM retired, box labels, M10.10, M10.11, the kit), P3 (F1 by biotype, LUSC histone rule, box order determinism), V20 (746), 692 UpSet, 705 (wording read against GDC's own page, the project accession named, added to M7), K.6 release mode.
- Files: src/ctx/clinical/**, src/ctx/ncrna/**, src/ctx/gather/gdc*.py, conf/datasets/tcga_*.toml (clinical blocks), src/ctx/report/content/clinical.toml, ncrna.toml, ui.toml, glossary.yaml, src/ctx/dataset.py, scripts/python/build_kit.py, check_kit.py, scripts/python/preview_clinical.py, their tests. Not src/ctx/proteo (C1): Dunn's is imported, not copied.
- Proof: Liu Table 3 counts for KIRC, LUAD, LUSC; two runs of BRCA `ctx ncrna` give the same box order and their RSS recorded with the cause; the release kit refuses a made-up root; mutations.

**G. Guide and course text** (wave 1; Phase 2G, 3, 5 docs). Size M.
- Items: G.6/U.12, G.7/U.13, M1a/M1b content (677), teacher guide (715), CONVENTIONS amendments the D-entries already require (lifelines boundary, D-055; pytximport, D-054, both still in CONVENTIONS.md:29, :33).
- Files: CONVENTIONS.md, src/ctx/report/content/guide.toml, gather.toml, a new content file for the teacher guide and its page wiring in src/ctx/report/course.py, tests/test_guide.py, tests/test_prose.py. Marks only to existing glossary entries.
- Proof: `pixi run lint-prose` and the glossary lint clean; every Copy clicked at 1280 and 375 px in both themes; the teacher guide's times read from Reproduce blocks, none typed.

**D. M5 and M6 text, the threshold kind, DE knobs** (wave 2; work order 4, Phase 2D, 2P, 2Q). Size XL, two merges.
- D1 items: P2S.4, P2P.3, P2S.7, the M6 GDC acknowledgement. D2 items: P2S.5 (threshold kind), GS.3, the QC exclusion knob and the volcano's raw or shrunk x (665), P2Q.7.
- Files: src/ctx/de/**, src/ctx/function/**, src/ctx/report/hub.py, builder.py, src/ctx/viz/annotate.py, content/de.toml, function.toml, conf/datasets/tier1.toml (plain contrast labels only), plus wave 2's hot files.
- Proof: the 92 recorded Step-1 metrics reproduce value by value (D-046); the flat-dispersion fixture's caption names the mean-based trend; the threshold kind's JS twin equals the Python definition; the unthinned volcano page measured under 15 MB; every planted mutation caught.

**E. Gene selection and figure controls** (wave 3; work order 5 to 10, GS). Size XL, merged module by module.
- Items: WO5 (subset and view kinds, M1 and M2, item 2's two leftovers), WO6, WO7, WO8, GS.2, GS.4, GS.5, GS.6, M9.8 (GS.9), M4.14, P7.10, M8 F2, F3, F7 controls, M10's Cox screen on the live selection and live captions.
- Files: src/ctx/report/templates/** (new `_ora.js.j2`), src/ctx/report/genes.py, each module's page.py and figures.py for its controls, scripts/python/build_kit.py and check_kit.py for entry pages (after C2), plus wave 3's hot files.
- Proof: each control's JS twin against Python; entry pages reached by the kit's link walk, a missing one refused; live ORA within 1.52e-10 of `lite.ora` (D-049); a real lasso and box at 1280 and 375 px; the M9 page at the byte cap with 57 cards.

**F. The reuse layer and JupyterLite** (wave 4; Phase 2D, 5, P2C.4). Size XL.
- Items: 650 and 738 (de and function in the byte-equal proof), 652 cohort bundles, 654 FigureSpec, 655 contour templates, 656 timestamps, the `figures` Dockerfile target (653's build half), P2C.4, 711 lite-build and notebooks 05 to 08.
- Files: src/ctx/pipeline.py, tests/test_pipeline.py, conf/cohorts/, a new src/ctx/figures/ package, src/ctx/viz/save.py, contour.py, containers/analysis.Dockerfile, notebooks/, plus wave 4's hot files.
- Proof: two renders of one figure byte-equal; the cohort schema refuses a bad bundle, every error listed at once; the Reproduce panel run in a clean container on fixtures; the Track B bundle under 100 MB compressed run from an empty root with the repo unreachable.

## Questions for Rytis, recommended answer first

1. **Cluster: M3.11 with M4's Q2, then M3.13 and P2Q.8.** Yes: an agent writes one block, the 8 trimmed runs aligned with STAR to GENCODE 50 and Salmon run on the same reads; M3.13 and P2Q.8 are planned from what it measures. Cost: Q2 is one job of 4 cores, 16 GB, 1 h (M4 plan); the STAR index's RAM is measured there, not assumed. Unblocks M4 F5 and F8 and M3's full-reference side.
2. **G.5 / U.11 rehearsal on Kebnekaise.** Yes, after Pkg G's CONVENTIONS line. Cost: your session.
3. **Docker memory on this Mac for U.10.** Raise Docker Desktop to 20 GiB of the 24 GiB. Now 12,525,977,600 B (11.67 GiB) against the 16.8 GiB Salmon index peak. Cost: less memory for the rest of the Mac while it runs.
4. **P2S.8 vertical axis titles.** Each contract's room computed from its own drawn axis height, replacing `AXIS_TITLE_SHARE = 0.55` (src/ctx/viz/style.py:139). Cost: some titles rewrap; figures re-measured.
5. **TCGA mRNA sanity rules (rows 739, 740).** Register them now, labelled NOT blind as `[ncrna.sanity]` was (conf/datasets/tcga_brca.toml:143), since M5 has already run on these counts; KIRC names NDUFA4L2 by v36. Then the P2Q.10 bands, after Pkg B fixes the definitions. Alternative: drop the rules.
6. **Downloads.** Approve CPTAC-2 breast miRNA and isoform files, 5,251,801 and 47,005,217 B, GDC open (M8's protein footprint, row 695). For SCAN-B (591,676,211 B, terms UNVERIFIED), TOIL and GTEx for row 688, CPTAC ccRCC, LUAD and LSCC for row 693, and the GTEx V11 junction file: approve an agent measuring bytes by HEAD and reading each licence, then ask again with the numbers.
7. **HPA's own tau.** Close it: v25's files carry none, and the page already labels tau as ours (D-056).
8. **Datasets not yet chosen.** Cassetta's breast monocytes and macrophages: an agent writes the dataset files with the paper's own comparisons and registers their sanity rules before `ctx de`, brought to you first. Step-2 (TCGA-LAML or DLBC) and the D28 treatment dataset: defer past the v1.0 pilot; cost, v1 has no Step-2 and no treatment effect.
9. **GHCR.** Yes: once Pkg F adds the `figures` target, dispatch the containers workflow with version 0.1.0, which pins digests in containers/digests.lock. CI jobs run again (2026-10-08), so billing no longer blocks. Unblocks K.5, rows 653, 719, 737.
10. **Release.** Ship the kit as a GitHub release asset now and hold Zenodo until the Step-1 redistribution terms are settled (the IDF names no licence). The release unblocks K.4's Windows and download checks, K.6, the newcomer test (G.8 / U.14), Pages for JupyterLite and the pilot.

## Not done

No pixi, prek or test run here; the CI log was read with `gh`. No browser check. Plan files under docs/plans are proposals as written and were read for their step lists, not ticked.
