# Decisions — Cancer Transcriptomics Teaching Suite

Analysis decisions that a reader could reasonably have made differently, with the
evidence that settled them. Not a changelog: record the ones where the outcome
would change if the decision changed.

Newest first. One entry per decision, and never delete one — a superseded
decision is marked superseded, with a pointer to the entry that replaced it.

Template for a new entry:

## D-001 — <one-line statement of the decision>

**Date:** YYYY-MM-DD · **Status:** decided | superseded by D-00N

**Context.** What made this a question.

**Options.** What was actually considered, including the one not taken.

**Decision.** What was chosen.

**Evidence.** The number, plot or run that settled it. A decision with no
evidence line is a preference — say so rather than dressing it up.

**Consequence.** What now depends on this, and what would have to be redone if
it were reversed.

---

## D-050 — Read-level pages show the full Step-1 cohort; the 8-run subset is an optional laptop exercise; made-up pages leave at the release

**Date:** 2026-10-01 · **Status:** decided

**Context.** Rytis asked why the course needs a teaching subset: at the end it should hold full real data and no made-up results. The read-level modules (M2's read stage, M3) were designed for Track A: 8 real Step-1 runs, reads capped, aligned to chr22 on a laptop within 15 min and 6 GB on 4 cores, because the 210 runs need a cluster (a full human STAR index needs about 32 GB; TODO.md, Track A design). The count-level modules (M1, M5, M6) already show all 210 samples. Made-up previews (P2V.2) give M3, M4 and M7 to M10 a page before their real data exist; a kit holds them only through `build_kit.py --made-up`.

**Options.**
- Read level: the full cohort on every page with the subset kept as an optional laptop exercise; the full cohort only, the subset retired, so no reader can rerun alignment on a laptop; or the pages on the subset, as designed.
- Made-up pages: out of every kit from v5, or kept in development kits and out of the release kit.

**Decision.** Rytis, 2026-10-01: the first of each. Every page a reader gets for M2's read stage and M3 shows all 210 Step-1 runs, aligned to the full GENCODE 50 genome on Kebnekaise. The 8-run subset stays as the laptop exercise, under Track A's budget, and its pages say it is the optional exercise. Development kits keep the made-up previews, marked Made-up; the release kit is built with no made-up root.

**Evidence.** None: a teaching choice by the course's author. Its cost is stated, not measured: 210 full-genome alignments on the cluster, with RAM, wall time and disk still to be measured on M3.11's 8-run pilot.

**Consequence.** M3.6 is unchanged, because `ctx process` and its page run on any dataset. M3.9 still measures the laptop exercise. M3.11 becomes the pilot of M3.13, the full-cohort M3; P2Q.8 runs M2's read stage on the full cohort; K.6 builds the release kit with no made-up root. Students still download only the subset, never the full FASTQs.

---

## D-049 — One gene selection drives M6, M9 and M10: carried in links, never in browser storage

**Date:** 2026-10-01 · **Status:** decided (Rytis's three answers, 2026-10-01); built by GS.1 to GS.10

**Context.** Rytis, 2026-09-30: "The M6, M9, M10 should be flexible and interactive. There should be an option to select a list of genes (maybe selecting them from a table, or in a Volcano plot - like in "Figure 2 Samples on the PCA, beside their labels and QC status") which would then automatically change the M6, M9, M10 shown results". Pages are separate self-contained files opened from disk, at most 15 MB each.

**Options.** Three designs were prototyped and measured in Chrome 138, Edge 154, Chromium 152 and Firefox 156 (Safari untested): (A) the selection in the URL fragment of ordinary links between pages; (B) one selection in browser storage with a bar on every page; (C) one combined gene-explorer page per dataset. Two judges scored them; the plan (scratchpad `genesel/PLAN.md`, outside the repo; the GS items in TODO.md are the record) builds A with C's unthinned volcano, survival checks and symbol rules, and B's entry pages, genes-by-contrasts view and build-list notice.

**Decision.**
- A gene is its unversioned Ensembl number, with the sign of its raw log2 FC in the source contrast; a page holds one selection.
- The selection travels in the fragment `#v=1&src=<dataset>:<contrast>&n=<count>&g=<numbers>&end` (plain decimal, `-` for a down gene), decoded strictly against `n`, rewritten with `history.replaceState` on every change. One definition, `ctx.report.genes`, with its JS twin `ctxSelect`. The end mark, added at GS.1 after review: a link cut anywhere, inside its last gene number included, loses it and is refused as cut rather than read as another gene.
- Browser storage is not used: in Firefox at its default policy each local file has its own storage, and an open M6 tab kept a stale list with no warning. "Open beside" windows follow each change through `postMessage`.
- Entry pages (Rytis, decision 3): each selection-consuming module writes `outputs/NN_<module>/index.html` on every build, forwarding to its newest page with the fragment; the naming lint exempts that one name in module folders. The kit's link walk reaches them, and a listed module that was not built refuses the kit.
- M6 recomputes ORA live (Hallmark; direction and background knobs, the trap labelled), reproducing `ctx.function.lite.ora` to within 1.52e-10 relative in p on Step-1. GSEA stays static.
- M10 (Rytis, decision 1): a Step-1 gene list is scored on TCGA-BRCA tumours, labelled cross-tissue: 1,089 patients and 145 PFI events after CDR redactions. The score is one definition (`ctx.clinical.score`, JS twin) in a fixed order: log2(CPM + 1), int16 storage, GENCODE v36 by unversioned ID, the expression filter, z per gene, the signed mean, scaled to SD 1.
- M10 (Rytis, decision 2): Cox for a selection's score is fitted live in the page, univariable and adjusted for age and stage, validated against R and Python in tests; a monotone likelihood is refused. This diverges from the M10 plan's "Cox fits precomputed", because an arbitrary selection cannot be precomputed. The per-gene Cox screen stays precomputed.
- M9: evidence across all contrasts for any tested gene; full cards only for the build's card list, packed to the byte cap; a gene without a card shows the command that builds one.
- Every live result carries its test, n, adjusted p and effect size; the PDF and SVG show the build list, and the page says when the live selection differs.

**Not done here.** Safari, Windows paths and touch lasso are untested; tabs not opened through Open beside do not update; GSEA, Reactome and the other contrasts' volcanoes are not live; no correction across the lists a student tries, which the page counts and says.

---

## D-048 — Rytis's draft v3.12 merged: the glossary, term tooltips and the glossary page; three v3.11 prose rules applied

**Date:** 2026-09-30 · **Status:** decided (the draft's rules; the defaults below set by the lead of the merge, open to Rytis)

**Context.** Draft v3.12 (`TEMP/glossary.yaml`, 179 entries, 8 marked `confirm`; the glossary sections of `TEMP/CLAUDE10.md` and `TEMP/todo10.md`) adds a glossary and hover tooltips on jargon, planned as Phase 2L. An audit of the v3.11 and v3.12 drafts against the repo found the glossary unbuilt and three v3.11 rules unapplied: no rhetorical questions, the figurative verbs in every form, and Reproduce comments at most 12 words per line.

**Decision.**
- `src/ctx/report/glossary.yaml` is the draft's file with the edits below. `ctx.report.glossary` loads it (PyYAML), looks an entry up by id, term or alias in any case, and `explain(term)` prints one entry in a notebook. `ctx.report.marks` holds the one definition of the markup; `ctx.report.glossary_lint` is the lint (`python -m ctx.report.glossary_lint`, and `tests/test_glossary.py`).
- Marks: `[[term]]` shows the words as written, `[[term|display]]` shows the display words, `[[term!]]` gives a later mention its icon. They are allowed in a page file's sections, questions, pitfalls and figure captions and in the cards, the text `explain_terms` renders; anywhere else the lint refuses them, since that text is printed without the filter. The first mention per section carries the icon.
- Off the page a mark is its display text: `annotate.place_caption` wraps a caption's display text into the PDF and SVG, and `--help` prints a card's. A figure refuses a mark anywhere but its caption (its title, an axis label, a hover text, any drawn text) before any file is written, since the plotly figure and the page heading would print the brackets (`ctx.viz.save`). A mark inside backticks is code: `marks.pieces()` is the one reading of the markup, so the lint, the tooltips a page embeds and the renderer agree.
- The glossary page lives at `outputs/09_report/<date>_report_glossary_v1.html` in each ctx root, beside the hub. `builder.write` writes it before any page that marks a term, and before the hub, unless the same bytes are there; the page links it as `../09_report/<name>` from the top bar and from every tooltip. Every constraint holds: the name passes the naming lint; the file exists before the page that links it; a kit copies it because a page links it (and `current()` also takes the newest); `START.html` links the newest glossary of any root. A page with marks written outside `outputs/NN_<module>/` keeps its tooltips and links no glossary: from Phase 2D, `ctx figure` run from a Reproduce panel writes into the reader's own folder (CONVENTIONS.md, Outputs), where no glossary sits, and a refusal there left the PDF and SVG written and the page missing.
- The tooltip script is `templates/_glossary.js.j2`, inlined like the other page scripts: hover (mouse only) and focus open; a click or tap pins it, or closes one that is pinned or was opened by focus; Escape, an outside click, blur and a resize close it, and Escape gives focus back to the button; one is open at a time; the panel has `role="tooltip"`, the button `aria-expanded` and, while open, `aria-describedby`; the panel sits in the term, so Tab reaches its link, but is `position: fixed`, placed in viewport coordinates and moved with its term on any scroll, so no scrolling card table clips it (no ancestor of a term may set transform, filter, backdrop-filter, perspective, contain or will-change; a test reads the stylesheet); it opens below unless it would leave the page's own box and fits above, and moves sideways into that box, measured as `documentElement.clientWidth` and `clientHeight`, since `innerWidth` counts a scrollbar and on a phone grows with the open panel. Keyboard focus keeps a hover-opened tooltip open when the pointer leaves; a press on the panel's text keeps it open (the panel takes focus, `tabindex="-1"`); a press that ends in no click does not stop the next keyboard focus from opening a tooltip. A page embeds only the entries it marks.
- The prose lint (`ctx.report.prose`) now bans know, bring, want, tell and arrive in all 21 forms, -ing included, and refuses a `?` anywhere in a content file outside the Questions (titles, headings, leads, lists and prose). A strand end (5', 3') is not a number. A sentence may start in lower case (fastp, lncRNAs, nf-core, pytximport): a stop before a lower-case word ends one, unless it closes e.g., i.e., vs., cf. or et al.
- The glossary schema requires `aliases` (a list, which may be empty), as the draft's header marks it, and refuses a blank required field. At a release tag the lint fails on `confirm` and on unsourced numbers; only a `v*` tag is a release, from GitHub's ref or from `git describe`.

**Divergences from the draft, and why.**
1. The icon uses the link colour, #1D4E89 light and #8CB8F0 dark, not #0072B2 and #56B4E9: D-024 keeps data colours off page furniture.
2. The tooltip's link reads "Glossary", with no arrow: the pages carry no arrows.
3. A number in `short`, `long` or `typical` of an entry with no `source` is a warning now and an error at a `v*` release tag, like `confirm: true`. No source was invented: 33 entries give a number with no source.
4. Definitions follow the prose rules, banned words included, except the typed number: a definition states numbers and names their source.
5. The loader is `ctx.report.glossary`, not `ctx.glossary`; `explain()` prints.
6. The lint is `python -m ctx.report.glossary_lint`, not `ctx lint-glossary`: a `ctx` subcommand must be a course module (D-045's reason for `lint-prose`). The `pixi run lint-glossary` task is added at the merge, with the PyYAML pin; this branch does not touch `pixi.toml`.
7. The page is `outputs/09_report/<date>_report_glossary_v1.html`, not `reports/glossary.html`: pages live in dated module folders here (D-047's kit copies what they link).
8. `tooltip.js` in `assets/` is `templates/_glossary.js.j2`: every page script is a template.
9. The Reproduce comments rule adds nothing: the Reproduce section prints one command built from the command line, with no comment text, and the Phase 2D panel that would hold comments is not built.
10. Question-form group leads on the gene card page (3) became imperatives, and "told apart" (3 page sentences) became "distinguished".
11. `[[term]]` shows the words as written, which the draft left open: the writer marks the natural word, and `[[term|words]]` covers any other.

**Edits to the draft's glossary** (before and after of each are in the merge report):
- Three aliases that equal another entry's term, so one word named two entries: "per-base sequence content" left `random_hexamer_bias`, "Salmon" left `pseudoalignment`, "tissue source site" left `tcga_barcode`.
- For the repo's cohort guard (`tests/test_repo.py::test_no_cohort_identifiers`): the aliases TCGA-BRCA, TCGA-KIRC, TCGA-LUAD and TCGA-LUSC left `tcga`, and `tcga_barcode.short` writes the barcode as TCGA-<site>-<patient>-01A, not TCGA-XX-YYYY-01A.
- 55 definition fields of 49 entries differ from the draft: 7 fields of the 4 entries corrected to their sources (below), 1 for the cohort guard, and 47 reworded to the prose rules (sentences over 25 words, paragraphs of 4 sentences, intensifiers, figurative verbs), meaning kept. A first version of this entry said 5 rewordings were made only because the lint read a sentence starting in lower case as part of the one before; review measured 3 (adapter_content.long, polyg.long, expression_filter.long), and the fault was the lint's sentence splitter, now fixed. The first two are Rytis's text again; expression_filter.long keeps his lower-case start but not "the first casualties", figurative language about data. mapping_rate.typical and tx2gene.long exceed 25 words under either splitter; they keep his "nf-core/rnaseq" and "pytximport" starts, and tx2gene.long keeps "changes in abundance rather than in isoform length". Where "known" (a banned form) was replaced, splice_junction.long reads "confirm annotated isoforms" and pseudogene.long "a recognised source of noise", not "existing" and "common", which changed the meaning.
- The 8 `confirm` entries, checked against sources read on 2026-09-30 by the merge's source research, each quoted: per_sequence_quality, n_content and gc_content against FastQC 0.12.1's `Configuration/limits.txt` and module documentation (and Lander et al. 2001 for 41 % GC), mapping_rate against Conesa et al. 2016 and nf-core/rnaseq 3.27.0, are sourced as written. Four are corrected to their sources: q30 (Illumina's specification pages give 70 to 90 % above Q30, not 75 to 90 %, measured on the PhiX control; Q30 means an error chance of 0.1 % or less, not below it), overrepresented (FastQC lists every sequence over 0.1 % of reads, tracked in the first 100,000, with no chance model), small_rna_seq (the quoted 1 to 2 M reads was mapped reads correlated in one sample, not a depth for differential expression; Campbell et al. 2015 found fold changes reproducible at about 5 M reads per sample) and star (STAR 2.7.11b: a read over 10 loci gets no alignment, matched bases and score must each reach 66 %, 16 GB of RAM at least and 32 ideally). Each now gives its source, and none is flagged `confirm`.

**Evidence.**
- The lint on the edited file: 0 errors; 218 warnings: 33 numbers with no source, no `confirm` flag, and 185 glossary terms left unmarked in the page files; no section's first mention of a term is bare where a later one is marked (4 were, and are moved: lncRNA and BH in M8's Questions, the annotation GTF in M3's Input, shared peptides in M7's Pitfalls). The prose lint: no violation.
- 154 marks in the 9 page files and 6 in the cards, the first mentions that read naturally. The course preview rebuilt from the final tree (`scripts/python/build_preview.py`, 93.5 s, peak RSS 618,348,544 bytes): 10 of the 11 module pages carry icons (M6's text is still in Python and its cards mark nothing), 25 of 69 pages embed tooltips and every glossary link and anchor resolves, no page shows a raw mark and no figure's plotly text holds one, the naming lint passes on all 520 outputs, and none of the 112 PDFs or 112 SVGs holds a mark. A kit built from it (`build_kit.py`, 526 files) passes `check_kit.py`, and its `START.html` links the glossary.
- Sizes: the glossary page 285,752 bytes (each entry carries its search text); `ctx setup`'s page 76,307 bytes on HEAD's code and 89,009 on this branch (the script, the styles and 3 terms).
- In Chromium (the Claude Browser pane) on the preview's pages: hover opens and leaving closes; Tab opens, Tab reaches the link, Escape closes and returns focus, Enter toggles; a click pins and an outside click closes; near the bottom the panel opens above and stays inside the view; the button measures 24 by 24 px; 375 px and dark theme; the glossary page's anchors, search and topics. Taps in the pane arrive as mouse clicks, so touch is tested under node. After review, on the rebuilt preview: at 375 by 812 px the FASTQ tooltip spans 16 to 367 px and the page's scrollWidth stays 375, and qc_bands, its phrase at 778 to 801 px, opens above at 606 to 771; at 1280 by 800 px the M5 card tooltips draw past their table's edge, a hit test at their corners finds the panel, the table grows no scrollbar, and the right edge stops 8 px inside the 1265 px page, not under the scrollbar; the panel follows a page scroll; after a drag off the MADs icon, Tab opens qc_bands and then mann_whitney; a click on a pinned tooltip's text keeps it open; the heading is the text colour in every M8 section, Pitfalls included; the glossary search finds 5, 10, 2 and 0 entries for "also", "source", "typical" and "verified".
- Tests: `tests/test_glossary.py` 45, `tests/test_tooltip.py` 17 and `tests/test_prose.py` 37 (14 before D-048); `test_cli.py::test_hub_links_resolve_after_move` now checks the hub's glossary link, and `test_report.py::test_copy_has_no_dashes_or_arrows` reads the YAML too. 82 planted mutations, each applied once to a copy of the tree and restored: every one fails the test named for it (the hub mutation also fails the `test_cli` test). The full suite after review: 632 passed, 1 skipped, 7 deselected; 631 in one run of 369.6 s, peak memory footprint 584,058,800 bytes (`/usr/bin/time -l`), where `test_text.py::test_an_edited_copy_changes_the_page` found no `ctx` on the PATH outside pixi, and passed rerun with the environment's `bin` on the PATH.

**Consequence.** L.5 is done: the 8 `confirm` entries carry sources. Before a release tag, each of the 33 numbers the lint lists needs its source. Every new page text is marked as it is written. The lint runs at a release only where a `v*` tag is checked out; `ci.yml` runs on branches.

---

## D-047 — The student kit: the newest outputs and what their pages link, START.html, one checker

**Date:** 2026-09-30 · **Status:** decided (Rytis's K.0; the rules below are the builder's)

**Context.** Rytis's experience for a student (Phase 2K): download the course's generated results and reports, open the folder, double-click START, and reach every report from one page. The kit is a release asset of the private repo, and on Zenodo with a DOI at the first session (K.0).

**Decision.**
- `scripts/python/build_kit.py --out <kit> --real <folder>=<ctx root> --made-up <folder>=<ctx root> [--zip]` builds it in `<kit>.partial` and renames it only when every check passes; the kit, its partial folder and its zip are never written over.
- What ships from each root's `outputs/` (resolved, so a symlinked folder works): the newest dated copy of every output in each module folder, plus every file the shipped pages link to, followed through linked pages (MultiQC's report sits in an undated subfolder). Undated files no page links stay behind.
- Refused before anything is written: reads, alignments and their indexes by extension; any file over 50,000,000 bytes (the largest output on 2026-09-30 is Step-1's counts CSV, 39.5 MB); any link in a source page that does not resolve or that fetches from off the machine. Every attribute a browser follows counts (`href`, `src`, `srcset`, `data`, `poster`, `xlink:href`, CSS `url()` and `@import`); a clicked `http`, `https` or `mailto` link is allowed, a `file:` link never.
- Checked after writing: the same link rule over every page of the kit, and every name against Windows: the full path after Explorer's "Extract all" into `C:\Users\<20 characters>\Downloads\<kit>\<kit>\` within 259 characters, no forbidden character or reserved name, no two names differing only by case.
- The checks are `scripts/python/check_kit.py`, standard library only, which the builder imports and copies into the kit: `python check_kit.py` in a copy checks every file against `MANIFEST.json`, every link and every name.
- The home page is `START.html` (`ctx.report.course`), built from the pages the kit holds: one card per module M0 to M10, one link per dataset's newest overview page, each marked Real or Made-up by its root, and each dataset's hub. M0 checks the computer, the same in every root, so it is listed from real roots only, and the count of modules on real data runs over M1 to M10. The page tells a reader whose links fail to unzip first. Text is `content/course.toml`, under the prose lint, which a test now runs.
- START is the HTML file, because a double-clicked `.command` is typed into the user's login shell, and on Rytis's Mac oh-my-zsh's overdue update prompt read its first character. `START.command`, `START.cmd` and `START.sh` stay beside it for readers who look for them, and open it. Amended 2026-10-01 (Rytis): the three scripts are gone. Each only opened START.html, which by then warned by itself when copied out of the zip (`kit.js`), and each brought what the page avoids: the reader's shell and its startup, and the prompts macOS and Windows show for a downloaded program. `check_kit.ps1` lost its START.cmd step; `tests/test_kit.py::test_the_page_is_the_only_way_in` fails if a script or any executable file returns to the kit.
- `MANIFEST.json` lists every other file with its size and SHA-256, the commit that built the kit, marked `-dirty` when any code path (`src`, `scripts`, `conf`, the pixi and pyproject files) has a change, an untracked file included, and each root's folder and kind.
- The made-up dataset is named "Made-up study" on every page (`ui.toml` `preview.dataset`, also its dataset file's label), not its slug.

**Evidence.**
- The first kit, built before the review from the uncommitted tree: 1,064 files, 785.1 MB, zipped 253.5 MB, in 17.8 s at 280 MB peak RSS; longest zip entry 128 characters, 193 with a 20-character user name after "Extract all", against 259; no forbidden character, reserved name, non-ASCII name or case collision.
- An adversarial review of that build (6 agents) confirmed 9 defects, all fixed here: the web-load check saw only `src`; a failed check left a half-written kit that blocked every rerun; M0 counted as real data; a symlinked `outputs/` refused every link; the home page gave no sign that the kit must be unzipped; the widened unused-key test passed three planted dead keys; nothing ran the prose lint; one test file was not formatted; the tests missed the newest-page rule, the size cap, broken-link refusal and the START target; untracked builder code did not mark the commit dirty. Five other claims were refuted.
- `tests/test_kit.py`, 14 tests; each of 10 planted mutations fails one (the newest-page rule, the cap, source-link refusal, the resolved outputs folder, M0 in the data count, the START target, manifest hashes, clicked links, the Windows length, the rename). Not tested: the dirty mark, which depends on git state.

**Not done.** K.4 on Windows (a PC with a checklist, by Rytis's choice) and the macOS download prompt, which waits for the release; K.5, blocked on CI billing and the Phase 2D image.

---

## D-046 — The Step-1 record: 7 of 8 registered rules pass; R4 fails on its monocyte condition

**Date:** 2026-09-30 · **Status:** decided (the record, as registered)

**Context.** The counts of nf-core/rnaseq job 47485533 reached the Mac after the last registration (3620021). Nine writers were editing the page code in the main tree at the time, so the record ran from a clean checkout of 891ab62, whose frozen files are identical to the registration's; the checkout's branch, `step1-record`, carries the record.

**Decision.** The verdicts stand as the hub decided them:

| Rule | Measured | Verdict |
|---|---|---|
| R1 | PERMANOVA R² of cell type 0.684; condition 0.089; batch 0.026 | PASS |
| R2 | p < 0.05: 6,370 against 2 × 488.8 (ndn_gcsf); 3,426 against 2 × 988.8 (ndn_hsct_early) | PASS |
| R3 | 6 genes at padj < 0.05 in the null split, against 1 % of tested, 140 | PASS |
| R4 | NDN: NES 3.18, FDR 0. Monocytes: NES 3.31, not below NDN's | FAIL |
| R5 | NES −2.44, FDR 0 | PASS |
| R6 | PERMANOVA p for batch 0.002 | PASS |
| R7 | 6,290 genes paired against 4,570 unpaired | PASS |
| R8 | 0 samples disagree with their recorded sex | PASS |
| PDAC | NES 2.48 adjusted, 2.26 within the long batch | RECORDED |

R4's first two conditions hold: HSCT neutrophils carry the interferon-alpha signature. Its third, monocytes less so, does not: in these bulk data the monocyte NES is 3.31 against the neutrophils' 3.18. The rule compared NES across two contrasts of different size (8 against 4 monocytes, 13 against 5 neutrophils), and NES is not calibrated across such contrasts; that is a limit of the rule as registered, stated here, not a reason to change it. The course teaches the failure beside the paper's claim.

**Evidence.**
- Gather with counts (P2C.2): exit 0, 5.47 s, 929,497,088 bytes peak RSS; the 634 files of the rnaseq copy hash the same before and after; 210 samples in the sheet and the counts, every one joined; 78,470 genes as Int32 from 0 to 20,313,260, none missing; `[expect]` met; nf-core/rnaseq 3.27.0, commit a1fcddd.
- `ctx de` 131.2 s at 1,551,384,576 bytes; `ctx function` 16.3 s at 836,239,360 bytes; `ctx report` 0.63 s at 121,356,288 bytes. Pages: M1 4,963,184 bytes, M5 12,773,948, M6 5,365,560, hub 67,666, all under 15 MB.
- `ctx function` first refused the ndn_hsct_early running-score figure: a 96-character caption line ended 0.07 in past a double page, because `style.caption_chars` budgeted the full width while `place_caption` draws at a 0.02 inset. The fix (d3d4f39) budgets the width inside the inset; layout only, and `style.py` is not frozen.
- The record, written by `scripts/python/record_sanity.py` and committed as ac4a0b5: `2026-09-30_de_tier1_sanity_v1.csv`, 15 rows, SHA-256 `ebb82fd7aa98aca35581cf77304b7a90ea026a24945bd024d423916c4a59812e`; `2026-09-30_de_tier1_metrics_v1.csv`, 92 rows, SHA-256 `641121260c9a0663a77d7228743800a5a248b2173f8749b99045c538d9a90453`.
- `tests/test_sanity_record.py` and `tests/test_record_sanity.py` on the branch: 22 passed; the freeze test skips once a record exists, by design, and `test_rule_unchanged` holds the freeze from there.

**Consequence.** The frozen paths are free: P2S.4 and P2P.3 can move M5, M6 and the hub to content files and the prose rules, and `annotate.py`'s legend layout can be fixed. The hub counts 11 modules (D-044).

---

## D-045 — Report prose: nine fixed labels, word budgets and a lint (draft v3.11)

**Date:** 2026-09-30 · **Status:** decided (Rytis)

**Context.** Rytis's draft v3.11 (`TEMP/todo9.md`, `TEMP/CLAUDE9.md`, `TEMP/lessons9.md`) corrects the pages' prose: narrative headings ("The question", "The idea") and figurative phrasing ("a finished download can still be incomplete") came from the plan's own skeleton wording. Its one change from v3.10 is the prose rules. The rest of `CLAUDE9.md` predates decisions made here and does not replace them: Arial and 12/10/9 pt (D-024), not Helvetica at 14 pt; frozen dataclasses, not pydantic; Liberation Sans, not TeX Gyre Heros.

**Options.** v3.11 puts prose in `src/ctx/report/prose/<module>.yaml` and runs `ctx lint-prose`. The repo already keeps page text in `content/<page>.toml` (D-024), read through `ctx.report.text` with tests for unused keys and unfilled placeholders; a `ctx` subcommand must be a module (`tests/test_cli.py::test_subcommands_in_modules`).

**Decision.** Rytis chose to keep `content/*.toml`, and to rewrite the frozen M5, M6 and hub pages right after the Step-1 record (P2S.4), not register again for them.
- Nine sections, labels verbatim: Purpose · Input · Method · Parameters · Results · Interpretation · Questions · Pitfalls · Reproduce. A page file names them `purpose`, `input`, `method`, `interpretation`, `questions`, `pitfalls`. `builder.FIELDS` maps them onto the `ModulePage` fields the frozen pages still pass, until P2S.4 renames the fields; `input` is optional until then, and a page without one shows no Input section.
- `ctx.report.prose` holds the rules: Purpose 30 words, Input 40, Method 90, each figure caption 40, Interpretation 120, Questions and Pitfalls 3 items of 20, a card's meaning and effect 12 each; sentences at most 25 words, paragraphs at most 3 sentences; a banned list; a quantity only from a `{placeholder}` or in `code` (a digit inside a name, Step-1 or SHA-256, is not a quantity). `pixi run lint-prose` runs it, and a test runs it in CI.

**Evidence.** `tests/test_prose.py`: each rule fails its planted violation, and names with digits, placeholders and code pass. `tests/test_report.py::test_skeleton_order` reads the labels from CONVENTIONS and compares them, in order, with the builder's. On the day it landed the lint listed 63 violations in the setup, gather and QC pages and the cards; the test that every page file passes lands with their rewrite.

**Consequence.** Every new page, the previews included, is written to the rules from the start. M5, M6 and the hub keep their text, under the new headings, until P2S.4.

---

## D-044 — The hub counts the course modules M0 to M10; registered again

**Date:** 2026-09-30 · **Status:** decided (Rytis asked for every module in the preview)

**Context.** The hub's "Modules run" came from `naming.MODULES` minus the report: M0 to M8 and the JupyterLite build, 10 in all. Gene cards (M9) and clinical (M10), which TODO.md files under `outputs/11_genecard` and `outputs/12_clinical`, were missing, and the JupyterLite build, which writes a site and never a module page, was counted. The preview hub read "4 of 10"; Rytis asked for every module. `hub.py` is frozen, so the fix registers again.

**Decision.** `MODULES` gains `genecard` (11) and `clinical` (12), the folder numbers TODO.md already gives them. The hub skips `NOT_MODULES`, the report and the JupyterLite build, and titles M9 and M10. It now counts 11 modules, M0 to M10, in course order. `tests/test_report.py::test_hub_lists_the_course_modules_in_order` asserts the order and each title's number, not only the count.

**Evidence.** The hub, naming and CLI tests pass: 50, plus the order test. The hub's verdicts, the only part of it the Step-1 record copies, read nothing that changed. At registration (2026-09-30 11:04) the Mac held no Step-1 count matrix; nf-core/rnaseq job 47485533 was still running.

**Consequence.** A Step-1 hub reads "N of 11". Pages for M9 and M10 get names the lint accepts (`*_genecard_*`, `*_clinical_*`), which the preview needs before those modules exist.

---

## D-043 — The read-length batch is named by its read lengths; registered again

**Date:** 2026-09-30 · **Status:** decided (Rytis)

**Context.** Step-1's batch, from read length, was named `short` (75 and 76 nt) and `long` (101 nt). Rytis: "long" reads as long-read sequencing, whose reads run to thousands of nt. He asked for "Reads ≤ 75 nt (min 35 nt)" and "Reads > 75 nt (max 500 nt)". `conf/datasets/tier1.toml` is frozen, so any rename registers again, and only while no Step-1 count is on the Mac.

**Options.**
- **The boundary at 75 nt, as first written.** Step-1 has 13 runs at 75 nt, 98 at 76 nt and 99 at 101 nt, all on one Illumina NovaSeq 6000, so read length is the only batch signal. The registered selection code on the gathered sheet gives, per contrast (reference/tested):

  | Contrast | Registered | Lower group | Upper group |
  |---|---|---|---|
  | ndn_gcsf | 13/13 | 1/1 | 16/16 |
  | mono_gcsf | 8/8 | 1/1 | 9/10 |
  | ndn_hsct_early | 13/5 | 1/1 | 16/6 |
  | mono_hsct_early | 8/4 | 1/1 | 9/5 |
  | ndn_gcsf_paired | 9/9 | 1/1 | 10/10 |
  | ndn_gcsf_unpaired9 | 9/9 | 1/1 | 8/8 |
  | ndn_pdac_within | 4/14 | none | 16/14 |
  | ndn_null_split | 7/6 | 1 sample | 8/8 |

  Every within-batch fit would mix 76 and 101 nt, the PDAC batch lesson and rule R6 would lose their meaning, and every contrast would need a new registered n.
- **The boundary at 76 nt,** Rytis's names with 76 in place of 75: the registered split, renamed.

**Decision.** Rytis chose 76 nt once the table was in front of him. The levels are "Reads ≤ 76 nt (min 35 nt)" (75 and 76 nt) and "Reads > 76 nt (max 500 nt)" (101 nt); min and max are words in the name, not enforced limits, so `dataset.py` and `curation.py` are unchanged. `tier1.toml` changes its bins, the 8 contrast filters that name a level, the sources of R6 and the recorded PDAC rule, and `registered`; `tier1_teach.toml` takes the same bins (`test_subsets_curate_as_their_parent`). The synthetic test data keep `short` and `long`: no student reads them. The preview renames them when it is built.

**Evidence.**
- Gathered again with the renamed file: 111 runs in the lower group and 99 in the upper, as before.
- All 10 contrasts select the same samples in the same groups as the committed file did on the 2026-09-29 sheet, compared sample by sample.
- pydeseq2 accepts the names in a `~batch + group` design: `ctx de` on the synthetic data with its batch renamed exited 0 in 38.1 s and wrote the adjusted contrast's results.
- **What was on the Mac at registration:** no count matrix. A copy taken while nf-core/rnaseq job 47485533 was still quantifying brought the 210 Trim Galore reports and, for 208 samples, Salmon's `lib_format_counts.json` and `meta_info.json` (library type and mapping summaries; all 208 report the library type `U`). No sanity rule reads any of them.

**Consequence.** Pages, tables and legends say the read lengths instead of short and long. Anything written before this commit that filters on `short` or `long` would find no samples; `ctx de` stops on the registered n, so it cannot pass silently.

---

## D-042 — The Track A subset is its own dataset file, a declared subset of tier1, curated by tier1's rules

**Date:** 2026-09-30 · **Status:** decided

**Context.** Track A gathers the 8 teaching runs of P3.0 from their own fetchngs outdir. `tier1.toml` cannot curate it: its `[expect]` demands 210 runs. A second file repeating tier1's tables was the M3 plan (R4), and its premise failed when measured: curation checks the SDRF join both ways, so each of the 202 SDRF rows with no run in the outdir is a problem, and `curation.py` and `dataset.py` are frozen until the Step-1 record.

**Options.**
- **Change curation and register again** before step 7: an eleventh registration, and a change to a frozen file for a dataset the Step-1 record never reads.
- **Wait for the record**, then change curation: Track A waits on P2C.3.
- **Cut the SDRF in the caller.** `harmonise.py`, not frozen, reads the SDRF text and hands it to `curate`. A dataset declaring `[dataset] subset_of` gets only the rows whose join key names a run of its outdir. `dataset.py` keeps `[dataset]` whole in `info`, so the key needs no loader change.

**Decision.** The cut in the caller. `conf/datasets/tier1_teach.toml` declares `subset_of = "tier1"` and repeats tier1's `[guard]`, `[sdrf]`, `[columns]`, `[order]` and `[derived]`, its `[dataset]` provenance, and the SDRF and GTF resources; `[expect]` pins 8 runs and their read lengths, conditions, cell types, tissue and donors. `tests/test_repo.py::test_subsets_curate_as_their_parent` holds every declared subset equal to its parent in those tables. Everything else curation checks still applies to the kept rows: a run with no row, cross-checks, value maps, `[expect]`. The README of the sheet states the rows kept. `tier1.toml` is frozen and has no `subset_of`, so Step-1's gather cannot take the new path.

**Evidence.** Measured 2026-09-30 on a scratch copy of the metadata-only Step-1 outdir cut to the 8 runs, as fetchngs writes 8 run inputs (one `metadata/<run>.runinfo_ftp.tsv` each):
- the file without `subset_of`: exit 1, 202 problems, every one "SDRF row matches no run_accession";
- the same with that one check dropped in memory: exit 0, `[expect]` met; a donor count of 1 for UPN5 exits 1 with "donor = UPN5: 2 runs; expected 1";
- the committed file, no patch: exit 0, 4.05 s, 432,013,312 bytes peak RSS, README "8 of 210, the runs of this subset of tier1".

The equality test fails each of six planted edits to the subset file (a value map, a bin, a guard pattern, a resource hash, an extra table, the source line). The CLI tests fail both planted edits to the cut (keeping every row; reading another key).

**Consequence.** Track A gathers as soon as P3.0's outdir reaches the Mac, with no change to the registration. A second cohort's teaching subset is one more file with `subset_of`.

---

## D-040 — Agreement beside correlation: Spearman with a percentile bootstrap, Pearson with Fisher's z, Bland-Altman on log2(count + 1)

**Date:** 2026-10-01 · **Status:** decided; the caption half (M3 R6.1, approved by Rytis 2026-10-01) is recorded with the figures in M3.6

**Context.** M3 compares two counts of one gene: STAR's against featureCounts', and one featureCounts option against another. M7 will compare RNA with protein. A correlation near 1 says two methods order the genes alike, not that they give the same numbers (Bland and Altman 1986). No helper existed: `ctx.viz.annotate` is frozen and tests groups with Mann-Whitney only.

**Options.**
- **Correlation alone**, as TODO.md's figure table asks: Spearman with a bootstrap interval.
- **Correlation beside Bland-Altman's bias and limits of agreement** on log2(count + 1).
- **Lin's concordance correlation coefficient**, one number for both; kept as a think prompt.
- For Spearman's interval: BCa, percentile, or a closed form (Bonett-Wright, Fieller), the closed forms read only in secondary sources.

**Decision.** `ctx.concordance.agreement(a, b)`, one definition for every module:
- genes with a + b > 0, filtered before any transform; a gene zero in one method only is kept, because it is a disagreement;
- the size rules come first, for identical counts as for any: no correlation from fewer than two genes or a constant vector, no SD from fewer than two genes, no interval from fewer than four;
- both correlations come from one formula, Pearson's r as the sum of products of centred values over the root of the product of their sums of squares, clipped to [-1, 1] as `pearsonr` clips it: on average ranks it is Spearman's rho, for the estimate and every resample alike, and on log2(count + 1) Pearson's r. It is exactly 1 on identical counts, and rho exactly 1 on counts with identical ranks, so no short-circuit is needed;
- Spearman's rho as `scipy.stats.spearmanr` defines it; its interval from the draws of `scipy.stats.bootstrap` with `paired=True`, 1,000 resamples of genes, seed 42: the 2.5th and 97.5th percentiles of the finite draws. The statistic is the vectorised rank correlation above. A resample of tied values has no rank correlation: it is dropped and counted, where scipy's own interval turns one such resample into NaN. scipy's own interval is not read, so its `method="percentile"` only spares scipy the BCa jackknife;
- Pearson's r on log2(count + 1), with Fisher's z interval tanh(atanh(r) ± 1.959964 / √(n - 3)), the formula of `PearsonRResult.confidence_interval`, and the point itself when |r| = 1;
- d = log2(b + 1) - log2(a + 1); bias = mean(d) with a t interval; limits bias ± 1.96 SD, the SD with n - 1; the empirical 2.5th and 97.5th percentiles of d beside them; fold = 2^bias, a ratio of (count + 1) values, never called a log ratio (Euser et al. 2008);
- `per_unit` gives one row per sample and an `all` row on each gene's sums over samples, so the gene stays the resampled unit.

Percentile, not BCa, because BCa's jackknife returned NaN on degenerate input (M3 review 1), not for speed. No hypothesis is tested, so there is no BH. The intervals are nominal: the genes are the chromosome's whole annotated set, not a sample of genes, so a page leads with the spread across samples. Intervals for the limits of agreement wait until Bland and Altman 1999 is read in full.

**Evidence.** Measured 2026-10-01 on this Mac, scipy 1.18.1:
- coverage: 1,000 bivariate normal datasets, n = 60, Pearson 0.6, whose population Spearman is (6/π)·asin(0.3) = 0.5819. The 95 % interval covered it 960 times; a [5, 95] interval 913 times;
- identical counts: `spearmanr` of a 17-gene vector with itself gave 0.9999999999999998. Over 2,000 non-constant random count vectors of 2 to 399 genes, `spearmanr` of a vector with itself missed 1 on 500 and `pearsonr` on log2(count + 1) on 784; the one formula missed on none, and `agreement` gave r 1 on 1,999 and no r on the one left with a single gene;
- corrected after the M3 review, 2026-10-01: the first build answered identical counts with a short-circuit before the size rules, so one identical gene got rho 1 and the interval [1, 1], and it took the estimate from `spearmanr`, so counts one read apart with identical ranks gave rho 0.9999999999999998 below its own interval [1, 1]. Both came from the short-circuit and the second formula, which are gone;
- time at 1,747 genes, chr22's count in GENCODE 50 (synthetic Poisson counts, 1,682 kept), measured again after the review fix: `agreement()` took 0.155 to 0.167 s in five runs, and 0.150 to 0.154 s on identical counts, which now run the bootstrap; the whole process 2.66 s and 290,390,016 bytes peak RSS (`/usr/bin/time -l`), most of it the imports. Before the fix: 0.164 to 0.167 s;
- `tests/test_concordance.py`, 16 tests: 25 of 26 planted mutations fail their test. They are Pearson on the raw counts as the estimate, ordinal ranks, 1/√n, z = 1.96, a [5, 95] interval, a basic interval, resampling gene × sample points, a nonzero-in-both filter, genes zero in both kept, ddof 0, d as a - b, an interval from three genes, scipy's NaN interval, no check on negative counts, the estimate from `spearmanr`, r from `pearsonr`, r not clipped, the old short-circuit before the size rules and from four genes on, seed 43, 200 resamples, the seed argument ignored for 7 or for `SEED`, no guard on zero genes, and `per_unit` without its duplicate check. The 26th, `method="basic"` passed to scipy, changes no output, since scipy's own interval is not read. The first record here said 12; that sweep held 13.

**Consequence.** M3's concordance table and agreement figures read these fields, and M7 imports the same module. A change of interval method changes every agreement caption at once.

---

## D-039 — The strand is called from featureCounts' forward share under its own registry entry; the registry gains an "align" stage

**Date:** 2026-10-01 · **Status:** decided; the registry half landed in M3.2. The per-sample call, the stop on an undetermined call and the comparison with Salmon land in M3.5 and M3.6.

**Context.** M3 counts at the library's strand, and the sheet says "auto" (`harmonise`). On Track A, M3 runs before M4, so the subset has no Salmon output of its own; Salmon's call reaches M3 only through `--compare-dataset`. The registry knew four stages, and `test_qc_registry.py` requires every stage to have a metric. `salmon_strand_bias.computed` said SF / (SF + SR), which holds for single-end runs only.

**Options.**
- **Reuse `salmon_strand_bias`** for featureCounts' share: one entry for two measurements.
- **STAR's ReadsPerGene columns 3 and 4**: no call under HISAT2.
- **featureCounts' forward share** (Assigned at `-s 1` over `-s 1` plus `-s 2`) under its own entry with nf-core's thresholds, and Salmon's call shown beside it when `--compare-dataset` brings one.

**Decision.** The third. `conf/qc_metrics.toml` gains stage "align" (`registry.STAGES`) with two metrics:
- `gene_strand_bias`: kind strandedness, stranded 0.8, unstranded 0.1, on process's `strandedness` table, column `forward_share`. The thresholds are nf-core/rnaseq 3.27.0's, set on Salmon's fragments over transcripts, so applying them to gene-assigned reads is a course convention with that rationale;
- `reference_start_codons`: numeric, higher is better, good 0.99, usable 0.9, action exclude, on process's `reference` table, column `atg_share`. A genome and annotation from different assemblies make every later number meaningless.

`salmon_strand_bias.computed` now says SF / (SF + SR) for single-end runs and ISF / (ISF + ISR) for paired-end ones. `Metric.strandedness` serves both strand entries. Neither new metric is in MultiQC's General Statistics, so the generated MultiQC config is unchanged, and `ctx qc`'s judge skips both, since its tables carry neither column.

**Evidence.** Measured 2026-10-01:
- the 5 `lib_format_counts.json` of `tests/fixtures/rnaseq`: `strand_mapping_bias` equals SF / (SF + SR) for the 2 single-end reports (RAP1_UNINDUCED_REP1 0.1300) and ISF / (ISF + ISR) for the 3 paired-end ones (WT_REP1 0.1037), where SF / (SF + SR) gives 0.5016;
- start codons, by hand on GENCODE 50 and the GRCh38 chr22 sequence of D-035: 8,852 of 8,866 read ATG (99.84 %; 7 CTG, 6 GTG, 1 ACG; 12 split records skipped), and 0 of 8,866 with the sequence shifted by one base;
- `tests/test_qc_registry.py` fails each planted mutation: gene_strand_bias's unstranded set to 0.01, its stranded set to 0.85 or 0.75, no "align" stage, the single-end-only Salmon text, the start-codon good band at 0.995 and usable band at 0.96.

**Consequence.** M3.5 calls each sample from `gene_strand_bias`; M3.6 refuses a reference under `reference_start_codons`' usable band; M3.12 adds M2b's dataset-only align bands to the same stage.

---

## D-037 — M3 builds its own chr22 indices into content-keyed tool folders under `outputs/03_process/tools/`

**Date:** 2026-10-01 · **Status:** decided; the lint rule, the keyed folders and the index builds landed in M3.2 and M3.3. `--genome`, the ATG refusal and `tools/README.md` land with the command in M3.6.

**Context.** STAR, HISAT2, fastp and featureCounts write files under their own names (`SAindex`, `genome.1.ht2`, `ReadsPerGene.out.tab`), which the naming lint fails, and MultiQC finds a tool's logs only by those names. Where the indices are built was contradictory: `ctx setup` builds indices into `resources/` (TODO.md M0), while S0.4 and P3.0 say M3 builds the chr22 indices on the laptop. An index, a trimmed run or an alignment is costly to redo and cheap to reuse, and a rerun with one option changed must not overwrite what an earlier page describes.

**Options.**
- **Renamed copies** of every tool file under the naming rule: MultiQC no longer finds them, and every tool output is written twice.
- **`outputs/03_process/tools/` exempt at any depth**, as `multiqc` was: a `tools` folder anywhere in any module escapes the lint (M3 reviews 1 and 2).
- **Exempt only directly under the module folder**, with content-keyed subfolders.
- **Indices in `resources/` by `ctx setup`**: setup fetches published files with checksums, and a chr22 index has no URL.

**Decision.**
- `naming.NATIVE_DIRS` is {multiqc, tools}, matched only as the first folder under the module folder; a `tools` or `multiqc` folder deeper down is linted.
- `ctx.process.cache.keyed` names a folder `<prefix>_<first 12 of the SHA-256 of its inputs as canonical JSON>`, builds it into `.part`, writes `ctx_manifest.json` (key, inputs, each file's size) and renames it only after the build returns. It is reused only when the manifest names the key and every listed file has its listed size. A damaged folder (a listed file short or gone) or a leftover `.part` is renamed to `.stale<N>`, never deleted and never reused. The folder is resolved to an absolute path, because each builder runs its tool inside it: with a relative parent, STAR exited 104 and hisat2-build 1 (M3 review). One helper serves the reference, the indices, and later the trimmed runs, alignments and counts.
- M3 builds the indices (`ctx.process.reference`): `prepare` writes the genome uncompressed (STAR refuses gzip) and the GTF lines on its contigs, matched exactly, into `tools/index/reference_<key>/`; with `half`, each contig's first floor(L / 2) bases and the genes whose gene record ends inside them. `star_index` and `hisat2_index` key on the tool, its pinned version, the copied genome's and subset GTF's SHA-256, `half`, and for STAR `--genomeSAindexNbases` = min(14, floor(log2(L) / 2 - 1)) and `--sjdbOverhang` 100, STAR's default and the ideal for 101 nt reads. HISAT2's index is the plain one. `hisat2-build` joins `tools.TOOLS` at the hisat2 pin.
- The genome is given by path and recorded by the SHA-256 of the file as given, beside the start-codon check (D-039); the Salmon chr22 indices move to M4.

**Evidence.** Through the code on this Mac, 2026-10-01, 4 threads, `/usr/bin/time -l`, into a scratch folder outside the repo; genome `GRCh38.primary_assembly.chr22.fa.gz` (SHA-256 32f26397..., D-035), annotation GENCODE 50 (SHA-256 89bbad69..., as `tier1.toml` declares):

| Step | Wall | Peak RSS, bytes | On disk, bytes | What it found |
|---|---|---|---|---|
| reference | 6.45 s | 200,261,632 | 51,665,450 genome | 50,818,468 bases, 11,658,691 N; 1,747 genes, 278,849 lines; 8,852 of 8,866 start codons ATG |
| STAR index, SAindexNbases 11 | 13.04 s | 757,301,248 (STAR) | 434,837,566 | |
| HISAT2 index | 5.87 s | 300,236,800 (hisat2-build) | 64,166,924 | |
| half reference | 6.21 s | 211,566,592 | 25,832,729 genome | 25,409,234 bases; 713 genes; 2,242 of 2,243 ATG |
| half STAR index | 4.60 s | 539,181,056 | 174,690,876 | |
| half HISAT2 index | 2.72 s | 174,981,120 | 25,392,050 | |

A second run reused all three folders in 0.05 to 0.06 s each. On the fixture slice, the STAR index the Mac binary built into `.part` and renamed then mapped `se_a` under linux STAR 2.7.11b with GeneCounts (the D-036 image, `--pull never`): 1,000 input reads, 99.90 % unique; the index records the `.part` path in `genomeParameters.txt`, and mapping succeeded with that path gone. After the review fix, the same three builds from a relative parent folder: reference 6.54 s and 199,442,432 bytes, STAR 12.69 s and 756,793,344 bytes, HISAT2 5.73 s and 293,470,208 bytes, the same sizes on disk and 8,852 of 8,866 ATG. `tests/test_process_cache.py` and `tests/test_process_reference.py` fail each of 21 planted mutations, among them a key without one input, reuse without the size check, a missing file read without checking it is there, a leftover `.part` finished by renaming, `round` or `ceil` in the SA rule, prefix matching of contigs, a gene kept by its start, `half` left out of the reference's key, the annotation out of STAR's, the genome out of HISAT2's, and the folder left relative. The first record here said 14; that sweep held 16.

**Consequence.** `outputs/03_process/tools/` holds content-keyed folders that datasets and reruns share, and `lint(outputs/03_process)` skips it whole. A change of key fields rebuilds beside the old folder; deleting the folders is the reader's call.

---

## D-036 — STAR runs on Linux, HISAT2 on Apple silicon; a probe at the door, never a silent fallback

**Date:** 2026-09-30 · **Status:** decided (Rytis); the probe and its tests land in M3.4

**Context.** M3 aligns with STAR, pinned at 2.7.11b in the `reads` environment (D-027), whose check is `STAR --version`. The osx-arm64 build passes that check and does not align. Pixi on Apple silicon is a student platform, so this decides what a Mac student runs, not only what this Mac runs.

**Options.**
- **A container** for STAR only, behind a per-tool engine switch: Docker's VM joins the Track A memory budget, and the image is either a third-party one outside `pixi.lock` or a project image that needs linux-aarch64 added to the `reads` feature.
- **HISAT2 on Apple silicon, STAR on Linux** (CI, the cluster, Linux students).
- **A STAR build carrying the unmerged fix** (PR #2691 for issue #2632): off the pin, and not shown to cure the second fault.

**Decision.** HISAT2 on Apple silicon, STAR on Linux; `--aligner star` stays the default. Before any alignment, `ctx process` maps 10 exact reads through STAR with the production options (`--readFilesCommand gzip -dcf --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within`) and requires all 10 read and every ReadsPerGene column to sum to 10. A failing probe stops the run, names this entry and both routes (`--aligner hisat2`, or a working STAR first on `CTX_TOOLS`), and writes no page; nothing switches aligner on its own. A test needing STAR may skip only when the probe fails on darwin-arm64 with `CI` unset, and fails everywhere else, a departure from D-027 that waits for Rytis's confirmation (M3 R6).

**Evidence.** Reproduced here, 2026-09-30, `.pixi/envs/reads/bin/STAR` 2.7.11b on osx-arm64, the fixture slice index built by the same binary (`--genomeSAindexNbases 7`, `slice.gtf`), `tests/fixtures/chr22_se/se_a.fq`, 1,000 uncompressed reads:

| Run | Exit | What STAR reports |
|---|---|---|
| `--outSAMtype BAM Unsorted` | 0 | `Number of input reads` 0, uniquely mapped 0.00 % |
| the same with `--quantMode GeneCounts` | 109 | `could not open input file /geneInfo.tab`, although `geneInfo.tab` is in the index folder |

The input is plain FASTQ, so the first fault is not about gzip. The second path has lost its genome-directory prefix; why is not traced. The same reads under linux STAR 2.7.11b (a linux/arm64 image already on this Mac, run with `--pull never`) against a slice index it built itself, with GeneCounts: 1,000 input reads, 99.90 % uniquely mapped. Index building works on the Mac: the real chr22 index built by the Mac binary (12.07 s, 757 MB peak) aligned a proxy input under that linux STAR. An exit code and `--version` both pass on the failing build, so neither can serve as the check.

**Consequence.** The Mac's Track A page is a HISAT2 page (M3.9), without the STAR against featureCounts comparison, which needs STAR's own counts; the STAR page comes from Linux. If a fixed osx-arm64 STAR is pinned later, the probe passes and nothing else changes.

---

## D-035 — The Step-1 download relaunched fresh, ten at a time, retrying refused connections; chr22 cut from the cluster's genome

**Date:** 2026-09-30 · **Status:** decided

**Context.** P2A.3 step 1, fetchngs 1.12.0 on all 210 Step-1 runs through the lab launcher, ended "Pipeline completed with errors". `SRA_FASTQ_FTP` got `wget ... failed: Connection refused` (exit 4) again and again, and fetchngs killed 48 running transfers. Rytis pasted the state:
- 2 FASTQs and 2 md5 files published (2.0 GB), both hard-linked;
- no samplesheet;
- 47 GB in `work/`, mostly killed partial transfers.

The compute nodes did download data, so network access itself works. The refusals are consistent with ENA turning connections away under the launch's 50 concurrent downloads (`queueSize = 50`), but ENA's limit was not looked up and the cause is not proven.

**What made it fatal.**
- fetchngs passes `wget -t 5 -nv -c -T 60` (`modules/local/sra_fastq_ftp/nextflow.config`).
- GNU wget treats a refused connection as final unless given `--retry-connrefused`, so the five tries never happened.
- The process's `error_retry` label allows 2 retries (`conf/base.config`), then fetchngs ends the whole run.

**Options.**
- **`-resume` into the same outdir:** reuses 2 GB, and writes a second `params_*.json`, which `ctx gather` refuses by design, since two params files in one outdir normally mean two runs.
- **Adding `--retry-connrefused` to `ext.args` and resuming:** changes every download task's script and hash, so the resume re-downloads everything anyway.
- **A fresh run:** move the failed `nfcore_fetchngs/` aside, and relaunch with the download settings changed.

**Decision.** A fresh run. Moving the failed directory aside keeps it (never deleted) and gives the relaunch a new launch directory, so the launcher's `-resume` finds no cache and one params file is written. The launch-directory `nextflow.config` adds, for `SRA_FASTQ_FTP` only:
- `ext.args = '-t 5 -nv -c -T 60 --retry-connrefused --waitretry 60'`;
- `maxForks = 10`, my choice, not an ENA figure;
- `maxRetries = 6`, with exponential backoff: 1, 2, 4 and so on up to 32 minutes between attempts.

The P3.0 block takes the same settings. `ctx gather`'s one-params-file rule is unchanged. Rytis asked, at 08:35, for a 30-minute limit per download: `time = { 30.min * task.attempt }` replaces fetchngs's `process_low` 4 h (checked with `nextflow config`: `withName` wins over `withLabel`), and grows with each retry so a slow transfer is not killed for good. It applies to later launches: the running job 47484816 keeps 4 h, because its head job never rereads the config and a restart would resume into the same outdir. At about 08:45, 155 of 210 FASTQs were published; five downloads had failed once with wget exit 8, and ERR7425285's error was `ERROR 403: Forbidden` at 07:50:39. The same five URLs answered `200 OK` to a HEAD request from the Mac a few minutes later, and ENA's filereport lists ERR7425285 at the path fetchngs used, 4,152,628,581 bytes (the HEAD's Content-Length too), so the files exist and the 403 looks like a temporary block of the cluster's address; each gets up to 6 retries.

For M3, step 2b cuts chr22 out of the cluster's MD5-checked GENCODE 50 genome, so the Mac's STAR index aligns to the sequence the Step-1 counts were quantified against, with no new download.

**Evidence.** `nextflow config nf-core/fetchngs -r 1.12.0 -profile singularity`, run here with `NXF_SYNTAX_PARSER=v1` against this config, shows it merged over the pipeline's: the new `ext.args`, `maxForks 10`, `maxRetries 6`, the backoff closure, and the module's own `publishDir` kept. Rytis's paste, 2026-09-30 07:35: the cluster's wget 1.20.1 image lists both `--retry-connrefused` and `--waitretry`; the relaunch is job 47484816, and exactly 10 `SRA_FASTQ_FTP` jobs ran at once, so `maxForks` holds. Step 2b: one header, `>chr22 22`; 50,818,468 bases; 11,389,750 bytes gzipped; SHA-256 `32f2639716283f0d8599f540a9a6105d9e6577ea4599e9ea5c0fbc2365872771`, the same on the cluster and on the Mac after rsync, at `~/ctx_tier1/data/raw/reference/GRCh38.primary_assembly.chr22.fa.gz`.

**Consequence.** The relaunch re-downloads the 2 finished FASTQs. The failed run's 47 GB `work/` can be removed once the new run is checked; that is Rytis's call.

---

## D-034 — `ctx qc` runs three stages into one judgement: counts, the pipeline's reports, and the reads

**Date:** 2026-09-30 · **Status:** decided

**Context.** M2.6 turns the checkpoint (D-030), the read QC (D-031) and the figure contracts (D-032, D-033) into one command and one page. The command had to decide where the registry comes from, how the three stages' values meet, and how MultiQC's reports are laid out.

**Decision.**
- **The registry** is resolved from `--registry`, else `$CTX_QC_REGISTRY`, else the repository's `conf/qc_metrics.toml`. When none exists the run stops and says which to set; the page names the file used.
- **Counts** are always measured. The composition shares need the dataset file's `annotation_gtf` resource on disk; without it they are left out, and the page says which case it was.
- **`--rnaseq-dir`** runs MultiQC on the outdir's Trim Galore and Salmon reports, into `outputs/02_qc/multiqc/`, coloured by the generated config written beside the tables:
  - a paired sample takes read 1;
  - Trim Galore's fates must add up to its total, and Salmon must process exactly the reads written, or the run stops naming the samples;
  - each sample's reports are checked by MultiQC's table rows, because MultiQC 1.35's Salmon module lists no sources.
- **`--reads`** runs FastQC, fastp and seqkit on the sheet's FASTQ files:
  - MultiQC reads them into their own report, `multiqc/reads/`, because FastQC names samples by file, and one table would mix two naming rules that `by_unit` rightly refuses;
  - the reads are judged run by run, each run in its sample's stratum, and status gives each run its sample's flags and its own;
  - for the confounding tests, a sample takes its first run's read values and flags, the run M5's selection keeps;
  - the fastp grid runs on `--grid-run`, else the first run read, over `-q` 10 to 35 and `-l` 15 to 102. The largest `-l` passes the longest Step-1 read, so the map shows where the length floor bites. A grid where every cell keeps the same share is said in words, not drawn.
- **FastQC's plot bands**, good from Phred 28 and reasonable from 20, are the registry metric's `plot_bands`, not figure code. They were read from the bytecode of FastQC 0.12.1's `QualityBoxPlot` with `javap` (constants 28.0 and 20.0 beside its GOOD, BAD and UGLY colours), because FastQC's help names the colours without the numbers.
- **The fastp settings** live in a dependency-free `ctx.qc.Fastp`, so the CLI reads their defaults without importing polars before the runtime profile sets its thread limits.

**Evidence.** `tests/test_qc_run.py` runs `ctx` as a subprocess on synthetic gather outputs:
- a shallow library, a swapped sex label and an MT-rRNA-heavy sample are each flagged;
- `--rnaseq-dir` on the rnaseq fixture writes the fates, the bubble chart, the Sankey and MultiQC's report, and leaves the outdir untouched;
- `--reads` on the chr22_se samples flags se_c's adapter and se_d's per-base quality run by run, and draws the curves and a 36-cell grid that varies.

Two defects surfaced and were fixed at their source:
- `multiqc_io.run` ran MultiQC in its output folder, so a relative input pointed nowhere. A test with a relative path fails without the fix.
- An empty flags table had untyped columns and could not be joined. A test fails without the typed schema.

**Consequence.** M2.6d puts every figure on the page. On Step-1, `ctx qc` runs with `--rnaseq-dir` once the counts exist, and with `--reads` on the Track A subset.

---

## D-033 — Linked brushing: a scatter over parallel sets, painted by `ctxBrush`; page controls get a `kind`; a Copy button that printed a method

**Date:** 2026-09-30 · **Status:** decided

**Context.** P2Q.4 asks for linked brushing on the QC page: a lasso on a scatter, or a tap on a category of the parallel sets, paints those samples in both views. It also asks for a colour picker, Clear, a live summary, and a vector parallel-sets diagram for the PDF and SVG. The page had one kind of control, the volcano's live calls, which it found by `form.calls`.

**Decision.**
- **`LinkedData`** (`ctx.viz.linked`) is a scatter over a `go.Parcats` of the dimensions, in one plotly figure. Its static half is the parallel sets drawn through `FlowData`: a node per level, samples as reads, and the bands into the colour dimension (QC status) in that dimension's colours.
- **Page controls carry a `kind`.** `"calls"` is the volcano; `"brush"` renders a radio picker of the seven data colours, named in Okabe and Ito's words, a Clear button and an `aria-live` line.
- **`ctxBrush`** (`templates/_brush.js.j2`) holds the logic, and the page wiring only restyles:
  - one paint index per sample, and a repaint replaces the old one;
  - marker colours, with the reference grey for unpainted samples;
  - parcats codes, with 0 for unpainted;
  - the live line, per paint in palette order: its count and its counts by one dimension, most first.
- **After a paint** the page drops Plotly's own selection dimming and the lasso outline, so painted samples show at full strength.
- **One scatter view, not several.** A second view waits until a page needs one.
- **The Copy button.**
  - `ui.toml`'s `reproduce.copy` is now `reproduce.button`, and the new `brush.clear` is `brush.erase`. Jinja reads `ui.reproduce.copy` as the dict's own `copy` method before the key.
  - The Copy button on every module page built since `e662d1d` (S3a, 2026-09-29) read `<built-in method copy of dict object at 0x…>`. The browser checks after S3a missed it.
  - `tests/test_text.py::test_no_ui_key_shadows_a_dict_method` now fails on any `ui.toml` key named like a dict method.

**Evidence.**
- **Under node.** `tests/test_page_script.py::test_brush_paints_both_views_and_counts_as_python_does` runs `ctxBrush` against a Python oracle over four paints, including repaints and both event shapes, then Clear. Four mutations were each caught: codes off by one, the breakdown least first, no repaint, `pointIndex` read first.
- **In the browser**, at 1280 px:
  - A real tap on the LDN category sent `plotly_click` with 12 points of `{curveNumber: 1, pointNumber: <sample>}`, the shape `ctxBrush.indices` reads. It painted all 12 LDN samples in both views: "12 painted vermillion: LDN 12".
  - A `plotly_selected` event emitted through Plotly's own emitter painted "22 painted bluish green: NDN 21, LDN 1". Clear returned both traces to grey. The console showed no errors.
- **Not verified in the browser:** a real lasso or box drag. The pane's drags and synthetic mouse events drew no selection with points, so the scatter's own event shape rests on Plotly's documentation and the node test.
- **Page structure.** `tests/test_report.py::test_each_kind_of_control_gets_its_own_form`: the linked page has the brush form, seven paints and the empty line, and no calls form. The volcano has the calls form and no brush form.

**Consequence.** M2.6 puts a `LinkedData` of PCA scores and the five dimensions on the QC page. Every module page must be rebuilt for the Copy label.

**Amended 2026-10-01 (Rytis): the parallel sets across the card, and a PDF and SVG of the view as set.** Rytis found twice that nothing a reader changes on a page reached the PDF and SVG buttons: painting the linked view, and moving the thresholds of the M5 volcano (Step-1's Figure 4, "Volcano: gcsf over steady_state"). A bar's filters were the same. His choices: the parallel sets span the figure card, under the reading guide too; a moved view's buttons draw that view in the browser, at the journal sizes, with a small open-source PDF library embedded.
- **Layout.** The linked view is two plotly figures: the scatter beside the reading guide (`interactive()`), the parallel sets under both across the card (`below()`, `FigureBlock.below`, `.fig-below`). The brush addresses trace 0 of each; `ctxBrush` paints as before and gains `moved()`. The parallel sets get margins for their outer levels, measured in Liberation Sans advance widths (`ctx.viz.metrics`), since plotly sets those levels outside the plot. At 375 px the three stack; the parallel sets keep the least width their titles need (`layout.meta.min_width`) and scroll sideways under a line saying so, where a phone had squeezed four titles into 191 px.
- **One export path for every control kind.** Each kind's wiring tells the page script whether its view has moved from the defaults: `ctxBrush.moved`, `ctxCalls.moved` and `ctxFilter.moved`. A unit switch alone is no move. At the defaults the links open the matplotlib files, as before. Moved, they read "PDF of this view" and "SVG of this view" (`ui.toml` [export]) and link no file, so a middle click or "Save link as" cannot reach the default file under the view's label. `ctxExport` (`_export.js.j2`) draws the figure's plots as they stand:
  - marks back at the figure files' sizes, the page's own zoom divided out;
  - the journal template, `style.journal_template()`: the white theme at 9 pt ticks and legend, 10 pt axis titles, 3 pt ticks, 1 pt lines;
  - the figure files' size in pt from the contract's `inches()`; a painted linked view takes 170 by 153 mm (`ASPECT["linked"]`), the scatter over the parallel sets;
  - what the figure files draw and the page's plot lacks, from the contract's `export_layers()`: a thinned volcano's left-out points as float32 pairs in base64, its gene labels at the files' offsets, the parallel sets' outer margins in pt;
  - dashes as matplotlib draws them, `style.DASH_PATTERNS`, which matplotlib's rc reads too, where plotly drew "dash" 9 on 9 off; legend keys opaque, as in the files; colours written at four decimals, where jsPDF wrote two, up to 1/255 off;
  - the caption: the figure files' caption, never the page's thinning sentence; the live line; and `ui.toml` export.note, "Drawn in the browser from the view as set on the page, at the figure files' width.", each closed by a stop (the D-051 plan's question 1).
  The SVG is plotly's drawing sized in pt. The PDF is that SVG through svg2pdf.js 2.8.1 into jsPDF 4.2.1, with Liberation Sans 2.1.5 regular and bold embedded, each under its PostScript name. Back at the defaults, the links open the files again.
- **The libraries.** Vendored in `src/ctx/report/vendor/` with their licences and `MANIFEST.toml` (source, version, SHA-256, licence). jsPDF and svg2pdf.js came by `npm pack` from the official registry, each tarball's SHA-512 equal to npm's published integrity. The fonts were found on this Mac first, a copy inside the fastqc-rust 1.0.1 crate, and are byte for byte the official 2.1.5 release, whose tarball measured SHA-256 `7191c669…` on arrival; GitHub publishes no checksum for it. `ctx.report.export` refuses a file whose hash differs. A page inlines them once, and only when a figure on it has page controls, with the bundles' source-map line dropped. Each script follows its manifest's licence notices, its own and those of the packages its bundle holds: the svg2pdf.js bundle printed none, which MIT does not allow in a copy (review, 2026-10-01).
- **Who pays the bytes (the D-051 plan's question 2, Rytis 2026-10-01).** Whole, the two fonts were 1,100,224 of the export's 1,630,308 bytes on every page with controls. `export.subset_font` cuts them at build time with fontTools' subsetter, which ships with matplotlib. It keeps:
  - every character of the ctx package's text: `content/*.toml`, the glossary, and the code that still names axes and captions, `ctx.viz` among it;
  - every character the page's figures with controls carry (`figures.prepare`), their data's names among them;
  - a margin for a name or sign the package has not set yet: printable ASCII, Latin-1 and Latin Extended-A, Greek, general punctuation.

  That is 506 of each font's 2,327 characters. Outlines, hinting and all 30 name records stay, and no timestamp is renewed, so the same characters give the same bytes. The page prints the fonts' OFL text, its Reserved Font Name with it, once, beside the manifest's note that they are subset for this page's PDFs only. The OFL FAQ (1.1-update7) allows a subset embedded in a document (1.11, 1.12), and counts a subset webfont delivered to a browser as a Modified Version that may not keep a Reserved Font Name (2.6). A page carries these fonts only to embed them in the PDFs it draws of its own figures, so this treats it as embedding, an open question for Rytis. Per page with controls: fonts as base64 1,100,224 to 195,608 bytes; scripts with every notice 521,527 to 526,375; the export in all 1,630,308 to 730,540. In the same change a page draws a scatter's line through the points it needs (`scatter.needed_points`, Ramer-Douglas-Peucker): none of the line's own points lies over 0.25 px (`style.LINE_TOLERANCE_PX`) from the line drawn, at twice the widest page plot. The figure files draw every point. M5's dispersion trend carried every gene: 12,796 on Step-1's featured contrast. M6's running score thins by the same rule.
- **Where the file still differs from the figure files.** plotly edges each parallel-sets band in lightgray and its colour segments in black at 0.2 pt, which the export drops, as the files draw none. jsPDF opens each page at 0.200025 pt; no stroke draws at it. A painted linked view is 170 by 153 mm, its files 93.5 mm tall. Tick values are plotly's.

**Evidence (amendment).**
- **Under node**, `tests/test_page_script.py`: the caption wrap against textwrap; the caption's three parts, each closed by a stop; every mark at 6 pt and every line at 1 pt from a page at zoom 1 and 2; the SVG's size in pt, its caption inside the insets and escaped; the face set in the PDF's copy; `moved()` false at each default only; the thinned points, labels, dashes and opaque keys the export adds.
- **In headless Chromium** 152, `tests/test_export_browser.py` (29 tests, skipped without Chromium) on pages built from the demos: the buttons at the defaults, moved and reset; size, fonts, type sizes, strokes and colours of each file; the paint, the moved calls and lines, the filtered bars; a thinned volcano's every point; the labels; bold through the subset fonts; every text of every file on its page; no request and no error. On module pages `ctx` built from the test fixtures (M2 by `ctx qc`, M5 by `ctx de`, M1 by `ctx gather` on the curated set), pressed by `tests/browser.mjs` on 2026-10-02, before the rename; its Fonts column is superseded under Evidence (rename), below:

| Fixture page | Moved | PDF, mm | Fonts | Change in the files |
|---|---|---|---|---|
| M2, linked view | 6 samples orange, 3 vermillion | 170 by 153 | LiberationSans | 6 and 3 points in the paints, both paints in the bands |
| M5, volcano | log2 FC 2, BH p below 0.01 | 170 by 93.5 | LiberationSans | 48 up, 46 down, as the page; cutoffs read off the SVG's axes at 2.0 and -log10 p 3.0101, the page's line 3.0119; the six labelled genes of the files, G00008 to G00083 |
| M1, library sizes | cell_type type_a | 170 by 85 | LiberationSans | 8 bars of 14, named as the page |

  Each file holds 9 and 10 pt text only, strokes at 1 pt and none, no colour off `ctx.viz.style`, the live line and the note; dashes 3.7 on 1.6 off, and 1 on 1.65 for the M1 dotted line, as the files; each embedded font program 17,028 to 18,856 bytes. Each page made no request but itself and no console error.
- **Bytes**, fixture pages at the base commit, with the whole fonts, and now:

| Page | 424459f | Whole fonts | Now |
|---|---|---|---|
| M1 overview, a filter | 4,939,171 | 6,591,597 | 5,692,220 |
| M2 overview, the brush | 5,080,356 | 6,743,179 | 5,843,816 |
| M5 overview, live calls | 6,663,023 | 8,318,324 | 7,382,162 |
| M6 overview, no controls | 4,996,386 | 5,017,862 | 4,993,457 |

  On the fixture, M5's dispersion trend went from 997 points to 32 and M6's running score from 997 to 15. The Track B budget test at 60,000 genes by 171 samples (`-m budget`): the M5 page 13,293,641 bytes with the whole fonts and every trend point (59,727), now 10,120,283 with the subset and 44 trend points; the M6 page 6,689,232 to 5,003,367, its running score from 59,949 points to 195; cap 15,000,000. `ctx de` took 324.4 s before and 447.4 s after, the second under a load average near 21 from another worktree's tests, both inside the 600 s budget.
- **Mutations**, 81 planted and 81 caught, each stated in its test: 29 on the four requests, 37 on the review's issues (among them a scale of 1.5, the parallel sets back beside the guide, a pt literal in a template, the notices, the thinned points, the labels, the dashes, the colours' decimals, the href kept while moved, the narrow band map, the measured title drop), and 15 on this round's changes: the page's characters, the package's text, the margin, the name records, the fonts' notice or the whole fonts in place of the subset; the page line drawn whole, thinned in axis units or in the files too; the note left out, no stop added, numpy imported where the CLI starts, numpy text skipped by the walker. The whole fonts inlined passed its test at first, which measured the function rather than the page; the test now reads the page's own fonts. Three are caught by an error, not an assertion: the brush or the filter never reporting a move, and a throw in `draw()`, stop the browser fixture.

**Not done.** The Step-1 pages under `~/ctx_tier1` were not rebuilt, so their bytes after this change are not measured. A real lasso drag was not driven by these tests; the review drove one by DevTools mouse events and found both views painted and both paints in the export (2026-10-01).

**Amended 2026-10-02 (Rytis): the subset is renamed Ctx Sans, which settles the reserved-name question.** The open question above, embedding or a Modified Version, is settled by renaming. The subset is a Modified Version, which OFL 1.1 condition 3 forbids the Reserved Font Name Liberation; the licence's first lines also reserve Arimo, Tinos and Cousine for Google. The restriction covers the primary font name as presented to users.
- **The names.** `export.rename` writes them into the subset's name table at build time, on the Mac and Windows records alike: family Ctx Sans (name ID 1); full Ctx Sans Regular and Ctx Sans Bold (4); unique Ctx Sans Regular; Version 2.1.5 (3); PostScript CtxSans-Regular and CtxSans-Bold (6). Every other record keeps Liberation Sans 2.1.5's text: the copyright (0), the version (5), the trademark notice (7), the description (10), the licence and its URL (13, 14).
- **The page.** It prints the OFL text once, as before. The fonts' note (`vendor/MANIFEST.toml`) now says that Ctx Sans is a Modified Version of Liberation Sans 2.1.5, cut to the page's characters and renamed as the licence requires, under the OFL. The old note's "not for use outside this page" is gone: condition 5 puts a Modified Version under the OFL alone, with no added restriction.
- **The PDF.** jsPDF registers the face as `CtxSans` and each style under its PostScript name (`export.FACE`, `export.NAMES`). Its subsetter copies the `name` table raw, so each embedded program carries the renamed records too.
- The vendored files stay Liberation Sans 2.1.5 byte for byte; the manifest's hashes are unchanged.

Evidence (rename):
- `tests/test_export.py::test_the_subset_is_renamed_and_keeps_its_notices` reads the page's own two fonts. On both platforms' records: none of the four Reserved Font Names in IDs 1, 3, 4 or 6, nor in 16, 17, 18, 21, 22 or 25, which these fonts lack; family Ctx Sans and the two PostScript names; since review (2026-10-02) the unique and full names too, as literals: Ctx Sans Regular; Version 2.1.5 and Ctx Sans Bold; Version 2.1.5, Ctx Sans Regular and Ctx Sans Bold; every other record, copyright, licence and URL among them, the original's. The OFL's title and "with Reserved Font Name Liberation" appear on the page once each, and the note once.
- `test_two_builds_give_the_same_fonts`: two fresh interpreters with hash seeds 1 and 2 give the same SHA-256 for both fonts.
- In headless Chromium, `tests/test_export_browser.py::test_bold_text_embeds_the_renamed_faces_from_the_subset`: the PDF names CtxSans-Regular and CtxSans-Bold. Each embedded program, parsed by fontTools, has family Ctx Sans, the PDF's PostScript name, no Liberation in IDs 1, 3, 4 or 6, and Red Hat's copyright line; "Regular −2 ≥ ρ" and "Bold ρ" extract as text. `test_files_take_the_figure_files_size_and_embed_ctx_sans`: CtxSans-Regular is the only font of the linked, volcano and bars files.
- Bytes, the demo module page (linked view, volcano, bars), measured with the rename and with `rename` switched off: regular 73,136 against 73,092 bytes, bold 73,552 against 73,612, the page 5,893,854 against 5,893,874.
- Mutations, 9 planted and 9 caught: `rename` never called, caught by the name records and by the PDF's embedded programs; only the Windows records renamed; the name records cut to the subsetter's default; the fonts' note left off the page; the unique name stamped with the build's time; the whole fonts inlined; a manifest hash one digit off; the PostScript names not set in jsPDF, which named both faces CtxSans.

Pressed again after the rename (review, 2026-10-02), on module pages `ctx` built from the test fixtures with this branch and the review's fixes, by `tests/browser.mjs` in headless Chromium 152, each view moved at a 1,440 and at a 375 px window. A view's two PDFs have the same content stream and text, differing only in jsPDF's creation date; its two SVGs differ only in plotly's random clip ids and, in M1's, the place of an empty style attribute:

| Fixture page | Moved | PDF, mm | PDF, bytes | Fonts | Change in the files |
|---|---|---|---|---|---|
| M1, library sizes | cell_type type_a | 170 by 85 | 14,099 | CtxSans-Regular | 8 bars of 14, named as the page |
| M2, linked view | 6 samples orange, 3 vermillion | 170 by 153 | 16,583 | CtxSans-Regular | 6 and 3 points in the paints, both paints in the bands |
| M5, volcano | log2 FC 2, BH p below 0.01 | 170 by 93.5 | 26,868 | CtxSans-Regular | 48 up, 46 down, as the page |

  Each embedded program, 17,260 to 18,900 bytes, names itself Ctx Sans (ID 1), Ctx Sans Regular; Version 2.1.5 (3), Ctx Sans Regular (4) and CtxSans-Regular (6); "Liberation" is left only in the trademark notice, ID 7. Each PDF's text holds the live line and ends with the note; each page made no request but itself and no console error. The overview pages now measure 5,693,443 bytes (M1), 5,845,160 (M2), 7,383,413 (M5) and 4,994,627 (M6, by `ctx function` on the M5 root), against 5,692,220, 5,843,816, 7,382,162 and 4,993,457 in the byte table above, measured before the rename.

**Fixed after review (2026-10-02): a view moved on a narrow window exported the page's wrapped x title.** The page wraps an x-axis title wider than the room under its plot (`_page.js.j2`, titles), and `ctxExport.figure` copied the layout as it stood, so a view moved at a phone's width broke its 170 mm files' title where the phone's plot broke it. The titles guard now keeps the title as built in `layout.meta.xtitle`, and `ctxExport.figure` draws that one. Measured in the SVG of the demo volcano moved at a 375 px window: "Fold change (log2 fold" over "change)" before the fix, "Fold change (log2 fold change)" after, as from 1,440 px. On the fixture pages above, after the fix, the M2 and M5 views moved at 375 px export "PC1, 42% of variance (VST units)" and "log2 fold change, case over control" on one line while their pages wrap them. Tests: `tests/test_export_browser.py::test_a_view_moved_on_a_phone_keeps_the_files_x_title` exports the moved volcano at 375 and 1,280 px and compares the SVGs' x titles; `tests/test_page_script.py::test_export_draws_the_x_title_as_built` runs `ctxExport.figure` under node. Mutations: the restore removed, caught by both; the guard never keeping the title, caught by the browser test. Also after review: the rename test stated only the family and PostScript names, so a full name set to the family alone, "Ctx Sans" for both styles, passed it, as did a unique name without its version; both fail the test now. `ctx.viz.metrics` said a page's PDF embeds Liberation Sans; it names Ctx Sans.

---

## D-032 — Seven figure contracts for M2; the density and cividis half of Phase 2D's contour helpers moves into M2

**Date:** 2026-09-30 · **Status:** decided

**Context.** M2.5 asks for the figures the QC page draws: strips with bands, a bubble chart, curves, heatmaps, the read-count Sankey, and contours. Each one sets a number of drawing rules a reader could have set otherwise. The Phase 2D item "contour helpers" held the density grid, the highest-density levels and the cividis tokens that M2's contours need now. It also held click-to-set sliders, which belong with FigureSpec's Reproduce panel.

**Decision.** Each contract is one frozen dataclass read by both renderers, and each has the rules below.
- **`StripData`.**
  - It draws points per category, spread by the golden ratio, and flagged points as triangles.
  - A metric's absolute bands (`Metric.spans()`, the ranges `band()` assigns) are shaded in the band colours at 0.15 opacity.
  - Each stratum's median ± k MAD is a grey bar with its median line. The bar's numbers go into the CSV.
- **`ScatterData.size`.** Bubble area grows linearly with the value, from 3 to 14 pt in diameter, with three legend keys.
- **`CurveData`.** It draws one line per series, and each series keeps one colour.
- **`HeatmapData`.**
  - Every cell holds a level and every level has a swatch, so one contract draws the QC band map and the expression map.
  - `diverging_levels` gives [a, b) bins in the two role colours at opacities 0.33, 0.67 and 1.0, with a blank middle.
  - `cluster_order` is average linkage on 1 − Pearson, with constant rows last.
  - The row statistic is carried as a number and formatted only where it is printed. A formatted string had made the CSV and the page disagree.
  - Above 60 rows the labels are left off and the drawn height stops growing.
- **`FlowData`.**
  - The Sankey is drawn by hand as vector Bézier bands for the PDF and SVG, and with plotly's own Sankey on the page.
  - It refuses a node whose reads in and out differ, a backward link and a doubled link.
- **`SurfaceData`.**
  - It fills eight equal cividis bands. The stops are matplotlib 3.11.2's colormap at the band midpoints, kept as style tokens.
  - Its isolines are in ink, and each label sits on a box of the grid grey. The stops run from 13.1:1 to 1.4:1 against white, so no ink reads on all of them.
  - On the page, isolines are black in both themes, because they sit on the fill and the theme does not change it.
  - The current setting is marked.
- **`DensityData`.**
  - It draws the 50 % and 90 % highest-density regions from one smoothed histogram (64 bins, σ 2 bins, a 10 % margin) on a grid every group shares.
  - A group needs at least 30 points to get a density. Smaller groups and highlighted samples are drawn as points.
  - Each sample's region is in the CSV.
- **Shared rules.**
  - A log axis spanning fewer than two decades labels 1, 2 and 5 of each (`log_ticks`).
  - `axis_title` now prints a unit that is only a symbol: "Reads kept (%)". It had dropped `%`, because a unit with no word in it counted as already named. No existing figure used such a unit.
- **Phase 2D's contour-helper item, split.** The density grid, the highest-density levels and the cividis tokens are done here. The two-knob and power-map templates with click-to-set sliders stay in Phase 2D, with FigureSpec.

**Evidence.**
- **Tests.** Each contract has a demo in `tests/synth.py`, enrolled in the measured PDF/SVG test (fonts, painted colours, Tf sizes, line widths, the 85 or 170 mm width, caption and axis titles), the overlap test and the page round trip. The suite passes with `UserWarning` raised as an error.
- **What the new unit tests pin:**
  - the spans against `band()` at every edge;
  - bubble area linear in the value;
  - clustering recovering two planted blocks, with a constant row last;
  - each diverging bin at its edges;
  - the flow refusing lost reads and a backward link;
  - the cividis stops against the colormap;
  - `hdr_levels` against its definition: the cells at or above the level hold at least the mass, and those above it less.
- **Coverage.** At n = 60 and 45, the 50 % region of the smoothed density holds 41 of 60 and 32 of 45 points. Smoothing spreads the mass, and each sample's region is in the CSV. At n = 5,000 the levels meet their definition exactly.

**Consequence.** M2.6 builds the QC page from these contracts. LinkedData (P2Q.4) remains.

**Amended 2026-10-01 (Rytis): the band map says what a column is, and each row its unit and cut-offs.** The x-axis read "Sample, as nf-core/rnaseq names it." and the rows gave titles alone. `HeatmapData` gains `column_title` and `row_display`. In both renderers the column title reads, for example, "Each column one run: 5, in name order; hover names it", with the columns now in name order. Each row is the metric's title over its unit and cut-offs, built from `conf/qc_metrics.toml`, the words in `qc.toml`: "Processed reads (reads; good from 10,000)", "Duplicated reads (%; good to 30, usable to 50)", "Adapter content (status; pass good, warn usable, fail bad)". A row is as tall as its label's lines. On the page the column labels stand upright, as in the PDF, and the plot grows by their length. The long title is right-aligned under them, wrapped to a phone's plot, since centred under cells that long row labels push right it ran off the plot. `metrics.tg_total_reads` is titled "Processed reads"; its meaning and computation stay, and no band moved.

Fixed after review (2026-10-01): the title's drop under the column labels was estimated at the average glyph and fell 38 px short of 19-character run names in capitals, so it ran into them; it is now measured from Liberation Sans advance widths (`ctx.viz.metrics`). At 375 px the long second lines started up to 132 px left of the plot; below the width the labels need, the page script swaps in a narrow layout the figure carries (`layout.meta.narrow`): every row label wrapped to half a phone's plot, the plot taller to hold it. Hover named every cell's value "None"; its last line now names the band, and a row by its one-line label. Measured in Chromium 152 on a band map of the registry's 11 banded rows over five runs named in capitals (2026-10-02): every row label starts 8.1 px or more inside the plot at 1,280 and 375 px and 12 px opened alone, and the title clears the lowest column label by 7.8 px on the page and 15 px alone (`tests/test_export_browser.py::test_band_map_title_and_labels_stay_clear`). Before, at 375 px, only every third row was labelled and the title was clipped.

---

## D-031 — `ctx qc` opens the FASTQ files the sheet names; fastp is titled by run; read 1 stands for a pair

**Date:** 2026-09-30 · **Status:** decided

**Context.** M2a runs FastQC, fastp and seqkit on the reads. CONVENTIONS made `ctx gather` the only module that reads the fetchngs outdir's layout, and every later module reads the harmonised sheet. The reads themselves are in that outdir, and MultiQC names each tool's samples its own way: FastQC and seqkit by file, fastp by its input or its report title. A paired run appears as a `_1` row and a `_2` row for FastQC and seqkit, and as one row for fastp.

**Options.** Copy the FASTQ files out of the outdir first, doubling the disk a 12.6 GB teaching set takes, or open them in place through the paths gather wrote. For the names, guess MultiQC's cleaning from file names, or give each fastp run a title.

**Decision.** `src/ctx/qc/reads.py`:
- **Which files.** It opens only the files the harmonised sheet names, as `fastq/<name>` relative to `--fetchngs-dir`, and only for runs whose every file gather marked `verified` or `no_md5`. A `mismatch`, `missing` or `not_downloaded` file leaves its run unread, and the run is listed with that status. A readable file that is not on disk stops the run, listing every such file. CONVENTIONS says so beside the gather rule.
- **Names.** Each fastp run is titled with its run accession (`-R`), and MultiQC names fastp's sample from the title. FastQC's and seqkit's rows map by file stem, and any name that maps to no file stops the run (`multiqc_io.map_names`).
- **Which row a run takes.** Each column comes from the run's own row, else its read-1 row. Read 2 is not judged.
- **fastp's settings.** fastp runs without `-o`, so it filters and reports without writing reads. Its HTML goes to the null device. Its settings are fastp 1.1.0's defaults unless changed: `-q` 15, `-u` 40, `-l` 15, and `--cut_tail` off with `-W` 4 and `-M` 20.
- **The grid.** Each `-q` × `-l` cell is titled `grid_q<q>_l<l>` and read from MultiQC's General Statistics, with read counts kept raw.
- A metric's absolute band and its MAD band can both flag it. The status counts that metric once and keeps both reasons.

**Evidence.**
- **Tests.** `tests/test_qc_reads.py` runs the five chr22_se samples and fetchngs's paired test run through FastQC 0.12.1, fastp 1.1.0, seqkit 2.13.0 and MultiQC 1.35.

  | Check | Result |
  |---|---|
  | Every single-end value | equals the tools' own files: 1,000 reads; seqkit's Q30 read from its TSV |
  | The paired run | takes read 1's FastQC count (91) and seqkit Q30, and fastp's pair row |
  | Flags on the planted metrics | exactly se_c for adapter content, se_d for per-base quality and Q30 (65 %), and the paired run's read 1 for Q30 (70 %, fetchngs's own data) |
  | Grid cells, se_d | (15, 15) keeps 1,000 reads; `-q` 38 keeps none, as every base is at most Q37; `-l` 102 keeps none |
  | Mutations | three caught: q and l swapped, `mismatch` read, read 2 taken for read 1 |

  Two mutations were equivalent and not caught, because each module writes a column on only one of a run's rows: reversing the row order, and admitting read-2 rows after read 1.
- **Measured and to be taught, not fixed:** fastp 1.1.0's adapter auto-detection did not find the adapter in se_c. 276 of its 1,000 reads run into it, but no `adapter_cutting` section was written and the mean length stayed at 101 nt. FastQC failed se_c's adapter content.
- **A unit gap:** for a paired run, fastp's `passed_filter_reads` counts both mates (172 of 182 here). nf-core's `min_trimmed_reads` band counts pairs, so the registry's meaning now says a paired run meets the band at half the pairs. The course's cohorts are single-end.

**Consequence.** M2.6 does three things with this: it writes the per-run values and flags beside the count checkpoint's, gives the fastp settings and the grid their cards, and draws the grid as a surface (M2.5).

---

## D-030 — The checkpoint on counts: bands per stratum, agreement per peer group, confounding on M5's own selections

**Date:** 2026-09-30 · **Status:** decided

**Context.** M2.3 judges every sample on its counts before M5 fits anything, and asks whether the judgement is confounded with the comparisons M5 makes. The registry marks all ten count metrics `dataset_only`, so each is judged against its own dataset's median ± 3 scaled MAD. Two things had to be chosen: which samples a band is computed over, and which samples a metric is tested between.

**Options.** Bands over the whole dataset (neutrophils and monocytes differ in genes detected by biology, so whole cell types would be flagged), per stratum, or per stratum and condition (a condition that shifts a metric would then never show as flags). For the confounding test, every condition pair in each stratum, or exactly the samples each registered contrast compares.

**Decision.** `src/ctx/qc/checkpoint.py`, in three steps:
- **Metrics.** Library size, genes detected (count ≥ 1), top-gene share, and four disjoint composition classes from the GENCODE gene types of the run's own GTF. The classes are: Mt_rRNA; the other chrM genes; nuclear rRNA and rRNA_pseudogene; and HBA1, HBA2 and HBB. A counted gene absent from the annotation stops the run. The sex check is `ctx.de.run.sex_check`, imported.
- **Agreement with the peer group.** The peer group is the stratum and the condition:
  - replicate agreement is the median Spearman correlation with each other peer, over the genes read in any of them;
  - PCA distance is the distance from the peers' median on PC1 and PC2 of one blind VST of every sample, the PCA M5 draws.
  - A sample with no peer has neither.
- **Bands.**
  - **Where computed.** Per stratum, not per condition: a condition that shifts a metric then shows as flags in one group, and the confounding test catches it rather than the band hiding it.
  - **Which side flags.** Only the bad side, from `direction`; both sides for `none`.
  - **When a stratum gets no band.** It needs at least `min_samples = 3` samples, now in `[registry]`: with two, each lies exactly one MAD from their median, so neither could ever be flagged. A MAD of 0 gives no band either, and says so.
  - **Metrics that never flag.** A metric with action keep is judged and never flagged.
  - **The judge is generic.** It also handles banded, status, category and strandedness metrics. A strand call flags when it is undetermined or differs from its stratum's usual call.
- **Status.** One row per run, from its sample's flags: pass, flag, or exclude where a flagged metric's action says exclude. None does, so nothing is excluded. Each reason is filled from `content/qc.toml`.
- **Confounding.**
  - **Which samples.** Each unpaired contrast's samples are exactly those of `ctx.de.deseq.select_samples`, one run per individual in each group. Contrasts that select the same samples are tested once, and the row names them all.
  - **Which tests.** Numeric metrics by Mann-Whitney through `annotate.compare_groups`, with the rank-biserial r. The other metrics' flags, and any flag at all, by Fisher's exact test with the odds ratio of tested over reference. Haldane's 0.5 is added only when a cell is 0.
  - **Correction.** BH runs over every test in the table.
  - **What is left out.** A paired contrast gets one row saying why it is not tested: its groups share individuals, and both tests assume they do not.
- `strata` and `peers` are arguments, not cohort logic. The CLI and its parameter cards come in M2.6, where Step-1 passes tissue × cell type.

**Evidence.**
- **The GENCODE 50 GTF**, `gene_types` on `gencode.v50.primary_assembly.annotation.gtf.gz`, took 6.0 s at 211 MB peak RSS. Every count below agrees with an independent `awk` count where one was run:

  | Class | Genes |
  |---|---|
  | All gene records | 78,941 (awk: 78,941), IDs all unique |
  | chrM | 37 (awk: 37): 2 Mt_rRNA, 22 Mt_tRNA, 13 others |
  | Nuclear rRNA | 53 (awk: 53) |
  | Nuclear rRNA_pseudogene | 497 |
  | HBA1, HBA2, HBB | all three present |

- **A correction to D-022:** the file holds no gene ID ending in `_PAR_Y`. D-022's "so gene IDs including `_PAR_Y` match the counts" was never measured. The IDs match because it is the file the cluster run used.
- **Tests.** `tests/test_qc_checkpoint.py` has 12 tests. Every value in the count test is computed by hand. Six mutations were each caught:
  - Mt_rRNA counted in the mitochondrial class;
  - no merging of identical selections;
  - Haldane's correction applied always;
  - both sides flagged;
  - a keep metric flagged;
  - a strand call not compared with its stratum's.

**Consequence.** M2.6 does four things with this:
- it writes the metrics, bands, status and confounding tables;
- it feeds the per-sample MultiQC values (read 1 rows, never the summed group rows) into the same `judge`, which is how the strandedness values arrive;
- it gives `--strata` and `--peers` their cards;
- its page states the 45S gap.

M5 joins the status table after the Step-1 record.

---

## D-029 — The QC registry: 30 metrics, their bands at their sources; MultiQC coloured by the same bands, on raw counts

**Date:** 2026-09-30 · **Status:** decided

**Context.** CONVENTIONS keeps QC bands in `conf/qc_metrics.toml` alone, each with a source and a source type, and M2's spec has MultiQC's General Statistics table coloured by them. Most RNA-seq QC metrics have no universal cut-off.

**Decision.**
- 30 metrics, across four stages:
  - reads: 13, from FastQC 0.12.1, fastp 1.1.0 and seqkit 2.13.0;
  - trimming: 4, from Trim Galore 2.3.0;
  - quantification: 3, from Salmon 2.7.0;
  - counts: 10, computed by ctx.
- Absolute bands only where a source gives one, 10 metrics in all:
  - FastQC's duplication limits from its `limits.txt`, and its own pass, warn and fail for six modules (`tool_code`);
  - nf-core/rnaseq 3.27.0's `min_trimmed_reads` 10,000 and strandedness thresholds 0.8 and 0.1 (`pipeline_default`; the 10,000, set for Trim Galore, labelled as transferred to fastp);
  - raw Q30 at 80 %, a course convention with its rationale.
- Every other metric is `bands = "dataset_only"` with a rationale, judged against its own dataset by median ± 3 scaled MAD. The rRNA metric says GENCODE 50 annotates no 45S gene, so it is a lower bound; mitochondrial rRNA is counted apart from the other chrM genes, so no read is counted twice.
- `ctx.qc.registry` loads the file in `ctx.dataset`'s pattern, every problem reported at once. It writes the MultiQC config: band colours from `style.BAND_COLOURS` (Okabe-Ito), one conditional rule per banded General Statistics column, a section comment per module from the metrics' meanings, and read and base count multipliers of 1.

**Evidence.**
- **Sources, read on 2026-09-30:** FastQC 0.12.1's installed `Configuration/limits.txt`, whose values match TODO.md's; `fastp --help` for 1.1.0's defaults (-q 15, -u 40, -l 15, window 4, mean 20, -n 5); the rnaseq fixture's params file for nf-core's thresholds.
- **MultiQC's rule semantics, measured on 1.35:** classes apply in the order pass, warn, fail, and the last that matches wins (a Q30 of 65 matched warn and fail and took fail's colour). Conditions in one list are alternatives. `tests/test_qc_registry.py` applies the generated rules that way and gets the registry's band for every value around every edge.
- **The count trap:** General Statistics scales read counts to millions, in its TSV too, and colours the scaled value. The fixture's Trim Galore totals of 50,000 to 200,000 reads read as 0.1 and were painted bad against the 10,000 band. With the multipliers at 1 they are raw and painted good; the real-run test checks the colour.
- **Also found:** General Statistics adds a summed row for each paired sample (100,000 for 50,000 pairs), so per-sample values are read from read 1's row, not the group row.

**Consequence.** `ctx qc` always runs MultiQC with the generated config. The counts-stage metrics are computed in M2.3. Band contexts (library type, read length, platform) wait until a source gives a band that differs by context.

---

## D-028 — MultiQC's output folder keeps MultiQC's names; `ctx.qc.multiqc_io` reads it

**Date:** 2026-09-30 · **Status:** decided

**Context.** CONVENTIONS makes MultiQC the only parser of tool QC logs, and M2 shows MultiQC's own report beside the atlas (TODO.md, M2). MultiQC names its files itself (`multiqc_report.html`, `multiqc_data/multiqc_<module>.txt`), which the naming lint rejects under `outputs/`.

**Options.** Rename MultiQC's files after each run, which breaks its report's links to its own data; or exempt its folder, as AGENTS.md exempts nf-core's `results/`.

**Decision.** `ctx.io.naming.NATIVE_DIRS = {"multiqc"}`: files inside a folder of that name are not linted; every other file under `outputs/` still is. `ctx.qc.multiqc_io` runs MultiQC 1.35 by subprocess from the `reads` environment into `outputs/02_qc/multiqc/`, reads its per-module tables as text, maps each sample name through an explicit alias table and refuses any name it cannot map, and stops a run when an expected sample has no source report (D-022).

**Evidence.** `tests/test_multiqc_io.py` on the rnaseq fixture: MultiQC's `multiqc_trim_galore.txt` equals the JSON reports read without MultiQC, field by field; every report's total equals written + too short + too long + too many N + discarded untrimmed, exactly; the reads Trim Galore wrote (read 1 for paired samples) equal Salmon's `num_processed` for all five samples; MultiQC drops `_trimmed` from single-end names and keeps `_trimmed_1`/`_2`, which the alias table covers; a report removed from the sources stops the run naming the sample; the lint passes the folder and still flags a badly named file beside it. One MultiQC run takes about 4 s.

**Consequence.** Every M2 table that comes from a tool log comes through `multiqc_io`.

---

## D-027 — The read-level tools live in their own pixi environment, `reads`; fastp 1.1.0 and samtools 1.22.1

**Date:** 2026-09-30 · **Status:** decided

**Context.** Track A (M2-M4) and `ctx qc` run FastQC, fastp, seqkit, MultiQC, STAR, samtools, featureCounts, Salmon, Trim Galore and HISAT2 on laptops as well as on the cluster. The design review found that adding them to the analysis environment would move it: MultiQC 1.35's bioconda recipe caps plotly below 7 (the lock holds 7.1.0) and switches polars to its compat runtime, and htslib 1.21 needs a libdeflate the locked libtiff excludes. The analysis environment is what `ctx de`, `function` and `report` run in, after the Step-1 rules were registered.

**Options.** One shared solve, moving plotly, polars and libtiff; or a separate environment, as `fetch` already is, with ctx running the tools by subprocess.

**Decision.**
- A `reads` feature and environment in `pixi.toml`, both platforms, outside the `main` solve group. `ctx.tools` finds each binary in `$CTX_TOOLS`, else the sibling `.pixi/envs/reads/bin`, else PATH, runs it with its own bin folder first on PATH, and refuses any version other than the pin; `tests/test_tools.py` holds its table equal to `pixi.toml`.
- Versions are nf-core/rnaseq 3.27.0's module pins, except two. On osx-arm64 every STAR 2.7.11b build needs htslib 1.21 or 1.22.x (libdeflate 1.20 or 1.22), while fastp 1.2.0 and later need libdeflate 1.25, so no solve holds nf-core's fastp 1.3.6 beside STAR: fastp is pinned at 1.1.0, the newest build on libdeflate 1.22, and samtools and htslib at 1.22.1 with STAR's `_8` build. The cohort was trimmed with Trim Galore, so fastp's version changes no comparison with it. seqtk stays instructor-side, on the cluster.
- CI installs `reads` beside `default`; a test that needs a tool fails with the install command rather than skipping.

**Evidence.**
- `pixi lock`: every existing environment's package list is byte-identical before and after, both platforms (default 337 osx-arm64 and 389 linux-64 packages, laptop 313 and 364, prod 214 linux-64, fetch 29 and 77); `reads` adds 132 packages on osx-arm64 and 165 on linux-64, with its own plotly 6.9.0.
- Rejected solves, as pixi reported them: STAR with samtools 1.21 and fastp 1.3.6 (libdeflate 1.20/1.22 against 1.25); samtools 1.22.1 with fastp 1.3.6 (1.22 against 1.25).
- `pixi install --locked -e reads`: 13.6 s, 1.6 GB on disk. Every tool reports its pin: FastQC v0.12.1, fastp 1.1.0, seqkit v2.13.0, STAR 2.7.11b, samtools 1.22.1, featureCounts v2.1.1, salmon 2.7.0, trim_galore 2.3.0, hisat2 2.2.3, multiqc 1.35.
- Through `ctx.tools` from the `default` environment, which holds no Java: FastQC wrote its report and MultiQC parsed it, per-base quality table included.

**Consequence.** Every Track A module and `ctx qc` run their tools through `ctx.tools`. pytximport, a Python library, is decided with M4.

---

## D-026 — The Step-1 record is written by a script: the hub's verdicts, byte for byte, and the M5 and M6 metrics together

**Date:** 2026-09-30 · **Status:** decided

**Context.** The design review for the Step-1 record found no writer for the files the record tests read: `tests/test_sanity_record.py` globs `deliverables/*_de_tier1_sanity_v*.csv` and `*_de_tier1_metrics_v*.csv`, while `ctx report` writes `outputs/09_report/<date>_report_<slug>_sanity_v<N>`. `test_decision_follows_rule` decides every rule again from the metrics CSV alone, and rules R4, R5 and PDAC read `gsea.*` metrics that only M6 writes, so a copy of M5's metrics alone would re-decide them as FAIL against the recorded PASS. If either file is missing, the four record tests skip silently.

**Options.** Change the record tests to read both modules' metrics tables, or write the files the tests already read. `tests/` is not frozen, but the tests are the registered contract; a writer changes nothing that was registered.

**Decision.** `scripts/python/record_sanity.py` writes both files, CSV with README. It refuses a second record, a metric written by both modules, a hub CSV that is not the CSV of its Parquet, a hub verdict table whose README does not list the metrics tables as they are now (a rerun of M5 or M6 after the hub), and a metrics union that decides differently from the hub; in each case it writes nothing. The verdict CSV must be byte-identical to the hub's, or what was written is removed.

**Evidence.** `tests/test_record_sanity.py` runs the frozen `ctx.report.hub.sanity` on synthetic M5 and M6 metrics and the synthetic rules S1-S3: the record's verdicts equal the hub's bytes; deciding again from the record's metrics gives S1 PASS, S2 PASS (only because the M6 rows are present) and S3 FAIL (no module writes its metric); each refusal fires and leaves nothing behind. Reducing the union to M5 alone fails 4 of the 6 tests. No frozen file changed.

**Consequence.** P2C.3 runs the script after `ctx report`; the freeze test is run before it, because it skips once the record exists.

---

## D-025 — The Trim Galore reports leave the cluster as JSON; the cluster run publishes them; registered again

**Date:** 2026-09-30 · **Status:** decided

**Context.** D-022 copied `trimgalore/*_trimming_report.txt` for the read fates of `ctx qc`. Before any M2 code, the design review ran MultiQC 1.35 on those files: it built no Trim Galore table. Its search patterns say why. The cutadapt module skips a report reading "Trim Galore version: 2 or later" (`exclude_contents_re`, commented "TG v2+ has custom JSON reports"), and the `trim_galore` module reads only `*_trimming_report.json`. nf-core/rnaseq 3.27.0 publishes only TRIMGALORE's `*.txt` (`conf/modules/qc_trim_filter.config`, third `publishDir` entry), so the JSON never reached `results/`. With D-022 as it stood, the read-fates chart Rytis asked for would have had nothing to draw.

**Options.** Parse the text reports in ctx, against CONVENTIONS' rule that MultiQC is the only parser of tool logs; or publish the JSON on the cluster and copy it. The cluster run had not started, so the second was still possible.

**Decision.**
- The cluster run's `rnaseq_kebnekaise.config` repeats TRIMGALORE's three `publishDir` entries, the last widened to `*_trimming_report.{txt,json}`. `scripts/python/make_rnaseq_fixture.py` holds the same text as `PUBLISH`, and `test_repo::test_publish_override_single_source` holds the two equal.
- `COPIED_FILES` and step 7's rsync take `trimgalore/*_trimming_report.json` in place of the `.txt`, which nothing reads.
- The rules are registered again, at 2026-09-30 00:04.

**Evidence.**
- **MultiQC 1.35 on the committed `.txt` fixture:** modules salmon, cutadapt and trim_galore requested; only `multiqc_salmon.txt` written.
- **The override, on the P2A.2 rehearsal:** re-run with `-resume` and the override into a new outdir (`out_json`, so its params file records that outdir): 117 s against the first run's 705 s, every task cached. Its `trimgalore/` holds the 8 `.txt` and 8 JSON reports; `salmon/` and the `.txt` reports are byte-identical to the first run's. MultiQC 1.35 on that output writes `multiqc_trim_galore.txt` with `total_reads`, `reads_with_adapter`, `reads_written`, `reads_too_short`, `reads_too_long`, `reads_too_many_n` and the pair columns (WT_REP1 read 1: 100,000 total, 99,566 written, 434 too short).
- **The fixture,** rebuilt from that run with `--from-run`: the 8 JSON reports replace the 8 `.txt`; the params file and trace are the resumed run's; the other 14 files are byte-identical. `test_rnaseq_fixture`, `test_rnaseq` and `test_repo` pass; changing one character of the override in TODO.md fails the single-source test.
- **No Step-1 count exists:** `~/ctx_tier1/data/raw` held only `tier1_meta` at 2026-09-30 00:04.

**Consequence.** Cluster step 5 must use the new `rnaseq_kebnekaise.config`; step 7 waits for this commit. The M2 read fates come from `multiqc_trim_galore.txt`.

**Amendment, 2026-09-30: the override published symlinks.** On the cluster, job 47485533 published every Trim Galore report as an absolute symlink into `work/`, and a mid-run copy left 117 dangling links on the Mac. The run's params file says `publish_dir_mode = copy`; no launch-directory or `$NXF_HOME` config exists. The cause is the override's `mode: params.publish_dir_mode`, read when the `-c` file is parsed: under `NXF_SYNTAX_PARSER=v2` a param set only in the pipeline's own config is not visible there, so the mode is empty and Nextflow's default, symlink, applies. The `path` and `saveAs` closures read their params at task time and were right. Reproduced with Nextflow 26.04.6 on a one-process pipeline:

| Parser | Override reads | Params file sets the mode | Published as |
|---|---|---|---|
| v1 | `params.publish_dir_mode` | no | run refused: unknown config attribute |
| v2 | `params.publish_dir_mode` | no | symlink |
| v2 | `params.publish_dir_mode` | yes | file |
| v2 | literal `'copy'` | no | file |

The fix is `publish_dir_mode: copy` in step 5's params file, one definition that every config file sees; `tests/test_repo.py::test_rnaseq_override_reads_only_given_params` requires every `params.X` the override reads outside a closure to be set there. Job 47485533 keeps its links: step 7 copies their targets (`rsync -L`), and its `work/` must stay until its `results/` is copied out with links followed.

---

## D-024 — Figures take the Candida Ribo-seq style; the Helvetica and 14 pt overrides are dropped

**Date:** 2026-09-29 · **Status:** decided

**Context.** Rytis asked for the figure styling and parameters of the Candida FINAL figures (`riboseq_public_candida/deliverables/FINAL/figures/pdf`). This project had overridden the global figure rules with Helvetica and a 14 pt base (TODO.md, Phase 1 table); the Candida figures follow the global rules. The same day he asked for the ALI report's page design and the riboseq status report's split of page text into files (Phase 2S).

**Options.** Copy Candida's mechanism too, one TOML read by one module with an environment override, or change the values in `style.py`, which is already the one place a size or colour is set. The TOML adds a loader and its validation for one reader; `style.py` it is.

**Decision.** Taken from Candida, measured on its 32 PDFs and its `settings/figure_style.toml`:
- Arial, Liberation Sans, DejaVu Sans; 12 pt panel letters (bold capitals), 10 pt axis titles, 9 pt ticks, legends, annotations and captions, nothing under 8 pt.
- 3 pt ticks pointing out, minor ticks at the axis width; 6 pt points at alpha 0.5; grid #E6E6E6 at 1 pt (#39424A dark, from Candida's own dark scripts).
- No seaborn `paper` context: every rc is set in `style.py`.

Not taken, and why:
- Candida's light-only output: CONVENTIONS requires `_white` and `_black`.
- Its condition pair (blue and purple) and its domain colours: here the role pair reference #0072B2 and comparison #D55E00 matches the volcano's down and up.
- Free page heights: heights stay derived from aspect ratios.
- Creation timestamps: from Phase 2D figure files carry none.
- dpi 1000: only PDF and SVG are written, so it would change nothing.
- A panel-letter helper: no figure has panels yet; the style is set in `PANEL_LABEL`.

Two defects found on the way, fixed at their source:
- **Log-axis labels.** Matplotlib labels a log axis in mathtext, whose exponent is 0.7 of the tick size (8.4 pt at 12 pt, 6.3 pt at 9 pt) and whose SVG names one font instead of the chain. The MA and dispersion figures carried both. `scatter._plain_log_ticks` labels 100, 1,000, 10,000. The measured-figure test had no log axis; an MA contract joins it, and the SVG check is now equality per text element.
- **Screen sizes (mine, 2026-09-28).** `on_page` added 4 to a size in pt and used it as px, so the plots on a page ran 0 pt larger, not 3. It now returns `(pt + 3) * 96/72` px, serves the Plotly templates only, and a test holds the numbers; page text is set in px on its own (18 px base, unchanged).

**Evidence.** Every figure of the synthetic M5 run, 20 files in both themes, measured from the PDF and SVG before and after:

| Measure | Before | After |
|---|---|---|
| PDF font | Helvetica | ArialMT |
| PDF text sizes (Tf) | 8.4, 12, 14 | 9, 10 |
| PDF line widths | 0, 0.8, 1, 1.5 | 0, 1, 1.5 |
| SVG font-family | the chain, and `'Helvetica'` alone on log labels | `'Arial', 'Liberation Sans', 'DejaVu Sans', sans-serif` only |
| Page widths | 240.945, 481.890 pt | the same |

The M6 figures (10 files) measure the same as M5's after. Layout constants measured at 14 pt, measured again: at 9 pt a box plot fits 4 features in 85 mm and overlaps at 5 (`SINGLE_COLUMN_MAX_FEATURES` 3 to 4); 20 horizontal bars stop overlapping from 3.0 mm per row, and `BAR_ROW_MM` goes from 7 to 4.5 mm, the 4.1 mm line pitch of a 9 pt label plus a gap. Mutations that each fail a test: the log labels back to mathtext, a 14 pt annotation, a two-font SVG chain.

**Page design (P2S.2), same day.** The pages take the ALI modular report's design (`scAnalysis/knowledge/2026-06-25_ALI_report_modular`): its page, surface and inset colours, hairlines, 14 and 9 px radii, two-layer shadow, frosted sticky top bar with a three-colour dot, eyebrow kicker, monospace section numbers, stat tiles, table wraps, tinted callouts and bordered `details`. Changed from ALI, and why:
- Its accent #0072B2 is a data colour, so links keep #1D4E89 / #8CB8F0.
- Its faint grey measures 3.14:1 on white and 3.00:1 on its inset: not used for text. Every other text colour measures 4.5:1 or more on every surface it sits on in both themes; the lowest is the warn label in its callout, 4.58:1.
- Its pills put white on #009E73, about 3.4:1 at 12.5 px: the verdict pills keep their own colours (6.83:1 or more).
- Its secondary text at 13 to 15 px would undo the on-screen increase Rytis asked for on 2026-09-28, so the sizes stay (16 and 17 px secondary, 18 px body).
- A skip link and focus rings are added; ALI has neither. Pages stay offline; ALI loads Plotly from a CDN.

Checked in the browser on the synthetic M5 page and hub: 1280 px light, the toggle clicked to dark (page #101316, plot text #F2F2F2, choice stored), 375 px dark with no horizontal scroll, the skip link on the first Tab, no console errors, no request beyond the page. Each page grows by about 5.6 kB (M5 overview 6,636,384 to 6,642,016 bytes).

**Page text (P2S.3), same day.** The riboseq status report keeps every sentence in content files that an engine assembles, with a test that the page holds its load-bearing strings. Here: one TOML file per page in `src/ctx/report/content/`, plus `ui.toml` (template words and the skeleton headings) and `cards.toml` (option meanings, also the `--help` text), read by `ctx.report.text`. Strict: a missing key or an unfilled placeholder raises, naming the file and key; `tests/test_text.py` fails on a key no code names; `CTX_CONTENT` points at an edited copy with no fallback to the package. A sentence that varies with the data is two keys, and the code picks one. Moved now: setup, gather, the cards and every template word. Not moved: M5, M6, the hub and the captions, whose code is frozen until the record (Rytis's choice), as P2S.4. Mutations that each fail a test: a key nothing reads, a placeholder the code does not fill.

**Consequence.** No frozen file changed; `annotate.py` reads the new sizes through the names it already imports. Phase 1's records of the 14 pt figures stay as history.

**Open alone, at twice the size (amended 2026-10-01, Rytis).** "Open alone" drew a figure as small as the module page did. A standalone figure is now the page figure at `style.ALONE_SCALE` = 2, a uniform zoom:
- `style.zoom` doubles every length on screen in its plotly JSON: layout height, margins, every font, marker and line width, tick lengths, node sizes and annotation shifts; data, positions and shares stay.
- Its width is twice the width the module page gives it at full layout, `style.PAGE_PLOT_PX`: 686 px beside the reading guide and 1,022 px across the card. They are computed from the layout tokens the stylesheet reads; `--fig-pad` and `--fig-gap` are new tokens for the card's padding and the gap beside the plot.
- The theme toggle swaps templates at the same scale, and the legend's key length (`style.PLOTLY_LEGEND_KEY_PX`) is in them, so it doubles too. The band map's legend swatches, which had no size, take the page's marker size.
- The standalone page takes the window's width and puts the reading guide under the plot. A plot wider than its card scrolls inside it, under a line saying so (`ui.toml` figures.scroll), shown only while it scrolls.

What the zoom cannot reach: plotly draws a marker's legend key at most 16 px across and a bar's at 12 px, so alone the page's 10 px keys reach 16 px and a bar's stays 12 px; key lines double, 30 to 60 px long and 1 to 2 px thick. The trade-off: at twice the width a plot is wider than its card on most laptop windows. Measured in Chromium 152 on the demo pages (2026-10-02), the 1,372 px plot against its box: 1,167 px at a 1,280 px window, 1,327 px at 1,440, fitting from 1,512 (box 1,399); the parallel sets, 2,044 px alone, scroll at every width up to 1,600. Whether the 200 % should keep type, marks and lines at twice the size within the card's width instead is Rytis's call (review, 2026-10-01).

Evidence: `tests/test_export.py::test_alone_is_the_page_figure_at_twice_its_size` compares every length of five figure types' standalone JSON with the page's, and the templates, against a literal 2, Rytis's number, never `ALONE_SCALE` itself; `ALONE_SCALE` 1.5 had passed every test that read the constant as its own oracle (review, 2026-10-01). In Chromium at a 1,440 px window, the M2 fixture page drew the scatter 686 px wide with 20 px ticks and 10 px markers, the parallel sets 1,022 px wide; opened alone, 40 px ticks and 20 px markers. `tests/test_export_browser.py` measures the same on the demos, the legend keys and the scroll line among them. Mutations caught: the templates left unzoomed; the width left to the window; marker sizes left out of the zoom; `ALONE_SCALE` 1 and 1.5; the card's border left out of the width; the legend's key length left out of the template; the band map's swatches at plotly's 6 px; the scroll line never shown.

**Capped at the window (amended 2026-10-02, Rytis: cap the width, keep the 2x type and marks).** His answer to the trade-off above. Opened alone, every length `style.zoom` doubles stays doubled, and a plot is twice the page's width or its box's, whichever is less (`figures.figure_html`). The plot's box takes the CSS width `min(100%, 2 x PAGE_PLOT_PX)` and the layout no width, so plotly draws it at its box's width and redraws it as the window changes (`responsive` in `PLOTLY_CONFIG`). The cap changed four things:
- **The scroll hint.** The page script checked it when its box resized. plotly redraws a narrowed plot a moment later and resizes no box, so the hint stayed over a plot that fit. It is checked again on `plotly_afterplot`.
- **The parallel sets' least width.** `layout.meta.min_width` was doubled whole by the zoom, plotly's own 16 px bars, 40 px insets and 5 px label offsets with it: 1,572 px for the demo, where twice the type needs 1,392. It is now `meta.least_width`, pairs of px that grow with the type and px that do not (`LinkedData._least`); at zoom s the plot needs the largest s x scaled + fixed, 786 px on the page as before.
- **A heatmap's narrow layout.** `meta.narrow` holds the width below which the page wraps the row labels and the height it then sets, both at the page's size. Opened alone at a 768 px window, 7 of the band map's 11 row labels started up to 198 px left of its plot. `style.zoom` now doubles both (`style.NARROW_LENGTHS`).
- **A floor.** No plot opened alone is narrower than twice a phone's plot, 2 x `PHONE_PLOT_PX` = 634 px, the narrowest a page lays a figure out for. Capped alone at 375 px, the band map's 11 row labels all started left of its 317 px plot, by up to 103 px, and its title left the page. Below the floor the plot scrolls inside the card under the hint.

Measured in Chromium 152 on the demo pages (2026-10-02), the plot drawn against its box, in px:

| Window | Volcano alone | Scatter alone | Parallel sets alone |
|---|---|---|---|
| 1,280 | 1,167 in 1,167 | 1,167 in 1,167 | 1,392 in 1,167, scrolls under the hint |
| 1,440 | 1,327 in 1,327 | 1,327 in 1,327 | 1,392 in 1,327, scrolls under the hint |
| 1,512 | 1,372 in 1,399 | 1,372 in 1,399 | 1,399 in 1,399 |
| 2,560 | 1,372 in 2,447 | 1,372 in 2,447 | 2,044 in 2,447 |

Tick labels are 40 px alone and 20 px on the module page; the parallel sets' titles 43.3 and 21.7 px. The page itself scrolls sideways at none of these widths. Opened at 2,560 px and narrowed to 1,280, the volcano is redrawn at 1,167 px with no hint. The band map alone keeps every row label inside its plot at 768 and 1,440 px, and at 375 px is drawn 634 px wide under the hint.

**Open for Rytis (P2Q.13): the parallel sets below a 1,505 px window, and the floor on a phone.** "The parallel-sets figure fits too" holds from a 1,505 px window, where the card holds the 1,392 px of their least width, not at 1,280 or 1,440 (measured at 1,500, 1,504, 1,505 and 1,506 px). At twice the type the demo's "steady_state" is 227 px wide, and plotly sets a middle dimension's levels left of its bar, so the step between bars needs 280 px. Narrower, the parallel sets keep their least width and scroll inside the card under the hint, the rule they keep on a phone. Drawn at the card's width instead, levels run over bars only in the narrower of these windows: this paragraph said they would at 1,280 and 1,440 alike, and at 1,440 none does (corrected after review, 2026-10-02; the table below). A third way sets their type below twice the page's where the window is too narrow, which breaks "every text size stays 2x". The floor departs from the decision too: at a 375 px window every plot opened alone is drawn 634 px wide in a 317 px card, the linked view's parallel sets 1,392 px, and each scrolls inside its card under the hint, where the decision reads "nothing scrolls sideways"; the page itself does not scroll. Which one is his call.

The demo's parallel sets opened alone, measured in Chromium 152 (review, 2026-10-02): drawn at their card's width, the overlaps among dimension titles, levels and bars, and the least gap from a level to the bars of the dimension before it, negative where it runs over them; and as built, at their least width.

| Window | Card, px | At the card's width | As built |
|---|---|---|---|
| 1,280 | 1,167 | a level over a level and a level over a bar; −24.1 px | 1,392 px under the hint; 32.2 px |
| 1,360 | 1,247 | a level over a bar; −4.1 px | the same |
| 1,380 | 1,267 | none; 0.9 px | the same |
| 1,440 | 1,327 | none; 15.9 px | the same |
| 1,505 | 1,392 | none; 32.2 px | 1,392 px, no hint; 32.2 px |

Measured at every 20 px from 1,280 to 1,500, the least gap at the card's width grows 5.0 px a step. The least width keeps `PLOTLY_PARCATS` gap_px, 16 px at the page's type, twice over at twice the type: 32.2 px measured. At 1,440 the card's width keeps 15.9 px, about the 16 px the page's least width keeps at the page's type.

Evidence: `tests/test_export.py::test_alone_is_the_page_figure_at_twice_its_size` reads each plot box's CSS width off the standalone page. In headless Chromium, `tests/test_export_browser.py`: `test_alone_keeps_twice_the_type_within_the_window` at 1,280 and 1,512 px for the volcano and the scatter; `test_alone_takes_twice_the_page_width_in_a_wide_window` at 2,560; `test_alone_reflows_when_the_window_narrows`; `test_alone_parallel_sets_fit_where_their_least_width_allows`, which also finds no title or level over another or a level over a bar; `test_band_map_title_and_labels_stay_clear` alone at 768 and 375; `test_alone_on_a_phone_scrolls_at_twice_a_phone_plot`. `test_alone_says_when_the_figure_scrolls` is gone: at 1,280 the volcano no longer scrolls, and the hint's test moved to the parallel sets. Mutations, 11 planted and 11 caught: the width left to the window alone; the width fixed at twice the page's; the type left at the page's; `ALONE_SCALE` 1.5; plotly's responsive redraw off; the hint checked only on its box's resize; the least width doubled whole; the least width dropped alone; the hint never shown; the narrow layout left unzoomed; no floor.

**Fixed after review (2026-10-02): the band map redraws as the window changes.** "Capped at the window" above says plotly redraws a plot opened alone as the window changes; it did not redraw the band map, on its module page either. Its narrow layout sets a height through `Plotly.relayout`, and plotly.js 4.1.1 then turns autosize off and fixes `layout.width` at the width it drew; its own resize handler does nothing for a plot with both a width and a height. The page script swaps that layout in or out once at load, so every band map kept the width it loaded at. The relayout now carries `autosize: true` and a null width (`_page.js.j2`, narrow). Measured in Chromium 152 on the band map of the registry's 11 banded rows over five runs named in capitals, the width plotly drew in the width its div was given, in px:

| Page, window then window | Before | After |
|---|---|---|
| Alone, 375 then 1,280 | 634 in 1,167; 7 of 11 row labels start 201.3 px (six) and 219 px left of the plot | 1,167 in 1,167 |
| Alone, 1,280 then 2,560 | 1,167 in 1,372 | 1,372 in 1,372 |
| Alone, 2,560 then 375 | 1,372 in 634 | 634 in 634 |
| Module page, 375 then 1,440 | 317 in 686; 8 of 11 row labels start 23.8 to 140.5 px left of the plot | 686 in 686 |
| Module page, 1,440 then 375 | 686 in 317 | 317 in 317 |

After the fix no row label starts left of the plot at any step. `tests/test_export_browser.py::test_band_map_reflows_as_the_window_changes` opens the band map alone at 375 px and resizes the window to 1,280, 2,560 and 375 px, through the floor, the box and the cap in turn, and the module page from 375 to 1,440 px; at each step the plot is drawn at the lesser of its box and 1,372 px and never under 634 px alone, at its box's width on the module page, with every row label inside it. The driver, `tests/browser.mjs`, takes a list of resize steps, and `SETTLED` waits until plotly draws each plot at its div's width. Mutations: the relayout without autosize and a null width, caught (634 px in a 1,167 px card at 1,280); autosize without the null width, caught, since the height's own edit fixes the width again; a null width without autosize, not caught: plotly then draws the plot at the width it first drew it at until its own resize handler, 100 ms after the window's resize, redraws it at its box's width, and the test does not time those 100 ms.

---

## D-023 — Rytis's draft v3.10 merged: contour figures and cividis for continuous colour

**Date:** 2026-09-29 · **Status:** decided

**Context.** Rytis put draft v3.10 in `TEMP/` (`todo7.md`, `CLAUDE7.md`): contour figures in M2, M5, M6, M7 and M10 (two-knob maps, power maps, density contours, a likelihood concept figure) and cividis for continuous colour in both themes.

**Options.** Carry the draft as written, or item by item as in D-020 and D-021.

**Decision.** Item by item; what was carried and why is in TODO.md, "Merge of draft v3.10".
- Cividis is for sequential values, with isolines and labels in ink. The P2Q.6 heatmap keeps its diverging two-colour scale.
- The contour helpers and the M5 contour figures go in Phase 2D, after the record, because `src/ctx/de` is frozen until then.

**Evidence.**
- Cividis's ends, measured with the WCAG formula on matplotlib 3.11.2's map: #FEE838 is 1.25:1 against white and #00224E is 1.34:1 against black. A line or label in either end colour vanishes on that theme's page.
- The draft's claims checked on 2026-09-29, one verdict each in the merge record: cividis in both renderers and its source (PMID 30067751), the 1.92 and 3.00 likelihood drops (χ² quantiles halved: 1.9207, 2.9957), Schoenfeld (PMID 6354290) and the continuous-covariate form (Hsieh and Lavori, PMID 11146149), and pydeseq2 0.5.4's dispersion attributes, which answer the draft's one verify item.
- No frozen file changed; `~/ctx_tier1/data/raw` still held only `tier1_meta`.

**Consequence.** The style tokens gain cividis's stops when the helpers are built, so the measured-colour tests accept it. The drafts moved to `../RNAseq_tutorial_archive/TEMP_2026-09-29c/`.

---

## D-022 — The Trim Galore reports leave the cluster; the GENCODE 50 GTF is a resource; registered again

**Date:** 2026-09-29 · **Status:** decided

**Context.** The read-fates chart of Phase 2Q (P2Q.3) needs, per sample, the reads Trim Galore removed and why, and the rRNA share needs every gene's type. Neither was on the Mac: the reports stayed on the cluster, and the counts carry gene IDs without types. Both changes touch frozen files, `src/ctx/gather/rnaseq.py` and `conf/datasets/tier1.toml`, so they go in before cluster step 7 copies anything.

**Options.** Parse Trim Galore's numbers out of MultiQC's cluster report (`multiqc/`), which step 7 does not copy either and whose layout changes between MultiQC releases; or copy the per-sample reports (about 2.6 kB each) and parse them locally with MultiQC (Rytis's choice). For the gene types, a second annotation from another source, or the cluster's own GTF.

**Decision.**
- `COPIED_FILES` gains `trimgalore/*_trimming_report.txt`, and step 7's rsync the same `--include`.
- `test_repo::test_rsync_matches_copied_files` fails when the rsync and `COPIED_FILES` disagree in either direction, so the two lists of one thing stay one.
- `conf/datasets/tier1.toml` gains `annotation_gtf`, the file the cluster run used (`gencode.v50.primary_assembly.annotation.gtf.gz`), so gene IDs including `_PAR_Y` match the counts.
- The rules are registered again, at 2026-09-29 15:36.

**Evidence.**
- **The GTF,** downloaded on Rytis's approval on 2026-09-29: 124,650,284 bytes, MD5 289b91e5e95e8b0450d223246f10a12e, equal to the publisher's `release_50/MD5SUMS`; SHA-256 89bbad69a8c89fee5fadec0a2f14a098752d27a4ef7f24e4de9bac681e1b18f4.
- **The fixture,** rebuilt with `make_rnaseq_fixture.py --from-run` from the nf-core/rnaseq 3.27.0 test run of 2026-09-27 still in the scratchpad: the 14 files it held are byte-identical, 8 Trim Galore reports (32 kB) join them, and `tests/test_rnaseq_fixture.py` passes with the new pattern. The README now gives the run's date from its params file and the copy date apart.
- **No Step-1 count exists,** so registering again is allowed: `~/ctx_tier1/data/raw` held only `tier1_meta` at 15:36.

**Consequence.** Cluster step 7 waits for this commit. A sample whose report the rsync missed makes `ctx qc` fail on the Mac, not draw a chart without it.

---

## D-021 — Rytis's drafts v3.7–v3.9 merged: M10 clinical and survival, D27–D28, D3 and D6–D9; Phase 2Q planned

**Date:** 2026-09-29 · **Status:** decided

**Context.** Rytis put a third set of drafts in `TEMP/`: v3.7 (`todo4.md`, with `CLAUDE4.md`), which adds M10 on clinical covariates, treatment and survival, and v3.9 (`todo5.md`, byte-identical to `todo6.md`), which accepts the defaults of D3 and D6 to D9. He added `lessons5.md`, the lessons file of his planning sessions. The same day he approved the Phase 2Q plan: the read-fates bubble chart, pathway heatmaps and linked brushing.

**Options.** Copy the drafts over the repo files, which would undo D-020 and everything since; or merge item by item as in D-020.

**Decision.** Merge item by item. What was carried, what was not and why is in TODO.md, "Merge of drafts v3.7–v3.9".
- The drafts' D21 and D22 become D27 and D28; the drafts had numbered them after their own D20, which the repo had already used.
- M10 gets output number 12 and `ctx clinical`, so no output directory moves.
- The registration Phase 2Q needs is D-022, not D-021 as the approved plan numbered it: this merge came first.
- Phase 2Q brings forward, for Step-1, parts of M2 that D-020 put in Phase 3 (MultiQC's parse, the read-fates metrics, a flag-only status table, the `ctx qc` page) and the heatmaps TODO.md had in Phase 4. This amends D-020's "from Phase 3" for those items. (Added after the review of bb43cd5.)
- M10 fits Cox models with lifelines, which takes a pandas DataFrame, so CONVENTIONS.md's pandas boundary gains `ctx.clinical`. (Added after the review of bb43cd5.)

**Evidence.**
- **Claims checked at their sources** on 2026-09-29 (notes in the session scratchpad, `v39_facts/`), one verdict per claim in the merge record:
  - Verified: TCGA-CDR's endpoints and Table 3 recommendations, BRCA OS 151 events among 1,097 patients, the competing-risk tab, TCGA's rule against systemic therapy before procurement (42 neoadjuvant exceptions in the CDR, 102 GDC notifications; treatment for a prior malignancy may predate sampling), the CDR's ER+ vs ER− breast comparison, HPA's and bc-GenExMiner's best-cut-off splits, Lausen and Schumacher's corrected p, and GTEx's open donor fields.
  - Corrected: SCAN-B GSE96058 holds 3,273 patients with median follow-up 52 months, where the draft gave 3,217 and 75 months, which belong to a later QC subset (Brueffer et al. 2020) that GEO does not mark. GDC's neoadjuvant and prior-malignancy annotations are notifications, not exclusions. lifelines 0.30.3 does not install in Pyodide 314.0.7 (it pins pandas < 3; Pyodide ships 3.0.2), while 0.30.0 installs and runs.
  - Qualified: the survival times' unit in days rests on a secondary README; a distant-metastasis interval depends on a field TCGA did not require (648 patients who died with tumour have no new tumour event recorded).
- **Lessons.** Six of `lessons5.md`'s seven were already in `knowledge/lessons.md`; the seventh, counting open decisions from the file, was added, and the three-question rule gained a note that two check-ins today asked four.
- **Open decisions,** counted from the table after the merge: none. D1–D9 and D11–D28 are resolved and D10 is superseded.
- **No frozen file changed,** so the rules are not registered again here; `~/ctx_tier1/data/raw` still held only `tier1_meta` at 2026-09-29 15:27.

**Consequence.** M10 is planned for Phase 4, breast first, and its linked swimmer plot uses the linked-view controls built in P2Q.4. The drafts moved to `../RNAseq_tutorial_archive/TEMP_2026-09-29b/`.

---

## D-020 — Rytis's drafts v3.4 and v3.6 merged; the reuse layer follows the Step-1 record

**Date:** 2026-09-29 · **Status:** decided

**Context.** Rytis put two parallel drafts of the plan in `TEMP/`, draft v3.4 and draft v3.6, and asked for "more adjustability in selecting result visualisation, and an ability to easily repeat this whole analysis and generate these reports by providing different RNA-seq datasets". Both drafts branched from TODO.md before Phase 1b and Phase 2 were recorded.
- **Draft v3.4** adds:
  - M9 gene cards;
  - an M2 QC atlas;
  - the nf-core/rnaseq reference run.
- **Draft v3.6** adds:
  - cohort bundles with `--cohort`;
  - a Reproduce panel per figure, with `ctx figure`.
- **What both drafts restate.** They restate rules the repo has since changed, and they reuse the decision numbers D14 to D20 for different decisions.

**Options.**
- Copy a draft over TODO.md and CONVENTIONS.md. That would undo D5, D14, D16, D18, D-011, D-012, D-019, the 18 px base and the Step names.
- Merge item by item, keeping every later repo decision, with Rytis's choices where a draft changes direction.

**Decision.** Merge item by item. What was carried, what was not and why is recorded in TODO.md, "Merge of drafts v3.4 and v3.6". Rytis's choices on 2026-09-29:

| Question | Choice |
|---|---|
| M9 gene cards | Yes, Phase 4, output number 11. Step-1's worked examples are the `box_genes` registered before any count. |
| M2 QC atlas | Yes. The checkpoint runs in Track B too, from Phase 3. |
| Step-1 cohort QC | Sample-level checks only. MultiQC stays on the cluster (critique M3). |
| The drafts in `TEMP/` | Moved out of the repo. |
| When to build the reuse layer: cohort bundles, `ctx run`, `ctx figure`, Reproduce panels | After the Step-1 record, as Phase 2D. |
| The volcano's `lfc_type` | The x axis may show the shrunk estimate for display. Colours always follow the call on the raw estimate (D-012), and the page says so. This reverses the page control D-019 turned down, by Rytis's choice. On a shrunk axis the vertical fold-change lines are not drawn, because the call's threshold applies to the raw estimate; `x_estimate = "shrunk"` is refused, with the reason, when `--no-lfc-shrink` left the shrunk columns empty. |
| Figure parameter models | Frozen dataclasses, as `src/ctx/dataset.py` validates, not pydantic. |

- **Renumbering.** Draft decisions get new numbers:
  - D14 to D16 of v3.4 become D21 to D23;
  - D18 to D20 of v3.6 become D24 to D26;
  - D17 of v3.4 is already covered by D18 and D-006.
- **`ctx run`.** This is not in either draft. It is added as the direct answer to "easily repeat this whole analysis with a different dataset".

**Evidence.**
- **Diff and plan.** Six section agents, one fact-check pass, a planner and a critic.
  - The diff produced 181 items.
  - 89 claims were checked at their sources: 86 verified, 2 wrong and corrected, 1 unverified.
  - The critique found 21 problems in the first plan, all fixed before any edit.
  - The v3.6 claims: 21 of its new claims were checked at their sources: 16 verified, 3 wrong in detail and corrected, 2 partly verified and marked so.
- **Frozen rule sources.** Registered rules R1, R2 and R7 cited "TODO.md:360" and "TODO.md:364". The merge moved those lines, so the sources now name their target ("TODO.md, Verification: …") and cannot go stale again.
  - That changes the frozen `conf/datasets/tier1.toml`, so the rules are registered again.
  - At 2026-09-29 14:00, `~/ctx_tier1/data/raw` held only `tier1_meta`.
- **Line references inside TODO.md** were remapped to the lines they pointed at.
- **Tests.** `pixi run preflight` passes on the merged files.

**Consequence.**
- Phase 2D sits between the Step-1 record and Phase 3.
- Any draft item marked UNVERIFIED stays so until it is read at its source.
- Step 7 of the cluster block waits for this registration commit.

---

## D-019 — LFC shrinkage is a switch, on by default

**Date:** 2026-09-29 · **Status:** decided; its rejected page control is partly taken up by D-020 (the volcano's x axis may show the shrunk estimate for display)

**Context.** Rytis asked for an option to turn off LFC shrinkage, on by default, with an explanation of what it does. `ctx de` always ran apeGLM shrinkage. Its output feeds two places only:
- the `log2FoldChange_shrunk` and `lfcSE_shrunk` columns;
- the MA plot.

The test, the p-values, the call, the volcano and the enrichment use the raw estimate (D-012).

**Options.**
- A switch on `ctx de`, with the parameter card that every CLI option carries.
- A control on the page that swaps the estimate the figures show. That would put the shrunk fold change beside a call made on the raw one, which D-012 decided against.

**Decision.** `ctx de --lfc-shrink/--no-lfc-shrink`, on by default.
- **Defined once.** The default is `ctx.de.deseq.LFC_SHRINK`, and the CLI's default is tested against it.
- **When off,** `fit(shrink=False)` skips the second `DeseqStats` and its `lfc_shrink`, and leaves the two shrunk columns empty. The MA plot draws the raw estimate alone, and its caption, lead and reading say that shrinkage was off.
- **The card.** It says what shrinkage does: it pulls each fold change toward 0 as far as its counts leave it uncertain, so a weakly expressed gene cannot top a ranking on noise. It also says what turning it off changes.
- **How a switch's default prints.** A switch's card and its test now print the default as "on" or "off", where "True" would have appeared.
- **CONVENTIONS.md** says shrunken and raw LFC are both available. It now names the switch.

**Evidence.**
- **Nothing tested or called changes.** Synthetic `ctx de`, 6 contrasts, run with and without the switch:
  - all six results tables are equal in every column but the two shrunk ones, which are empty when off;
  - the metrics, contrasts and sex-check tables are identical.
- **Time.** 32.71 s wall with shrinkage and 31.62 s without, both at 0.58 GB peak, on 1,000 genes. The card promises no speed-up.
- **Tests.**
  - `test_de.py::test_no_shrink_changes_nothing_tested` compares the two fits column by column.
  - `test_ma_draws_the_estimates_the_fit_has` checks the MA plot's series and its text.
  - `test_cli.py::test_lfc_shrink_reaches_the_fit` checks that the flag reaches `run` in both spellings and by default.
  - `test_lfc_shrink_default_single_source` checks the default against `LFC_SHRINK`.
  - Three mutants each fail one of these tests: shrinking whatever the switch says, the CLI dropping the flag, and the MA plot always drawing both series.
- **In the browser.** The shrinkage-off page shows the card with default "on", its effect text and `--no-lfc-shrink` in the reproduce command. The MA page draws one series, "raw estimate", for 997 genes. The default page still draws both series.
- **Registered again.** `src/ctx/de` is frozen. At 2026-09-29 10:20, `~/ctx_tier1/data/raw` held only `tier1_meta`.

**Consequence.** A table written with the switch off has empty shrunk columns. A reader ranking genes needs shrinkage on, and the page's "In practice" advice still says to rank by the shrunk estimate.

---

## D-018 — The M5 page says what BH promises and when it does not hold

**Date:** 2026-09-29 · **Status:** decided

**Context.** The M5 glossary said BH and BY "control the false discovery rate", with no condition attached. D-017 measured that DESeq2's Wald p-values run small because each gene's estimated dispersion is taken as known. On synthetic data, the share of false genes at padj < 0.05 is therefore about twice what BH promises. Rytis asked for the small-sample caveat on the page.

**Options.**
- A caveat that uses this run's own numbers. No run knows its false genes, so there are none to fill in.
- A fixed, measured example, stated with what it was measured on.

**Decision.** The measured example, on the M5 page in two places:
- **Glossary.** The BH and BY line now says they control the false discovery rate "when the p-values are exact".
- **Common mistakes.** A new entry covers reading padj < 0.05 as at most 5 % false genes. BH bounds the average share, and only when the p-values of unchanged genes are exact. The entry gives the measured numbers below.

**Evidence.** Measured on the committed `tests/synth.py` `design_counts` defaults: 6 vs 6, 1,000 genes, 100 of them changed, seeds 0 to 199, one fit each.

| Share of false genes | Mean | Standard error | Max |
|---|---|---|---|
| BH's promise, 0.05 × m0/m | 0.045 | | |
| padj < 0.05, any fold change | 0.084 | 0.0019 | 0.149 |
| the page's call, which also needs abs(log2 FC) ≥ 1 | 0.032 | 0.0012 | 0.084 |

D-017 gives the cause. With the true dispersion put in place, the null p-values are nominal, and R DESeq2 1.46.0 gives the same excess as ctx.

**Consequence.**
- **Frozen file.** `src/ctx/de/page.py` is frozen, so this commit registers again. At 2026-09-29 10:10, `~/ctx_tier1/data/raw` held only `tier1_meta`.
- **Fixed numbers.** The numbers in the mistake are a fixed reference, not this run's. If `design_counts`, the call rule or pydeseq2 changes, they have to be measured again.

---

## D-017 — Synthetic dispersion follows a trend; the planted tests rest on measured spreads

**Date:** 2026-09-28 · **Status:** decided

**Context.** Preflight printed four kinds of warning. Two came from pydeseq2:
- "The dispersion trend curve fitting did not converge" on seven fits;
- "overflow encountered in exp" on one.

The first had a root cause in `tests/synth.py`: `design_counts` gave every gene the same dispersion, 0.1. With no trend to fit, pydeseq2's parametric fit fails and falls back to the mean. With zero true spread, the dispersion prior's variance also sits at its floor of 0.25. Giving the data a trend broke three tests at seed 7, and each of them had been passing by the draw.

**Options.**
- Keep the flat data and assert pydeseq2's warning where it appears. That is honest about the data, but it fixes nothing.
- Give the data a trend and re-measure every planted assertion over many seeds.

**Decision.** A trend, with each assertion set from measured spreads.

- **Generator.** `design_counts` draws dispersion as (0.05 + 10 / mean) × exp(N(0, 0.3)). That is setting C7, the best of 8 settings swept.
- **Filtered fixture.** `filtered_counts` adds 1,500 low null genes at mean 20, dispersed the same way, where it used to add 600 flat ones.
- **Planted FDP bound.** In `test_planted_recovered`, the FDP bound on the padj < 0.05 set is 0.20; the page's call keeps 0.10.
- **Null split.** The CLI test holds the null split to rule S1, at most 1 % of tested genes, where it used to require zero.
- **Significant against called.** The CLI assertion that significant differs from up plus down needed a chance false positive. `test_de.py::test_significant_counts_any_fold_change` now checks the same thing on a constructed gene.

**Evidence.**
- **Why BH's 0.05 is not the bound.**
  - Setup: all genes null, 6 vs 6, flat data at donor SD 0, seeds 0 to 24, 74,884 genes.
  - ctx's p-values fall below 0.05, 0.01 and 0.001 at rates of 0.0595, 0.0148 and 0.0021.
  - With the true dispersion put in place they are 0.0516, 0.0102 and 0.0012.
  - Plain pydeseq2 gives ctx's rates to four decimals. So does R DESeq2 1.46.0 on all 9,000 genes of the same matrices, seeds 0 to 2.
  - So the excess comes from the method, which takes an estimated dispersion as known. ctx does not cause it.
  - Two skeptics could not refute this. One reproduced it on seeds 10 to 14, and the ctx arm and the true-dispersion arm disagree beyond chance (McNemar p = 1.2e-4 at 0.001).
- **The planted FDP, under C7, 6 vs 6, 1,000 genes, 100 planted:**

| Set | Seeds | Mean | Max | Bound |
|---|---|---|---|---|
| padj < 0.05, any fold change | 0 to 199 | 0.083 to 0.084 | 0.149 | 0.20 |
| the page's call, with abs(log2 FC) ≥ 1 | 0 to 199 | not recorded | 0.084 | 0.10 |
| raw p < 0.05, no correction | 0 to 199 | 0.353, on seeds 0 to 99 | minimum 0.273 | caught by 0.20 |

- **Trend-warning fits per setting.** The eight settings, summed over the seven test designs:

| Setting | a0 | a1 | scatter | Fits with the warning |
|---|---|---|---|---|
| C0 | 0.1 | 0 | 0 | 313 of 480 |
| C7 | 0.05 | 10 | 0.3 | 22 of 480, all in the old 600-gene fixture; 0 with the new one |

  On fresh seeds 100 to 119, C7 raised 0 such warnings in every `test_de` and `test_stats` fit. It raised 1 in 160 fits of the CLI's contrasts: b_case, seed 107.
- **The new fixture, measured on the committed code.** Over seeds 0 to 29 and 100 to 119:
  - the reasons hold on 50 of 50 seeds;
  - independent filtering removes at least 141 genes;
  - no fit warns.
- **Seed 7, the seed the tests use.** The warning-free result is from a pytest run of the test files.

| Quantity | Value |
|---|---|
| padj-set FDP | 0.107 |
| called-set FDP | 0.048 |
| raw-only genes | 13 |
| genes removed by independent filtering | 748 |
| paired vs unpaired, true genes significant | 80 vs 0 |
| paired vs unpaired, true genes called | 46 vs 0 |
| pydeseq2 warnings in these fits | none |

- **Preflight and budget.**
  - Preflight: 344 passed, 4 skipped, 0 warnings, 148 s wall, 1.63 GB peak.
  - The 60,000-gene budget test passes: `ctx de` in 296 s with an 11.0 MB page, `ctx function` in 15.3 s.
- **The metric test discriminates.** Counting the significant metric with the fold-change filter fails the new unit test.
- **The overflow is upstream.** It is scverse/PyDESeq2 issue 169, open on 0.5.4, the newest release, at `utils.py:1088` inside `lfc_shrink`. On C7 it appears in 12 of 30 unpaired fits at donor SD 1.0, never at seed 7. It is not suppressed. A fix that rewrites `1 / (1 + size·exp(-η))` as `expit(η - log size)` is drafted in the session scratchpad, not posted. It leaves raw log fold changes and p-values bit-identical, and moves shrunken ones by at most 5.6e-12.

**Consequence.**
- **Synthetic numbers change.** Every synthetic count matrix from `design_counts` differs, so the example pages drawn from it change too.
- **The FDP bound is tied to this generator.** 0.20 was measured on C7 alone. Pooled over all 8 settings, the FDP reached 0.198.
- **Not done: the student page.** The page states that BH controls the FDR at 0.05, with no small-sample caveat. It is unchanged.
- **Not done: the layout-check crash.** On 3 of 20 fresh seeds, `ctx de` stopped at the box figure's layout check, on overlapping p and r labels. That is a separate fault, and it is not fixed here.

---

## D-016 — Registered again for the fixes an adversarial review of D-015 found

**Date:** 2026-09-28 · **Status:** decided

**Context.** An adversarial review of commit f95010b (D-015) gave 10 findings. One skeptic per finding tried to refute it; 8 survived, and one of those was reported twice. None changes a shipped page. The fixes touch `src/ctx/de/stats.py` and `conf/datasets/tier1.toml`, which are frozen (D-013).

**Options.**
- Leave the frozen files as they are and correct only the tests and documents.
- Register again, which is allowed while no Step-1 count exists.

**Decision.** Register again, and fix each finding:

| Finding | Fix |
|---|---|
| `p_line`'s "if and only if" breaks at the last bit | The docstring and D-015 say where it holds, and a `ponytail:` comment names the clamp that would close it |
| `tests/test_viz.py`'s live-volcano fixture drew the old line, with n_tested 7 against 8 tested rows | It computes its line with `p_line` over its tested rows; the test writes the bar out as 5 × 0.05 / 8 |
| Nothing checked that the static line equals the line the page script draws on load | `tests/test_page_script.py::test_static_line_is_the_line_the_page_draws_on_load` runs `run._volcano` on a fit with filtered genes through `figure_block` and the page script |
| A viz test still expected no line when nothing passes | It now expects the first bar; a second test keeps the contract's own no-line state |
| The registration comment, TODO.md and D-012 still described the old line or `p_cutoff` | Updated in place, with pointers to D-015 |
| D-015 overstated its test coverage | Reworded, and the tightness test now runs at 0.01 as well |

**Evidence.**
- **The parity test discriminates.** A producer that counts m over every gene, not the tested ones, passes on the default design, where every gene is tested. On `tests/synth.py::filtered_counts`, where Cook's and independent filtering remove genes, it fails. That helper now builds this matrix for `test_de`, `test_stats` and the parity test, in place of two inline copies.
- **No count exists.** At 2026-09-28 21:45, `~/ctx_tier1/data/raw` held only `tier1_meta`.
- **Refuted findings, not acted on.** Two findings were refuted:
  - `ScatterData.all_lines` raises on an empty frame. That happens before f95010b too, so this commit did not cause it.
  - `tests/synth.py`'s no-line sentence can show only when no gene was tested.

**Consequence.** Step 7 of the cluster block waits for this registration commit.

---

## D-015 — The volcano's cutoff line sits at the bar the weakest passing gene cleared

**Date:** 2026-09-28 · **Status:** decided

**Context.** Rytis reported that the volcano's p-value cutoff "doesn't seem to work properly": only some cutoffs moved the dashed line. The line sat at the largest raw p among passing genes (D-012). That value moves only when a gene crosses the cutoff, and no gene's BH-adjusted p lay between several of the preset cutoffs, so the line stood still between them.

**Options.**
- Keep the line at the weakest passing gene and explain on the page why it does not move.
- Draw the line at the bar the procedure itself sets, which moves with every cutoff and still separates passing from failing genes.

**Decision.** The procedure's own bar, defined once in `ctx.de.stats.p_line` and repeated in `ctxCalls.line` in `templates/_calls.js.j2`. With m tested genes, k of them passing and j = max(k, 1):

| Basis | Line on the raw-p scale |
|---|---|
| raw p | α |
| BH | j·α/m |
| BY | j·α/(m·c(m)), with c(m) = 1 + 1/2 + ... + 1/m |
| Holm | α/(m − j + 1) |
| Bonferroni | α/m |

A tested gene passes if and only if its raw p is below the line, in exact arithmetic. In floating point, a p within one unit in the last place of the bar can land on either side. The review behind D-016 built such cases from p-values on the bar's own double, and found 0 of them in 1,600 settings with continuous p-values. When no gene passes, the line is the first bar, the one the strongest gene would have had to clear. A line is missing only when no gene was tested.

**Evidence.**
- **Before and after, on the same page.** Synthetic M5 volcano, 999 tested genes, BH, measured in the browser pane:

| Cutoff | Genes passing | Old line, -log10 p | New line, -log10 p |
|---|---|---|---|
| 0.1 | 106 | 2.00 | 1.97 |
| 0.05 | 103 | 2.32 | 2.29 |
| 0.02 | 100 | 4.85 | 2.70 |
| 0.01 | 100 | 4.85 | 3.00 |
| 0.005 | 100 | 4.85 | 3.30 |
| 0.001 | 100 | 4.85 | 4.00 |

- **The line still separates.** At 30 settings, 5 bases times 6 cutoffs, 0 tested genes on the page sit on the wrong side of the drawn line.
- **Tested in `tests/test_stats.py`.**
  - At every basis and four cutoffs, the line splits passing from failing genes.
  - At every basis and four cutoffs, the weakest passing gene still passes when moved to 1 − 1e-9 of the line, and fails at 1 + 1e-9. The first commit tested three cutoffs here; D-016 added 0.01.
  - The reported case, BH only: cutoffs 0.05 and 0.02 call the same genes and give lines 2.5 times apart.
- **The page script is tested against Python.** `tests/test_page_script.py` checks the page script's line against `p_line` over 80 settings. The BY constant is summed pairwise by numpy and in order by the page, so the two agree to a relative 1e-12.
- **Five mutants fail the tests:**
  - k in place of max(k, 1);
  - (k + 1)·α/m for BH;
  - α/(m − k) for Holm;
  - BY without c(m);
  - the page's BH line times 1.01.
- **Registered again.** `src/ctx/de` is frozen (D-013), so this commit changes `[sanity] registered`. At 2026-09-28 21:24, `~/ctx_tier1/data/raw` held only `tier1_meta`. A search under `~/ctx_tier1` for counts, salmon, rnaseq or quant files found nothing.

**Consequence.** The line's height and the reading sentence "above the dashed line at -log10 p = X" change for every volcano, in the PDF and on the page. The calls, counts and tables do not change. The no-line sentence now covers only the case where no gene was tested.

---

## D-014 — The cohort tiers are called Steps; pages show a dataset's label

**Date:** 2026-09-28 · **Status:** decided

**Context.** Rytis asked for "tier1" and the other tiers to read "Step-1", "Step-2" and "Step-3". A dataset's slug names every output file, and file names are lowercase with underscores, so "Step-1" cannot be a slug.

**Options.**
- Change only the shown name and the prose.
- Rename the slug and every path to `step1`: the dataset file, output names, `~/ctx_tier1` and the cluster paths.
- Change the prose only.

**Decision.** The shown name and the prose; this is Rytis's choice of 2026-09-28.
- `[dataset] label = "Step-1"` in `conf/datasets/tier1.toml`, read through `Dataset.name`, is what page titles and the hub show. Without a label a page shows the slug.
- The documents, test comments and page text say Step-1, Step-2 and Step-3 wherever they meant a cohort tier: 99 replacements. "data tier", which means access level, is untouched.
- The slug `tier1`, every file name, `conf/datasets/tier1.toml`, `~/ctx_tier1` and the cluster paths keep `tier1`, so nothing already built or running is renamed under it.

**Evidence.**
- **The freeze, handled.** The dataset file is frozen (D-013), so the rules were registered again in the same commit that adds the label. `[sanity] registered` carries a time, so each registration changes that line.
- **No count exists.** On 2026-09-28 at the registration, `~/ctx_tier1/data/raw` held only `tier1_meta`, and a search for counts, salmon or rnaseq under `~/ctx_tier1` found nothing.
- **The label is tested.** `tests/test_curation.py::test_label_is_the_name_pages_show`.

**Consequence.** A future dataset can carry any label. Its slug stays the file-name token, and the two can differ.

---

## D-013 — The Step-1 sanity rules are registered again, before any count exists

**Date:** 2026-09-28 · **Status:** decided

**Context.** The rules were registered on 2026-09-27 in commit 8219d16. `tests/test_sanity_record.py::test_rule_unchanged` then froze every file a verdict depends on, until the record commit:
- `src/ctx/de` and `src/ctx/function`;
- `src/ctx/viz/annotate.py` and `src/ctx/report/hub.py`;
- the dataset files.

The page redesign had to change those files: the call rule (D-012), the page text and the PERMANOVA without a movable label. Two of the rules, R3 and R7, also moved to the `significant` metric.

**Options.**
- **Two phases.** Change the pages now, and M5, M6 and the hub only after the record, when the cluster run is done.
- **Narrow the freeze list,** which would still mean one edit to frozen files.
- **Register again,** before any count exists.

**Decision.** Register again. This is Rytis's choice of 2026-09-28.
- The `[sanity]` table of `conf/datasets/tier1.toml` now holds `registered = "2026-09-28"`.
- `rule_commit()` returns the newest commit that changed that line. The record must still come after it, and the freeze runs from it.
- R1 to R8 keep their thresholds. R3 and R7 read `calls.<c>.significant`, BH padj < 0.05 at any fold change, which is the meaning they had, now under an explicit name.

**Evidence that no count existed at registration.** On 2026-09-28, before this commit:
- `~/ctx_tier1/data/raw` held only `tier1_meta`, the metadata-only fetchngs run.
- `find ~/ctx_tier1 -iname '*counts*'` and `-iname '*salmon*'` found nothing.
- No cluster step had returned output to this Mac. The cluster block's step 7, the only step that copies counts here, waits for this commit.

**Consequence.**
- **The freeze.** Any change to a frozen path between this commit and the record fails `test_rule_unchanged` when the record is committed.
- **Registering once more** is allowed only while no count exists, and needs the same evidence line.
- **An early warning** was first left out, because it blocks every M5 or M6 fix before P2C.3 until the rules are registered again. On Rytis's request of 2026-09-28 it was added: `tests/test_sanity_record.py::test_frozen_since_registration`. Until the record exists, it fails as soon as a frozen path differs from the registration commit, committed, uncommitted or untracked, and names the two ways out.

---

## D-012 — A gene is called at padj < 0.05 and |log2 FC| at least 1; four adjustments on the page

**Date:** 2026-09-28 · **Status:** decided

**Context.** Rytis asked every figure to say in practice what it shows, with the thresholds, and the volcano to take the reader's own cutoffs. Until then a call was BH padj < 0.05 signed by the raw fold change, with no fold-change threshold, so the volcano had no x cutoff to name.

**Options.**
- Keep padj alone and say that the x cutoff is 0.
- Filter the fold change after the Wald test.
- Test against the threshold, as DESeq2's `lfcThreshold` does, which needs a Wald test per threshold.
- Guide lines at plus or minus 1 without changing the call.

**Decision.** These are Rytis's choices of 2026-09-28:
- **Call rule.** Called means BH padj < 0.05 and |raw log2 FC| at least 1, the fold change filtered after the test. It is `ctx.de.stats.classify`, the one definition.
- **Volcano controls.** The page volcano lets the reader choose:
  - the fold-change threshold, in log2 or as a percent;
  - the p basis: the BH, BY, Holm or Bonferroni adjusted p, or the raw p;
  - the cutoff, from 0.0001 to 0.1.
- **Adjustments.** They come from the same Wald p-values over the genes pydeseq2's BH kept, so no refit is needed. `deseq.py` stores `padj_by`, `padj_holm` and `padj_bonferroni`.
- **Sanity rules.** R3 and R7 test the model's false positives and power, which a fold-change filter would flatter. They read `calls.<c>.significant`, padj < 0.05 at any fold change.

**Evidence.**
- **BH.** Our BH equals pydeseq2's `padj` exactly: the maximum difference is 0.0 on a fit with 1,071 tested genes, 527 independently filtered and 1 Cook's outlier. BH over all 1,598 p-values differs by up to 0.0933, so the tested set matters.
- **Other methods.** All four methods equal statsmodels `multipletests` to 2.2e-16.
- **Planted truth, old rule against new:**

  | Measure | Old rule | New rule |
  |---|---|---|
  | planted calls / false calls | 107 / 7 | 101 / 1 |
  | naive batch share | 0.980 | 0.890 |
  | adjusted batch share | 0.010 | 0.000 |

  All planted-truth tests pass unchanged.
- **Mutation checks.** Four mutations of the rule and the adjustments each failed a test (`tests/test_stats.py`, 47 tests).
- **Page.** On the synthetic M5 page, the live counts and the cutoff-line height equal a Python recount from the results parquet at 8 settings: all four methods, raw p, log2 and percent thresholds, and a fold-change threshold of 10 that no gene reaches. The branch for a cutoff no gene passes cannot be reached on this fixture, whose minimum cutoff of 0.0001 still passes genes, so it ran on a scratch copy of the page with every embedded p-value set to 0.5: the line hid, 0 up and 0 down, and the reading switched to its no-line sentence. Superseded by D-015: with no gene passing, the line now sits at the first bar and stays drawn.

**Consequence.**
- **Changed counts.** Every "up" and "down" count, the ORA input and the call-count bars now need the fold change, so they differ from pages written before 2026-09-28.
- **The height is not the call.** The volcano's height is the raw p-value. Its horizontal line sits at the largest raw p whose adjusted p passes, which is exact because each adjustment is monotone in p. The page says so. Superseded by D-015: the line now sits at the bar the weakest passing gene cleared, which moves with the cutoff.
- **Percent thresholds.** A percent is asymmetric in log2: +20 % is log2 1.2 = 0.263, and the same threshold down is 16.7 %. The page shows both.
- **Reversing it.** Going back would change the rule, the metrics names and the sanity rules. It would need a new registration (D-013).

---

## D-011 — The statistics caption sits under each interactive figure, not inside it

**Date:** 2026-09-28 · **Status:** decided

**Context.** Each interactive figure drew its caption as a plotly annotation and the page printed it again below. On the MA plot the annotation overlapped the x-axis title, measured in the browser on 2026-09-28, and the fixed-width annotation does not reflow on a phone.

**Options.**
- Keep the caption in the plot with its spacing fixed, so a PNG export keeps the statistics.
- Print it once, below the plot, as page text.

**Decision.** Printed once, below the plot. This is Rytis's choice of 2026-09-28.
- The PDFs and SVGs keep their caption in the figure.
- The plotly toolbar, and with it the PNG export, is off. That export would now carry no statistics, and over a phone-width plot the toolbar covered the legend.
- CONVENTIONS.md (Figures, the Plotly rule) now says where the caption goes.

**Evidence.** `test_round_trip` checks the caption tokens:
- in the page's caption, the PDF and the SVG;
- in none of the plot's annotations;
- the box figure's per-gene p and r labels stay drawn in the plot.

**Consequence.** Anyone who copies a plotly figure out of a page takes no statistics with it. The page's caption, the standalone figure page, the PDF and the SVG carry them; the CSV's README describes the columns, not the test.

---

## D-006 — Raw reads enter only through nf-core/fetchngs 1.12.0; `ctx gather` is the only ctx reader of its layout (TODO.md D18)

**Date:** 2026-09-26 · **Status:** decided

**Context.** `knowledge/recommendation1` (Rytis, 2026-09-26; moved out of the repo on 2026-09-27 to `/Users/rytis/Scripts/Pipeline/OLD/RNAseq_tutorial_recommendation1_2026-09-26/`) replaced custom ENA/`fastq-dl` download logic with nf-core/fetchngs 1.12.0 and split M1 into find, fetch and make usable. Its layout description followed the fetchngs docs; parts of it do not match the 1.12.0 code.

**Options.** Custom download code in ctx; fetchngs 1.12.0, which the lab runs on Kebnekaise; fetchngs 1.13.0 (released 2026-09-15, needs Nextflow >= 25.10.4, renames the versions file, changes the samplesheet row order, untested in the lab).

**Decision.** fetchngs 1.12.0 with `--nf_core_pipeline rnaseq`, one fresh outdir per run under `data/raw/<slug>/`, read-only to ctx; `ctx gather` (package `ctx/gather/`) is the only ctx module that reads it; no launcher and no ENA/SRA access code in ctx. The reader stays in ctx under D-003.

**Evidence.** From the 1.12.0 code (raw.githubusercontent.com/nf-core/fetchngs/1.12.0, retrieved 2026-09-26):
- File names `<EXP>_<RUN>[_1|_2].fastq.gz` (`bin/sra_runinfo_to_ftp.py:143-145`); `.md5` files only from FTP and Aspera, written from ENA's md5 and checked with `md5sum -c` (`modules/local/sra_fastq_ftp/main.nf:29,44,52`); the sra-tools fallback is automatic per run and writes no md5 (`workflows/sra/main.nf:79-81`).
- `*.runinfo.tsv` is not published (`modules/local/sra_ids_to_runinfo/nextflow.config:3-6`); `metadata/<INPUT_ID>.runinfo_ftp.tsv` is, one per input line.
- `--nf_core_pipeline rnaseq` only inserts `strandedness` = `auto` (`modules/local/sra_to_samplesheet/main.nf:33-49`); sheet paths are `<outdir as given>/fastq/<name>` (same file, `:24-53`).
- The input regex rejects `E-MTAB-*` (`subworkflows/local/utils_nfcore_fetchngs_pipeline/main.nf:144-165`); an SRR pulls every run of its experiment.
- Minimum Nextflow `!>=23.04.0` (`nextflow.config:223`); container images are amd64-only, so Apple silicon needs `-profile docker,arm` (`nextflow.config:110-112`).
- On Kebnekaise the lab runs 1.12.0 with `-profile singularity` on Nextflow 26.04.4, which needs `NXF_SYNTAX_PARSER=v1` (`General_tools/cluster_utils/SHARED_run_fetchngs.sbatch:56-77`; `cluster_utils/README.md`).

Observed while building the Phase 1b fixtures on this Mac (2026-09-26, `-profile docker,arm`, Nextflow 26.04.6):
- `metadata/*.runinfo_ftp.tsv` rows come in no fixed order: a rebuild swapped two runs of GSM4907283. `samplesheet.csv` came back byte-identical. ctx reads metadata as a set of runs; the rebuild test compares those files by rows.
- The sra-tools fasterq-dump step asks for 36 GB by default and fails on a 24 GB laptop; the fixture script caps the run at `--max_cpus 4 --max_memory 8.GB`.
- In one rebuild, `prefetch` inside the amd64-emulated sra-tools container failed to open the run from NCBI's S3 bucket and then sat idle, while the host reached the same URL in 0.4 s. The FTP runs were unaffected. Class runs should prefer runs with an ENA FASTQ link.

**Consequence.** An upgrade to 1.13.0 is one deliberate change: `FETCHNGS_VERSION`, a new fixture, and the renamed versions file. `ctx gather` refuses any other release until then.

**Addendum, 2026-09-27: two Nextflow warnings, the logs flag, and link publishing.**

Two warnings appeared on fetchngs runs from this Mac's default `~/.nextflow`; runs with a fresh `NXF_HOME`, such as the fixture builder's, showed only the second. Neither changed an output: ten metadata-only runs of one run accession (SRR13191702), with and without the fixes below, wrote byte-identical `samplesheet/` and `metadata/`; a one-run input has no row order to vary.

- **Multi-revision warning.** Its cause is the local asset cache, not ctx. `~/.nextflow/assets/nf-core/fetchngs` was a legacy single-revision clone; the clone dates point to Nextflow 25.10.4, which is inferred, not recorded. Nextflow 26.04 warns whenever it runs from a legacy clone (`CmdRun.groovy:640-644`). On Rytis's instruction, `nextflow pull nf-core/fetchngs -r 1.12.0 -migrate` was run on this Mac. It exited 0 in 9.87 s and added `~/.nextflow/assets/.repos/nf-core/fetchngs/{bare,clones/8ec2d934f9...}`: 7,312,049 bytes in 164 files, 7,588 KiB on disk. A `du -sh` right after the migrate read 8.4M; a review measured 7.4M with no file changed since, so the 8.4 first written here was wrong. The legacy clone stays in place, byte-identical across its 163 files, so deleting `.repos` reverts. The `nf-core/riboseq` legacy clone belongs to other projects and was not touched. A fresh cache under Nextflow 26 is multi-revision already, so students do not see this warning.
- **`monochromeLogs` warning.** Its cause is upstream. fetchngs 1.12.0 pins nf-validation 1.1.3, whose `paramsSummaryLog` reads the undefined `params.monochromeLogs` whenever `monochrome_logs` is false (`SchemaValidator.groovy:880-886`, tag 1.1.3; reported upstream as nf-schema#14). Nextflow 25.10.4 warns on the same run. fetchngs 1.13.0 moves to nf-schema and should not warn, but it was not run. Rytis chose `--monochrome_logs true` over a monochromeLogs-only parameter and over the 1.13.0 upgrade. It is a real fetchngs parameter (`nextflow_schema.json:210`); it removes both the warning and the colour codes from the log; and it sets `monochrome_logs` to true in `params_*.json`, which always carries the key (false by default), so that file's SHA-256 in gather's provenance changes. `tests/test_repo.py::test_fetchngs_runs_monochrome` is a guard, not a proof. It checks the M1b course command, the CONVENTIONS rule, the TODO and CONVENTIONS lines that route a run through the Kebnekaise launcher with `-profile singularity`, every fenced TODO block that calls the launcher (its `printf` write of the config), and every file under `scripts/` whose code names `nf-core/fetchngs`. The test's comments name what escapes it.
- **The Kebnekaise launcher takes no extra flags**, so the cluster run carries the setting in `<project_dir>/nfcore_fetchngs/nextflow.config`, the directory the launcher `cd`s into before `nextflow run` (`SHARED_run_fetchngs.sbatch:110,201`). Tested on this Mac the same day, with the cluster's config file and fetchngs started without the flag, as the launcher starts it:
  - `monochrome_logs = true` there reaches the run: the params file records it, and the log has 0 ANSI escapes and no monochromeLogs warning. It had no multi-revision warning either, but only because the cache had been migrated about three minutes earlier (the migrate ended 19:42:07, the run started 19:45:15); the config has no bearing on that warning.
  - `publish_dir_mode = "link"` there does NOT reach `publishDir`, which the module configs resolve while fetchngs's own config is parsed. `nextflow config` from that directory prints `publish_dir_mode = 'link'` among the params, yet `mode:'copy'` for `SRA_FASTQ_FTP`, and the published `metadata/*.runinfo_ftp.tsv` had 1 link.
  - `--publish_dir_mode link` on the command line does hard-link: the same file had 2 links.
  - So the cluster FASTQs are copied, not linked, and the Phase 2 cluster block now budgets for two copies (critique M1 settled). One copy would need the launcher to pass the flag, which is Rytis's call.
- **2026-09-28: the launcher passes the flag.** On Rytis's instruction, `General_tools/cluster_utils/SHARED_run_fetchngs.sbatch` gained a `PUBLISH_DIR_MODE` setting, passed as `--publish_dir_mode`. It is not committed in General_tools. Empty keeps fetchngs's copy, and a value outside the schema's enum (`nextflow_schema.json:193`) fails the preflight. Both arms ran through the launcher on this Mac for ERR1160846, one FTP run of 18,077 bytes. Two local stand-ins replaced cluster settings: a launch-directory config set the local executor, since the launcher exports `NXF_EXECUTOR=slurm`, and laptop caps.
  - With `link`, the FASTQ and its md5 each had 2 links and shared the inode of the file in `work/`. Without the setting, each had 1 link and a different inode.
  - Every file a process publishes was linked, 4 of 13. The other 9 are written by the launcher, by Nextflow's `pipeline_info`, or by `collectFile(storeDir:)` (`workflows/sra/main.nf:154,163,182`), and none is a FASTQ.
  - The cluster block's step 1 now sets `PUBLISH_DIR_MODE=link` and budgets for one copy. Whether the cluster filesystem hard-links is measured in step 3, not assumed.
- **The committed fixtures predate the flag** and stay as built; their READMEs record the command that made them. With the flag, a metadata_only rebuild matched the committed fixture in all 13 samplesheet and metadata files, and in its versions file. The timestamped params file and trace differed, the params file only in `monochrome_logs`. So did README.md, whose command row now carries the flag and a new date, and MANIFEST.sha256, which hashes those files. Two download rebuilds, the second a retry, stalled in `SRATOOLS_PREFETCH` for 10 min 22 s and 9 min 14 s, from the first prefetch submit in the execution trace to the SIGTERM that stopped them. Each prefetch task wrote one informational line and nothing more, with an empty `.command.trace`. That is the symptom of the stall recorded above; its cause was not logged this time. (A first version of this sentence said 10 min 23 s and 9 min 15 s, counted from the NCBI-settings task.)


**Correction 2026-09-29.** The Options paragraph says fetchngs 1.13.0 "changes the samplesheet row order". It changes the sort key: 1.12.0 sorts the per-sample files by their `<EXP>_<RUN>` name (`workflows/sra/main.nf:152`, `sort: { it.baseName }`), 1.13.0 sorts on the row text (`collectFile(sort: true)`, PR #379). The two orders differ only when one accession is a string prefix of another; on both committed fixtures they are the same (simulated; 1.13.0 was not run). The consequence stands: `ctx gather` reads by column name and compares runs as sets, now a rule in CONVENTIONS.md (D-020).


## D-005 — JupyterLite DE: precomputed grids plus live scipy tests, not live pydeseq2 (TODO.md D17)

**Date:** 2026-09-26 · **Status:** decided

**Context.** D4 left open whether the M5 notebook re-runs pydeseq2 in the browser. P1.5 fixed the
rule as `decide()` in `scripts/python/spike_pyodide.py`, committed in 247d3f5 before anything ran.
The record is `deliverables/2026-09-26_lite_build_pyodide_spike_v1.csv`.

**Options.** Live pydeseq2 0.5.4 in JupyterLite; the D4 fallback, precomputed grids plus live
scipy tests; live pydeseq2 under an in-process joblib backend, which the rule does not cover.

**Decision.** The fallback at all three grid sizes, 1,000 x 24, 20,000 x 24 and 20,000 x 200, with
plotly for figures. Stage B, the browser timing, did not run, because no size passed stage A.

The rule's "the `~ cell_type + condition` fit" was read as pydeseq2's default call. That reading
was fixed in the harness commit 4a54ff0, after a scratch dry run had shown the failure, not in the
rule commit 247d3f5. Under pydeseq2's public backend argument (the shim below) the 1,000 x 24 fit
meets every stage-A criterion, so the verdict at that size rests on the reading; how it was
settled is TODO.md D17. (The "D-005a" in AGENTS.md is the template's numbering, not this entry.)

**Evidence.** Stage A ran npm `pyodide@314.0.6`, Python 3.14.2, under Node v24.15.0 on macOS
15.7.9 arm64. npm integrity:
`sha512-BKDTJyIqFxC4BExLqeRS3f5xvXZIjOt8C3zGLN/Cc7tFxSwvKVhVkchQQ2AGLOtR4YrVOIFxbV8poyDOOmWwxQ==`.
- micropip installed pydeseq2 0.5.4 as 47 files of 41,800,855 bytes. anndata 0.13.4 arrived with
  typing-extensions 4.16.0 from PyPI, because the bundled 4.15.0 was never loaded, so the known
  risk did not fire in Node.
- The fit failed at every size. First blocking error, verbatim:
  `AssertionError: ThreadingBackend does not accept setting the inner_max_num_threads argument.`
  pydeseq2 0.5.4, the newest release on PyPI (https://pypi.org/pypi/pydeseq2/json, uploaded
  2026-01-23, retrieved 2026-09-26), wraps every inference call in
  `parallel_backend(self._backend, inner_max_num_threads=1)` (`default_inference.py:96`), where
  `backend` is a public `DefaultInference` argument that defaults to `"loky"`. Pyodide has no
  loky, joblib 1.5.3 falls back to threading, and threading rejects that argument. No version pin
  reaches it; the backend argument does.
- Loaded versions, Pyodide then CPython: numpy 2.4.6 and 2.5.3, scipy 1.18.0 and 1.18.1,
  scikit-learn 1.8.0 and 1.9.1, pandas 3.0.2 and 3.0.6, pydeseq2 0.5.4 in both.
- At 20,000 x 200, wheels plus data are 56,385,918 bytes, over the 50,000,000 limit, so that size
  fails even with a working fit.
- Fallback arm: in a Pyodide where `import pydeseq2` raises `ModuleNotFoundError`, scipy 1.18.0,
  pandas 3.0.2 and plotly 7.1.0 installed as 10 files of 32,226,431 bytes. Mann-Whitney U with BH
  matched CPython at every size: max |dp| 2.2e-16, max |dpadj| 5.6e-16, and no gene crossed 0.05.
- Recorded beside the rule, which never reads it: with a sequential joblib backend registered as
  `in_process`, pydeseq2 fits in Pyodide in 6.6, 38.5 and 44.8 s of Node time. It agrees at
  1,000 x 24 (0 disagreements, max |dLFC| 2.2e-5) and at 20,000 x 200 (0, 2.8e-4) but not at
  20,000 x 24: 389 disagreements. 388 are low-count genes, baseMean 0.27 to 1.34, whose padj is
  NaN in Pyodide only (CPython has none NaN); one of them, g05499, has CPython padj 0.0307. The
  389th is well expressed: g15471, baseMean 947.8, p 0.01027 against 0.003241, padj 0.0876
  against 0.0318. The largest p difference is 0.0381 (g01437, baseMean 1208.6). So the fits
  differ beyond independent filtering. Recomputed on 2026-09-26 from the run's raw results files.
- The record's `wheels_downloaded_bytes` holds the bytes of unique package files by name, not
  the bytes downloaded: formulaic, interface_meta and python_dotenv were each fetched twice
  (41,957,904 bytes of wheels fetched against 41,800,855 recorded), and 785,537 bytes of index
  pages and `.metadata` files are excluded. The Pyodide runtime is excluded as well:
  `pyodide.asm.wasm` 9,598,218 and `python_stdlib.zip` 2,545,564 bytes (npm package 314.0.6).
- The fallback's own budget at 20,000 x 200: 32,226,431 bytes of wheels plus 14,585,063 bytes of
  input leave 3,188,506 bytes of the 50,000,000 for the precomputed grids, before the runtime.
- Known gap in the pre-registered rule, left unedited: when both pydeseq2 version rows are
  missing, the version criterion compares None with None and passes. The record carries both
  versions, so this verdict is unaffected; a rerun must fix the rule first, in its own commit.
- Sources, retrieved 2026-09-26: https://registry.npmjs.org/pyodide (314.0.6 integrity; npm
  `latest` is 314.0.7, published 2026-09-14);
  https://registry.npmjs.org/@jupyterlite/pyodide-kernel/-/pyodide-kernel-0.8.6.tgz (devDependency
  `pyodide` 314.0.6); https://pypi.org/pypi/anndata/0.13.4/json (`typing-extensions>=4.16`).

**Consequence.** Phase 5 bundles precomputed DE grids and runs the scipy tests live. Reversing this
needs a new decision: the backend shim changes the notebook's code path, 20,000 x 24 still
disagrees under it, its wheels alone take 41.8 MB of the 50 MB budget, and no browser timing
exists. A numpy change that alters `counts()` fails `tests/test_spike_record.py`; the spike is
then rerun. The proof that the rule predates the record needs full history, so it runs locally
and skips in CI's depth-1 checkout.

## D-004 — Output root: the laptop profile defaults to the working directory; hpc needs an explicit root

**Date:** 2026-09-26 · **Status:** decided

**Context.** `settings/paths.toml` refuses guessed roots, because on the cluster a fallback
writes where nobody looks. Students run `ctx` from a course folder on a laptop.

**Options.** Always require `--root` or `CTX_ROOT`; default to the working directory everywhere;
default to it on the laptop profile only.

**Decision.** `--root`, then `CTX_ROOT`, then the working directory on `laptop`. `hpc` raises
without an explicit root. The resolved root is printed and written into every `_README.md`.

**Evidence.** None measured; this is a preference. The cluster risk is the one the template's
`settings/paths.toml` comment records.

**Consequence.** P1.6 tests the precedence. Reversing it changes where every student output lands.

## D-003 — ctx does not import labcore at runtime (TODO.md D15)

**Date:** 2026-09-26 · **Status:** decided

**Context.** `AGENTS.md` routes every reusable helper through the pinned `labcore` package.

**Options.** Depend on labcore after upstream changes (inline plotly, a theme toggle, a unit
column), or keep ctx self-contained with labcore in the environment for `labdocs` only.

**Decision.** Self-contained. Chosen by Rytis at the 2026-09-26 check-in.

**Evidence.** In labcore v0.6.6: `viz/save.py:147` writes plotly HTML with
`include_plotlyjs="cdn"`, one file per theme; the `write_table` README dictionary is
`Column | Dtype | Description`, with no unit (`io.py:239-273`); a CSV write is refused unless
`settings/formats.toml` lists it or `override=True` is passed; the Okabe-Ito order at
`viz/palette.py:60-61` puts `#CC79A7` fourth and `#F0E442` seventh, the reverse of
`CONVENTIONS.md`.

**Consequence.** Conformance level 5 ("labcore used for viz/IO") is unreachable. P1.2 adds an
AST test that fails on any `labcore` import under `src/ctx`.

## D-002 — Plan files moved to the template's locations; project rules live in CONVENTIONS.md (TODO.md D14)

**Date:** 2026-09-26 · **Status:** decided

**Context.** Copier owns `CLAUDE.md` and does not skip it when it exists. `TODO.md` and
`knowledge/lessons.md` are skip-if-exists, and on this case-insensitive disk `TODO.md` already
resolved to the root `todo.md`.

**Options.** Move and import (chosen); merge the project rules into `AGENTS.md`; keep the files at
the root and hand-merge `CLAUDE.md`.

**Decision.** `todo.md` became `TODO.md`, `lessons.md` became `knowledge/lessons.md`, and
`CLAUDE.md` became `CONVENTIONS.md`, which the generated `CLAUDE.md` imports. Where
`CONVENTIONS.md` conflicts with `AGENTS.md`, `CONVENTIONS.md` wins. `/data/`, `/resources/` and a
`tests/fixtures/**` re-include were added to `.gitignore`.

**Evidence.** Commit 45632ed is exactly three `R100` renames. The `copier copy` log (v0.9.10) has
0 overwrites and `conflict`/`skip` for `TODO.md` and `knowledge/lessons.md` only; all three files
kept their SHA-256.

**Consequence.** `copier update` merges against two short project additions in `CLAUDE.md` and
`AGENTS.md` instead of the whole rule set. They sit at the top of each file, after the import and
after the intro paragraph, not at the end: a real update to a synthetic template release that
appended to both files left conflict markers in each while the lines were appended, and merged
cleanly once they were moved (P1.11, 2026-09-26). A template change to those top lines would still
conflict; that is the smallest surface available.

## D-001 — `src/ctx` is an installable package, an exception to "no shared helper module"

**Date:** 2026-09-26 · **Status:** decided

**Context.** `AGENTS.md` puts every runnable step in `scripts/<lang>/` and forbids a shared
helper module in the repo; `labdocs` walks only `scripts/`, `bin/` and `containers/`.

**Options.** A typer CLI package in `src/ctx`, or scripts only.

**Decision.** The teaching product is the `ctx` package in `src/` (`CONVENTIONS.md`, Stack).
One-off steps such as fixture builders and spikes stay in `scripts/python/` with codebase-meta
blocks.

**Evidence.** None measured; this is a requirement: students install and run one CLI.

**Consequence.** `labdocs` does not document `src/ctx`, `knowledge/api_index.md` does not index
it, and the file-length ratchet does not cover it. Conformance stays at 1 until raised on purpose.
