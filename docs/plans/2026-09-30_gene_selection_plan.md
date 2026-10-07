# Gene selection that drives M6, M9 and M10: final plan

2026-09-30. This plan builds Angle A. It adds C's unthinned volcano, survival oracle and symbol rules, and B's per-module entry pages, genes-by-contrasts evidence view and build-list notice. B's browser-storage channel is dropped. Nothing in the repository was changed.

**Source tags.**
- [m]: measured for this plan by the scripts in [genesel/plan/](/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/genesel/plan/): `measure.py`, `sizes2.py`, `ids.py`, `rcheck.R` and `phprev.py`.
- [map]: the mapping findings.
- [A], [B] and [C]: the three design prototypes, in [angleA](/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/genesel/angleA/), [angleB](/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/genesel/angleB/) and [angleC](/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/genesel/angleC/).
- [J1] and [J2]: the two judges, in [judge/](/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/genesel/judge/).
- "Arithmetic" means a sum of measured parts; no such page was built.
- UNVERIFIED means nobody measured it.

## Outcome

A student selects genes on the M5 volcano with a lasso, a box or a click, ticks them in a table of every tested gene, or pastes names. The M6, M9 and M10 pages then recompute from that selection: live over-representation analysis (ORA) on M6, per-gene evidence and cards on M9, and a signed score with Kaplan-Meier, log-rank, Cox and the cut-point lesson on M10.

The selection travels in the URL fragment of ordinary links. That is the one channel that behaved the same on file:// in Chrome 138, Edge 154, Chromium 152 and Firefox 156 at its default policy. An optional "Open beside" window follows every change through postMessage, which was also measured in all four. Safari is untested.

## How selection works

**Identity.** A gene is its unversioned Ensembl number, never a Plotly curveNumber or pointNumber, because `liveCalls` rebuilds the volcano traces on every cutoff change [map].
- On Step-1, the 17,830 tested IDs of the 10 contrasts give 17,830 distinct numbers, with no _PAR_Y and no non-ENSG ID [m].
- The page's customdata holds versioned IDs such as ENSG00000123119.14, so the version is stripped once, on the way in [B].

**Direction.** Each selected gene carries the sign of its raw log2 FC in the page's featured contrast, set when the gene enters the selection. On Step-1 ndn_gcsf, 6,958 of the 12,796 tested genes have log2 FC < 0, and none has exactly 0 [m]. M6 uses the sign for its up and down lists, and M10 for a signed score, so a mixed lasso does not cancel [A]. A list pasted on a page without log2 FC for its genes is unsigned, and the page says so.

**One list, not seven paints.** Figure 2 of the QC page paints samples in seven colours to compare groups [map]. Here M6 tests one list against a background and M10 splits one score, so a page holds one selection, and direction rides on the sign.

**Controls on M5.** These follow the controls-kind pattern of D-033: a new kind, `select`, beside `brush`, `calls` and `filter`. The code is a partial `_selection.html.j2`, a pure object `ctxSelect` in `_select.js.j2` tested under node, and wiring `selecting(form)` in `_page.js.j2`, run inside `guard()`.

- **Volcano.**
  - A radio group, "Drag to: Zoom, Box, Lasso", sets `layout.dragmode` through `Plotly.relayout`. The modebar stays hidden, as PLOTLY_CONFIG sets it [map].
  - A second radio group, "Replace, Add, Remove", decides what a drag does. `plotly_selected` applies the drag and `plotly_click` toggles one gene.
  - Selected genes are drawn as rings on an overlay trace keyed by gene, in a new style role, `selected`, at 3:1 or more on both themes. Ring and size mark a selection as well as colour does.
  - Plotly's own dimming is cleared after each pick, as `brushing()` does.
  - Real drags on the built Step-1 volcano returned gene IDs:

    | Browser | Lasso | Box | Source |
    |---|---|---|---|
    | Chrome 138 | 129 of 134 drawn points | 129 of 134 | [A] |
    | Firefox 156 | 134 of 134 | 134 of 134 | [A] |
    | Chrome, Edge and Firefox | added one gene | not checked | [C] |

    The cause of Chrome's gap is not known.
- **Every tested gene can be lassoed.** Today the page draws 9,112 of 12,796 tested genes: 7,270 with p < 0.1, plus 1,842 thinned by stride from the other 5,526 [map]. The volcano's trace JSON then ships with empty arrays, and the page fills x, y and customdata from one column table of all tested genes, as `liveCalls` already rebuilds them [C]. The PDF and SVG are unaffected: thinning applies only on the page (`scatter.py:337`, read).
- **Genes table.**
  - Drawn by JS from the same column table, 25 rows a page, with a search box, column sort and a "selected only" filter [C].
  - Columns: Select, Gene, log2 FC, padj, and the call at the live cutoffs. That is within CONVENTIONS' 8-column cap, and adds no bytes.
  - "Select the current up calls" and "Select the current down calls" take the calls from `ctxCalls` at the reader's thresholds. This is the DE-thresholds-feeding-ORA knob of M6 in TODO.md.
- **Paste box.**
  - Accepts symbols or Ensembl IDs with any separator, in any case, versioned or not.
  - An ambiguous symbol is refused and listed with its IDs: 27 symbols name 55 tested genes on ndn_gcsf [m]. C's first, symbol-keyed version picked the wrong gene 4 times in 3,703 [C].
  - Names the page cannot resolve are listed by name.
- **Undo, Clear and the count.**
  - Undo keeps 50 steps, and Clear can be undone.
  - One `aria-live="polite"` line reads, for example, "42 genes selected: 30 up, 8 down, 4 not called at the current cutoffs; 3 names not found: X, Y, Z".
- **Leaving the page.**
  - Copy list writes one ID and symbol per line, with the existing select-the-text fallback on file://.
  - Save list writes a `.txt` file with `#` header lines for source, time and count, which `--genes-from` reads [B]. Whether a download works from file:// is UNVERIFIED, so Copy is the primary path.
- **Onward links.**
  - One real `<a>` per selection-consuming module, for example "Open these 42 genes in M6". Its href is rewritten on every change.
  - Beside each link, an "Open beside" button.
  - All of these work with the keyboard and with middle-click.
- **Accessibility.**
  - The lasso needs a pointer, so the keyboard route is the table plus the paste box.
  - Gene names from a link or a message are inserted only with `textContent`.
  - Tables scroll inside their own wrapper and never hide overflow, per lessons.md.
- **Page words** go in `ui.toml [select]`, and the prose goes in content/*.toml, after P2S.4 moves M5's text there.

**Controls on M6, M9 and M10.** The same panel appears without the volcano: table, paste, Undo, Clear, count, links and Open beside. Each page adds its own picker:
- M6: clicking a gene-set row makes that set's tested members the selection.
- M9: a table of the selected genes; a clicked row shows that gene's card.
- M10: a row of the per-gene Cox screen, flagged as chosen on these patients.

**Build-list notice [B].** When the live selection differs from the list the PDF and SVG were drawn for, the page says so, for example: "The PDF and SVG show the build list, 7 genes; this view shows your list, 42 genes."

## What each module shows and recomputes

### M5 (Phase 2D)

**Live on each change:** the ring overlay, the table's checkboxes and "selected only" filter, the count line by live call, and the hrefs of the onward links. No new statistic is computed. The calls form keeps working as today, and the rings stay valid through its restyles, because a ring sits at the gene's log2 FC and -log10 p, which the call rule does not change [A].

**Scope:** only the featured contrast has a volcano. On Step-1 that is ndn_gcsf, with 12,796 tested genes, 1,683 up and 2,020 down [m]. A volcano for each of the 10 contrasts would need 8,888,916 bytes of column arrays [map], so it is not offered.

### M6 (Phase 4)

**Live, in `ctxOra` (`_ora.js.j2`).** It reproduces `ctx.function.lite.ora` in the same order:
1. The universe is the chosen background.
2. Hits are cut to the universe.
3. Each set is cut to the universe, then kept at 15 to 500 genes.
4. P(X ≥ k) is `hypergeom.sf(k-1, N, K, n)`.
5. BH runs over the kept sets [map].

Genes map to symbols as `_symbols()` does: 12,796 tested genes on ndn_gcsf give 12,768 symbols [m], and 48 of the 50 Hallmark sets pass the cut [map].

**Measured agreement with `lite.ora` on Step-1:**

| Prototype | Largest relative difference in p | Time per run |
|---|---|---|
| [A] | 1.52e-10 | 14 runs in 20.8 ms |
| [B] | 1.04e-10 | 0.7 to 3.1 ms |
| [C] | 1.09e-10 | 0.8 to 1.6 ms |

Sets, sizes, hits and overlap genes were identical in all three.

**Knobs.**
- Direction: all selected, selected up, selected down. Each direction is its own BH family, as `run.py` tests up and down calls separately [map].
- Background:
  - the tested genes of any of the 10 contrasts, from bitsets; the default is the featured contrast;
  - the gene-set union plus the tested genes, labelled as the background trap: 14,477 genes against 12,768 on ndn_gcsf [map].

**Shown:**
- a live bar chart of -log10 FDR for the top sets, with fold enrichment k/expected and k;
- a table of every kept set: size after the cut, hits, expected, fold enrichment, p, FDR and hit genes;
- the selection's make-up in the source contrast: n up, n down, n not called;
- a copyable TSV of every number shown.

**Default with no fragment:** the featured contrast's up and down calls, which is what the static ORA figure and CSV show today. A test asserts the two are equal [A].

**Not live:** GSEA and its running-score plots. They use the whole ranking, not a list [A]. The panel marks the sets that contain selected genes [C].

### M9 (Phase 4, after M10, in Phase 2M's order: M10, M9, M8, M7 [TODO.md Phase 2M])

**Tier 1, any gene tested in any Step-1 contrast: 17,830 genes [m].**
- An evidence view of genes by contrasts: log2 FC on the P2Q.6 diverging scale, called cells marked, at most 100 genes drawn and the rest in a table [B].
- One evidence row per gene: log2 FC raw and shrunk, padj, and the D-012 call in all 10 contrasts.
- For genes tested in the featured contrast, a box of the fit's 26 samples, 13 steady_state and 13 gcsf [C], following M5's box convention in TODO.md: Mann-Whitney U, BH across the selected genes that have a box, rank-biserial r, and DESeq2 padj beside it with the note on why they differ.

**Tier 2, full cards: only genes in the build's card list.**
- The card list is the dataset file's registered genes, the top 25 up and top 25 down calls of the featured contrast, and any `--genes-from` list [A].
- Cards are drawn on demand with `Plotly.react` from compact per-gene data, into one set of card figures [A].
- The builder packs cards until `MAX_BYTES` and names how many did not fit.
- Bytes per card are UNVERIFIED until the external-data spike reads real GTEx and HPA responses. The preview's made-up card is 282,031 bytes [B], which does not have the real shapes.

**A gene without a card** shows its tier-1 evidence, links opened only on click (GTEx Portal, versioned `v<NN>.proteinatlas.org`, Ensembl), and the command that builds its card, `ctx genecard --genes-from <file> --dataset <slug>`, with Copy.

D22 holds: the page fetches nothing, and reference data come from `ctx genecard`'s cache at build time [TODO.md M9].

### M10 (Phase 4, breast root first, then the Step-1 root on the cohort of decision 1)

**Patients.**
- TCGA-BRCA: 1,111 primary-tumour files from 1,095 cases.
- One tumour per patient, taken by first aliquot barcode, joined to TCGA-CDR.
- CDR-redacted cases removed, as CONVENTIONS requires: 1,089 patients with PFI, 145 events [A]. B and C kept 4 of the 5 redacted BRCA cases, all censored, in their 1,093 [J1].
- Endpoints per CDR Table 3: PFI by default, and DFI with 84 events [A, C].

**The score, one definition in `ctx.clinical.score` with a JS twin, in this order:**
1. log2(CPM + 1) from unstranded STAR counts. Library size is the sum over genes, without the N_ rows.
2. Store as int16 of the value × 1000. The Python build and the page read the same rounded values. At this precision no patient crossed the median, and HR moved by at most 1.3e-5; uint8 moved 4 to 12 patients [A].
3. Map to GENCODE v36 by unversioned ID.
4. Drop genes with zero variance, and genes below a minimum expression, `--min-log2cpm`. The dropped genes are listed by name. The default is proposed at a median log2 CPM of 1, as C prototyped. Of the 3,201 mapped featured calls, 802 fall below it and 54 have zero variance [C]. The build item measures what the default drops, including any registered gene, before fixing it.
5. Compute z per gene over the analysed patients.
6. The score is the mean over the present selected genes of sign × z.
7. Standardise the score to SD 1, so each HR is per SD.

**Live, in `ctxSurv` (`_surv.js.j2`):**
- a pre-specified median split, by numpy's linear percentile, with tertiles as an option;
- Kaplan-Meier with log(-log) Greenwood 95 % bands, censor marks and a synced number-at-risk table;
- median survival with its 95 % CI, and median follow-up by reverse Kaplan-Meier;
- log-rank chi-squared and p;
- RMST difference at the horizon, 60 months by default, as the preview uses [`preview_clinical.py:54`, read];
- univariable Cox per SD of the score, and Cox adjusted for age and AJCC stage I to IV coded 1 to 4: 1,065 complete cases, 136 events [A]. Subtype joins when the PAM50 source of Phase 4 exists;
- the proportional-hazards (PH) test as R's `cox.zph(transform = "rank")`;
- Harrell's C;
- the cut-point scan from the 20th to the 80th percentile, whose smallest p is shown only beside its corrected p of the maximally selected statistic: 1,000 permutations, mulberry32 with seed 42 shared by Python and JS, and log-rank scores by cumulative sums.

**Fitting choices.**
- Ties by Efron, R's `coxph` default [m]. Breslow and Efron differ by at most 6.9e-6 in β on these data [A].
- Newton steps are halved when the likelihood falls.
- A monotone likelihood is detected and refused, never reported as converged. The preview's Cox fit gives β = -256.8 for ENSG00000248769. R 4.4.2 with survival 3.7.0 warns "coefficient may be infinite" and returns β = -0.394 with SE 52.4 [m].

**Timings.**

| Step | Time | Source |
|---|---|---|
| Kaplan-Meier, log-rank and Cox together | 0.4 ms, node | [A] |
| Adjusted Cox | 4.6 ms, node | [A] |
| Scan, 1,000 permutations | 109 to 114 ms in node; 132 to 153 ms in Chrome, Edge and Firefox | [C] |

**Precomputed at build.**
- The per-gene Cox screen over the panel, univariable and adjusted, with BH q. A fitted 3,150 genes in 0.7 s [A]. The table shows the screen filtered to the selection.
- The power map over events and HR, with an 80 % isoline and the endpoint's events marked [TODO.md M10].

**Always shown.**
- The cross-study label, for example: "genes chosen in Step-1 blood cells, contrast ndn_gcsf, scored in breast tumour bulk; not chosen on these patients".
- A list picked from M10's own Cox screen is flagged as chosen on these patients.
- A counter of the lists scored since the page opened.

**Gene coverage.** Step-1 lists shrink in GENCODE v36, and the page states the coverage of every list [A]:

| Step-1 list | Genes | In GENCODE v36 |
|---|---|---|
| Tested | 17,830 | 15,881 |
| Featured calls | 3,703 | 3,201 |
| Top 50 by padj | 50 | 47 |
| Registered box genes | 7 | 7 |

The panel embedded for BRCA holds the featured calls that pass the filter. Any other gene in a selection is named as "no TCGA data on this page".

## Transport

**Fragment format, one definition in `ctx.report.genes` with its JS twin in `ctxSelect`:**

```
#v=1&src=<dataset>:<contrast>&n=<count>&g=<numbers>
```

`g` is the unversioned Ensembl numbers in decimal, sorted and comma-separated. A leading `-` marks a gene with log2 FC < 0 in `src`. This is plain decimal, not A's delta-varint codec, because plain numbers need only a split and a parse, and the largest list still fits by a wide margin [J2].

A receiving page:
- checks `v`;
- decodes strictly and compares the decoded count with `n`, refusing a cut-off link with "n expected, m decoded" rather than half-using it;
- refuses a fragment over 500,000 characters [A];
- writes its own fragment on every change with `history.replaceState`, so reload, bookmarks and Back keep the selection. A replaceState of a 5,708-character fragment worked, and Back restored it [A].

| List on Step-1 | Genes | Plain decimal, characters [m] | A's delta codec, characters [m], equal to [A] |
|---|---|---|---|
| Registered box genes | 7 | 48 | 23 |
| Top 50 by padj | 50 | 342 | 134 |
| Featured calls | 3,703 | 25,537 | 5,698 |
| All tested, ndn_gcsf | 12,796 | 88,134 | 17,438 |
| Union of 10 contrasts | 17,830 | 123,034 | 24,134 |

Signs add one character per down gene, at most 6,958 on ndn_gcsf [m].

**Entry pages.** Pages are dated and M5 is built before M6, so M5 cannot name M6's file. Each selection-consuming module therefore writes, on every build, `outputs/NN_<module>/index.html`. The page is a few lines: `location.replace(<its newest overview> + location.hash)`, plus a plain `<a href>` to the same page, which serves no-JS readers and gives the kit's link walk a target.
- It is always as new as the module, unlike a router written once by `ctx report`, which goes stale when a module is rebuilt after the report [J1].
- **Measured:** a 500-gene list arrived through an entry page in Chrome and Firefox, and `history.length` stayed 2, so Back returns to M5 [B]. A 24,149-character fragment arrived through a forwarding page in all four browsers, and Back returned [A].
- M5 renders an onward link to module X only when X consumes a selection, which is a fact of the code, and the dataset's cohort bundle of Phase 2D lists X, which is configuration.
- A listed module that was not built makes the kit build refuse. `build_kit.py` follows links through `reached()`, and `check_kit.link_problems` reports a missing target; the fragment is ignored when resolving [read; J1].
- The file name needs a naming-lint exemption: decision 3.

**Open beside.**
- A button opens the entry page with `window.open(url, <module name>)`. From then on, each change is sent by `postMessage` to every window that this page opened or that opened it, except the window the change came from. Windows opened this way form a tree, so no message loops.
- The message is `{ctx: "genes", v: 1, src, n, g}` with targetOrigin `*`, which file:// requires [map, MDN].
- A receiver accepts a message only from `window.opener` or a window it opened. It decodes strictly, applies the list as Replace with an Undo entry, and announces "Selection changed on M5: 42 genes".
- A tab not opened this way keeps its own selection and shows its own source line; nothing claims to be shared.
- Opener and popup exchanged `postMessage` both ways on file:// in all four browsers [map]. B's Open beside exchanged changes both ways in Firefox at its default policy [B].
- Whether the popup's opener link survives the entry page's `location.replace` is UNVERIFIED. The browser check tests it.

**Behaviour per browser on file://**

| Channel | Chrome 138 | Edge 154 | Chromium 152 | Firefox 156, default | Safari 26.4 | Source |
|---|---|---|---|---|---|---|
| Fragment through a link to another folder, byte-exact at 1,000,015 characters | yes | yes | yes | yes | UNVERIFIED | [A] |
| Longest fragment that arrived | 2,000,015 | 2,000,015 | 2,000,015 | 1,000,015; 1,200,014 refused [B]; 2,000,015 failed [A] | UNVERIFIED | [A], [B] |
| Fragment through a forwarding or entry page, Back returns | yes | yes | yes | yes | UNVERIFIED | [A], [B] |
| Reload keeps a 24,149-character fragment | yes | yes | yes | yes | UNVERIFIED | [A] |
| postMessage between opener and popup, both ways | yes | yes | yes | yes | UNVERIFIED | [map] |
| gzip block decoded by DecompressionStream | 37 to 41 ms | 37 to 41 ms | 37 to 41 ms | 37 to 41 ms | UNVERIFIED | [A]; 56 ms in Chrome 138 [J1] |
| localStorage across folders; not used | shared | shared | shared | one store per file | UNVERIFIED | [map] |
| BroadcastChannel and storage event; not used | delivered | delivered | delivered | not delivered | UNVERIFIED | [map] |
| fetch or XHR of a sibling file; not used | refused | refused | refused | refused | UNVERIFIED | [map] |

**Why storage is not used.** In Firefox at its default policy, an open M6 or M10 tab kept its old list after the list changed on another page, and no warning appeared [B's firefox_bar.json; J1]. The same class exercise would then give different results in Chrome and Firefox.

## Sizes

The cap is 15,000,000 bytes per page. Of that, plotly.js takes 4,815,817 bytes as embedded, against 4,815,814 on disk; the map attributes the difference to the U+FFFD escape [map].

**M5, Step-1**

| Part | Bytes | Source |
|---|---|---|
| Page today | 12,779,190 | [map], [m] |
| Volcano rows today: 9,112 row dicts, JSON body | -2,127,793 | [m] |
| Volcano trace JSON arrays, emptied: 838,908 down to 1,652 | -837,256 | [m] |
| All 12,796 tested genes as 8 column arrays at 4 significant digits; gene_id as a number | +721,676 | [m] |
| **Page before the new script and markup** | **10,535,817** | arithmetic |
| Headroom | 4,464,183 | arithmetic |

- The map gives the rows as 2,129,250 bytes, 1,457 more than the JSON body measured here; the cause was not checked.
- C's whole four-panel logic and wiring was 33,172 bytes [C], an upper reference for the new script.
- The standalone volcano figure page, 7,862,159 bytes [map], comes from the same object and is re-measured after the change.

**M6, Step-1**

| Part | Bytes | Source |
|---|---|---|
| Page today | 5,371,213 | [map] |
| Gene index: 17,830 numbers and symbols, JSON | +313,307 | [m] |
| Tested-gene bitsets for 10 contrasts | +29,751 | [m] |
| Hallmark v2026.1.Hs: 50 sets, 5,651 memberships, as index lists | +29,699 | [m] |
| Featured log2 FC and padj over the union at 4 significant digits | +244,241 | [m] |
| **Total before the new script** | **5,988,211** | arithmetic |

**M9, Step-1, on the preview's base: the real base is UNVERIFIED**

| Part | Bytes | Source |
|---|---|---|
| Preview page | 5,183,733 | [map] |
| Gene index | +313,307 | [m] |
| log2 FC at 3 and padj at 2 significant digits, 10 contrasts over the union | +1,983,413 | [m], equal to [B] |
| Box block: 26 samples, uint16, 12,796 genes | +1,080,743 | [C] |
| **Total before cards** | **8,561,196** | arithmetic |
| Left for cards | 6,438,804 | arithmetic; bytes per card UNVERIFIED |

**M10, Step-1 root scored in TCGA-BRCA, on the preview's base**

| Part | Bytes | Source |
|---|---|---|
| Preview page | 5,306,275 | [map] |
| 3,201 genes by 1,089 patients: int16 in two byte planes, gzip, base64 | +6,729,268 | [A] |
| Per-gene Cox screen | +101,490 | [A] |
| Gene index | +313,307 | [m] |
| Clinical columns: 4 endpoints, age and stage for 1,095 patients | +44,927 | [C] |
| **Total** | **12,495,267** | arithmetic |
| Headroom | 2,504,733 | arithmetic |

- The same block measured 9,295,704 bytes as plain base64 and 7,627,700 as plain gzip [A].
- The calls of all 10 contrasts, 5,096 genes in v36, would need 10,707,052 bytes and do not fit [A].
- The expression filter shrinks the block by an amount not yet measured.
- Without DecompressionStream the page says so and shows its build list from the static data.

**Build cost.**
- Reading the 1,095 TCGA files took 9.1 s at a peak of 1.44 GB [C]. The Track B budget is under 10 min and under 4 GB [CONVENTIONS.md].
- The kit is 785,336,079 bytes today [map]. M5 loses about 2.2 MB and M10 adds about 7 MB, both arithmetic.

## Statistics and annotation of live results

Live captions are templates in content/*.toml, embedded in the page and filled by JS from computed numbers. They pass `pixi run lint-prose` like any page text, and a test fails on an unfilled placeholder.

- **M5:** counts by call, descriptive; the caption says "no test".
- **M6:**
  - the one-sided hypergeometric test;
  - N, the background, named, with its size;
  - n, the selected genes in the background, and the selected genes outside it or without a symbol;
  - the number of sets kept after the cut to the background;
  - BH over the kept sets, and how many reach FDR < 0.05;
  - fold enrichment as the effect size;
  - the Hallmark release and licence;
  - the lists tried on this page.
- **M9:**
  - Own evidence keeps M5's inference: the pydeseq2 Wald test, BH within each contrast, and the D-012 call.
  - The box shows Mann-Whitney, BH across the selected genes with a box, rank-biserial r, and DESeq2 padj beside it.
  - Reference panels are descriptive only: n, median, IQR and tau, with release and retrieval date.
- **M10:**
  - the endpoint and its CDR Table 3 source;
  - n and events per group, and median survival with its CI;
  - log-rank chi-squared and p;
  - RMST difference and horizon;
  - HR per SD with its 95 % CI and Wald p, univariable and adjusted, with the covariates named;
  - PH p and C;
  - the split rule;
  - for the scan: the cuts, the smallest p and its cut, the corrected p, the permutations and the seed;
  - the cross-study label and the tries counter.
  - The power map shows power across events and pre-specified HRs, with the endpoint's events marked. It never prints power at the observed HR, which is a function of the observed p [J1].
- **Where a CONVENTIONS rule does not map onto the data:**
  - "Adjusted p (BH)" has no family for one model per list. M10 prints the raw p beside the count of lists tried and says it is not corrected across them. BH applies where a family exists: ORA sets, the Cox screen and the M9 boxes.
  - "Cross-study comparisons ... with batch diagnostics" assumes values compared across studies. Scoring a Step-1 list in TCGA compares no values across studies. The page shows gene coverage in v36 and expression in the scoring cohort instead, and says why.

## Static figures

The PDF and SVG pairs, white and black from one render, show the build's default list, drawn from the same data contracts as the live figures. Each keeps its full annotation inside the figure.

| Module | What the static files show |
|---|---|
| M5 | The volcano as registered, with no selection |
| M6 | The ORA of the featured up and down calls, as today; GSEA unchanged |
| M9 | The cards and the genes-by-contrasts view for the card list |
| M10 | For the build list, which on Step-1 is the 7 registered box genes of D23: KM with number-at-risk table, the Cox table, the scan and the power map |

A reader's own selection gets its journal pair, CSV and README from Python, using the same definitions and the same rounded values the page used:
- Copy list or Save list;
- `ctx function --genes-from`, `ctx clinical --genes-from` or `ctx genecard --genes-from`, which accept a plain text file of symbols or IDs;
- from Phase 2D, `ctx figure <id> --params` with a `genes_from` field, and `ctx verify-figure` comparing that figure's CSV with the TSV the page exports.

A browser cannot write a PDF without a server.

## Tests

Pure objects run under the pinned node 26.10.0 against Python references, as `tests/test_page_script.py` runs `ctxCalls` [map]. For each test, the mutation it must catch is named.

1. **`ctx.report.genes` and `ctxSelect`.** Round trips on the empty list, one gene, 17,830 genes, duplicates and versioned IDs.
   - The order test: strip the version, deduplicate, sort, then format.
   - Catches: the version not stripped; the sign dropped; deduplication missing; the `n` check skipped, so a truncated link is half-used; an ambiguous symbol resolved to its first gene, C's first error; a case-sensitive symbol match.
2. **`ctxOra` against `lite.ora` and `scipy.stats.hypergeom.sf`**, rtol 1e-9, where 1.5e-10 was measured.
   - Catches: the size filter before the cut to the background; hits not cut to the background; `sf(k)` in place of `sf(k-1)`; BH without its running minimum; BH rank off by one.
3. **`ctxSurv` against statsmodels 0.15.0**: PHReg with Efron ties, SurvfuncRight and survdiff. Also against a committed fixture from R 4.4.2 with survival 3.7.0, for `cox.zph` rank and concordance. R is not in pixi; the fixture records the R version and the script.
   - Catches: Breslow in place of Efron; the PH test ranking event times only, the preview's current definition. On J1's lists that p is 0.310 against R's 0.266 on the top 50, and 0.00143 against 0.000808 on the down calls [m].
   - Also catches: log-rank without the (n-d)/(n-1) factor; Greenwood bands on the plain scale; percentile "lower" in place of "linear"; an off-by-one scan bound.
4. **The score's order,** in `ctx.clinical.score` and its JS twin.
   - Catches: z computed before rounding; raw log2 CPM averaged without z; an unsigned mean on a signed list; standardising before the mean; the expression filter applied after z; redacted cases kept, the 1,093 of B and C; two tumours for one patient.
5. **Monotone likelihood.** A fixture with one non-zero value, censored, shaped like ENSG00000248769, must be refused. Catches: the check removed, when the fit reports β = -256.8 as converged [m].
6. **Permutations.** Python and JS mulberry32 with Fisher-Yates give the same corrected p, exactly. Catches: another generator, another seed or another shuffle order.
7. **Default state.** With no fragment, M6's live ORA equals its static ORA CSV, and M10's live results equal its KM and Cox CSV. Catches: a default list or background that differs from the build's.
8. **Unthinned volcano.** The page's column table holds every tested gene, and the volcano trace JSON carries no data arrays. Catches: thinning left on the page, which makes 3,684 of the 12,796 Step-1 genes unreachable by lasso. The existing budget test runs at Step-1 size, per lessons.md.
9. **Page structure.**
   - `test_each_kind_of_control_gets_its_own_form` gains `select`.
   - `test_html_offline` applies to entry pages as well, with no plotly and no external request.
   - Every new `ui.toml` key is used, and none shadows a dict method.
   - A card exists for each new CLI option.
   - `test_text_contrast` covers the `selected` ring colour.
   - The naming lint accepts `index.html` only in `outputs/NN_<module>/`.
   - `check_kit` refuses a missing entry page.
   - Catches: a page without its select form; an entry page that fetches; `index.html` accepted in another folder.
10. **Browser check.** An opt-in local node driver with no npm package: CDP for the Chromium family and WebDriver BiDi for Firefox, the method used for every measurement here.
    - The path: M5, then M6 through the entry page, then M10, then M9, then back to M5.
    - At each step it checks the fragment, Back and reload.
    - It runs a real lasso and a real box drag on the volcano.
    - It checks Open beside in both directions, through the entry page.
    - Firefox runs with the Remote Agent's recommended prefs off, or it silently runs with `strict_origin_policy=false`, and with `app.update.disabledForTesting` set [map].
    - Catches: a fragment lost in a redirect; a Back that loses the list; an opener link dropped by `location.replace`.
    - The DOM wiring has no CI test today, and this check does not add one to CI.

## Build order

Add these to TODO.md.

### Phase 2D, after the cohort-bundle item and P2S.4

This set replaces the existing line "the M5 page's gene picker feeding M9", and leaves the QC exclusion set and the volcano's x axis on their own lines.

- [ ] **GS.0 Record the contract.** D-048: the unversioned Ensembl number as identity, with a sign; one selection per page; the fragment format; entry pages; Open beside by postMessage; no browser storage. Rytis's answers to the three decisions below. Cox for a selection's score becomes live if decision 2 says so, a stated divergence from TODO.md M10's "Cox fits precomputed".
- [ ] **GS.1 `ctx.report.genes` and `ctxSelect`.** In `_select.js.j2`, with node tests 1 and the order test. `--genes-from` reads through the same parser. No page changes yet.
- [ ] **GS.2 Entry pages.** Each selection-consuming module writes `outputs/NN_<module>/index.html`. The naming lint exempts that one name in module folders. `build_kit` reaches entry pages by link, and `check_kit` refuses a missing one. Test 9's entry-page parts.
- [ ] **GS.3 Unthinned volcano.** `ScatterData` ships the page's volcano trace JSON without arrays, and the page fills them from one column table of every tested gene. `liveCalls` reads columns. Test 8. The Step-1 M5 page is rebuilt and measured against the arithmetic 10,535,817 bytes, and the rebuild is recorded as after the D-046 record.
- [ ] **GS.4 M5 picker.**
  - The `select` controls kind, the drag-mode and Replace/Add/Remove radios, the ring overlay, the genes table, the select-calls buttons, paste, Undo, Clear, the count line, Copy and Save list, the onward links and the build-list notice.
  - Words in `ui.toml [select]`.
  - The browser check with a real lasso and box, at 1280 and 375 px, in both themes.

### Phase 4

- [ ] **GS.5 M6 live ORA**, after P2S.4 moves M6's text to content/function.toml.
  - `_ora.js.j2` and test 2.
  - The embedded index, bitsets and Hallmark lists.
  - The direction and background knobs, including the trap.
  - A set row becomes the selection; the caption, TSV and tries counter.
  - `ctx function --genes-from`; test 7 for M6.
  - M6 writes its entry page, and M5's M6 link turns on.
- [ ] **GS.6 Open beside.** The postMessage relay in `ctxSelect` and its wiring. Test 10 covers both directions through an entry page in the four browsers. If the opener link does not survive the redirect, the button opens the dated page directly and says so.
- [ ] **GS.7 M10 selection, breast root first**, inside the existing M10 item.
  - `ctx.clinical.survival` and `ctx.clinical.score`, the preview's functions moved into `src` under one definition, with the PH test changed to R's rank definition.
  - `_surv.js.j2`; tests 3 to 7.
  - lifelines pinned in ctx-dev once 0.30.0 installs from conda-forge, which is UNVERIFIED.
  - The int16 two-plane gzip block with its DecompressionStream fallback.
  - `--min-log2cpm`, with its measured drops.
  - Monotone-likelihood refusal, the scan, the precomputed Cox screen and the power map without observed power.
  - The cross-study label, the double-dip flag and the tries counter.
  - Static figures for the build list; `ctx clinical --genes-from`.
- [ ] **GS.8 M10 in the Step-1 root.** The scoring cohort from decision 1 is named in the pancreas cohort bundle and declared `external: true` with its SHA-256. The panel is the featured calls in v36. The coverage table is on the page. M5's M10 link turns on.
- [ ] **GS.9 M9 selection**, after the external-data spike and the M10 items, in Phase 2M's order.
  - The tier-1 evidence view and the box.
  - Cards drawn on demand up to the measured byte cap.
  - The missing-card command; `ctx genecard --genes-from` taking a text file.
  - M9 writes its entry page.
- [ ] **GS.10 Proof and review.**
  - The kit rebuilt.
  - The selection carried M5 to M6 to M10 to M9 and back, in Chrome, Edge, Chromium and Firefox at its default policy, at 1280 and 375 px, in both themes, keyboard only.
  - Safari by Rytis, from a checklist.
  - The round trip added to K.4's Windows checklist.
  - Page sizes, wall time and peak RAM measured.
  - A review section in TODO.md; lessons.md updated.

## Not done, and why

**Not tested.**
- **Safari 26.4 is untested for every channel:** fragment length, replaceState, Back, postMessage and DecompressionStream. Automating it needs "Allow Remote Automation", a setting nobody changed [map].
- **Windows file paths are untested.** That waits for the K.4 PC.
- **A touch lasso** is untested.
- **The cause of Chrome's 129 of 134** lasso points was not found.

**Channels not built.**
- **Open tabs that were not opened through Open beside do not update.** Browser storage would do that in Chrome, Edge and Chromium but not in Firefox at its default policy [map], so it is not used.
- **The single explorer page of Angle C is not built.** Its budget is structural: survival data for 1,560 of 3,703 calls, and AIM2 cut until forced [C, J2]. It stays a later option.
- **Unresolved names do not travel in the link,** as A's `s=` would carry them. They are listed on the page where they were pasted, and one list keeps one parser.

**Genes outside a page's data.**
- **No live result exists for genes a page does not embed.** On Step-1, 1,949 of 17,830 tested genes and 502 of 3,703 featured calls have no v36 match [A]. On BRCA the M10 panel holds at most the 3,201 mapped featured calls. Cards exist only for the card list.
- **There are no volcanoes for the other 9 contrasts:** 8,888,916 bytes of column arrays [map].

**Analyses not recomputed or configured.**
- **GSEA is not recomputed for a selection,** because it uses the whole ranking.
- **Reactome is not configured.** Its GMT is 898,680 bytes by HEAD, and its sets are UNVERIFIED, so the live ORA is Hallmark only.
- **M6's two-knob contour,** which a click turns into thresholds, is enabled by the live ORA but not built here.

**Definitions left to their owners.**
- **The ORA background's 2,051 ENSG-named genes.** On ndn_gcsf, 2,051 tested genes carry an Ensembl ID as gene_name [m], and they sit in `lite.ora`'s background of 12,768, though no Hallmark set can hold them [C]. Whether they belong there is an M6 definition question, decided once in `lite.ora`. The live ORA copies whatever `lite.ora` does.
- **No correction across the lists a student tries.** The page shows the count and says so.

**Outside this change.**
- **The theme is lost between pages in Firefox**, a measured pre-existing fault [map]. It is not fixed here. A `t=` field in the same fragment would fix it on these links.
- **Rebuilt Step-1 M5 and M6 pages will not come from the code state recorded in D-046.** The recorded metrics are unaffected, and the rebuild is recorded [map, C].
- **lifelines is not installed.** Its PH test is UNVERIFIED against `cox.zph`, and so is its tie method.
- **The swimmer plot of TODO.md M10** is unchanged and separate.

**Tell Rytis.**
- Firefox on this Mac updated itself from 143.0.3 to 156.0.1 at 15:16 on 2026-09-30, during the mapping agent's first headless launch [map].
- The design scratch holds browser profiles and a 265,690,928-byte matrix [B], which can be removed.

## Decisions for Rytis

1. **Which tumours score a Step-1 gene list on M10.** Recommended: **TCGA-BRCA**, labelled cross-tissue.
   - It is on disk, with 1,089 patients and 145 PFI events after redactions [A].
   - M10 is built breast first.
   - The cohort is one bundle entry, not code, so it can change later.
   - Its cost is coverage: the embedded panel is at most the 3,201 mapped featured calls, about 6.7 MB gzipped [A].

   Alternatives:
   - **TCGA-PAAD**, Step-1's own cancer: 185 patients, 110 PFI events and 100 OS events in the CDR on disk [A, B]. It is not downloaded, and its download size is UNVERIFIED. At 185 patients every mapped tested gene would fit: 3,917,316 bytes at int8 [B]. At int16 that is about twice as much, by arithmetic.
   - **No M10 link from Step-1 pages.** Selection-driven survival would then exist only in the breast root.
2. **Cox for a selection's score, fitted live in the page.** Recommended: **yes**, univariable and adjusted, validated by tests 3 to 6. The per-gene Cox screen stays precomputed. This diverges from TODO.md M10's "Cox fits precomputed", because a score from an arbitrary selection cannot be precomputed. The alternative: Kaplan-Meier and log-rank live, and Cox only for the build list, so a selection shows no HR.
3. **The file name `outputs/NN_<module>/index.html`, exempt from the naming lint in module folders only.** Recommended: **yes**, recorded in D-048. Each module rewrites its entry page on every build, so links never point at an older page. The alternatives:
   - one undated router, `outputs/09_report/route.html`, written by `ctx report`. It needs the same kind of exemption and goes stale when a module is rebuilt after the report [J1];
   - no forward links, only Copy and Paste between pages. That works everywhere but loses "automatically".
