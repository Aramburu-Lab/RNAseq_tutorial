# Critique of the D-051 control proposals

Critic's pass over the 75 proposals of the six audit groups, 2026-10-01. Read only; every number below was measured on the kit at `$HOME/ctx_kit/2026-10-01_course_kit_v4/` or read from the repository at 2590650, unless it is marked as the proposer's or as an estimate.

## Outcome

| Verdict | Proposals |
|---|---|
| Approve as proposed | 27 |
| Approve with changes | 23 |
| Defer: measure first, or wait for P2Q.8 or a defect fix | 5 |
| Reject and replace with the control Rytis already named | 1 |
| "None" agreed | 5 |
| Already planned, cite only: M10 KM and cut-point scan | 2 |
| Planned figures with no page: record the control in the figure contract, nothing to approve | 12 |
| **Total** | **75** |

Every row and byte count I checked on the built pages is exact: 64 of 64 `ctx-data` row and byte pairs, and every derived count in section 1. The proposals are sound on data. Their weak points are elsewhere:

1. **No proposal can meet D-051's "the PDF and SVG buttons export the view as set" today, and none costs it.** No export code exists anywhere: the `qc-page-fixes` worktree holds four uncommitted files about band-map labels, and Plotly 4.1.1 exports png, jpeg, webp, svg and full-json, never PDF. This decision constrains every control, so it comes first. See 2.1.
2. **The control that already ships breaks the honesty rule.** The bar filter never refills the lead or the caption: on Step-1 with cell type NDN the plot shows 73 of 210 runs and a median of 26.2 M, while the lead keeps "135 of 210 runs have fewer than 30 million". "Words follow" has to be a contract of the kind, tested, before more controls copy the gap. See 2.2.
3. **Two M8 proposals contradict or miss what Rytis named on 2026-10-01** in [TODO.md](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md) lines 616 to 617. Five more overlap planned knobs without citing them. See 2.4.
4. **Three statistical traps**: a marginal batch PERMANOVA on a design where batch nests in condition, "in one level only" read as confounding when it is nesting, and live correlations on 12 samples with no p. See 2.5.
5. **Teaching values measured on made-up previews are planted values, not findings.** Byte costs measured at 33 samples do not hold at D-050's 210. See 2.6.
6. **Four control kinds cover all 50 approved proposals**, the existing brush included. The groups named six overlapping kinds between them. See section 3.

## 1. What I checked

Script: [critic/m.py](/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/interactivity/critic/m.py), written independently of the proposers' scripts. Output for every overview page: [critic/all_overviews.jsonl](/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/interactivity/critic/all_overviews.jsonl).

How the bytes were measured:
- Rows: the UTF-8 length of the `ctx-data-N` script text.
- Traces: the groups used two methods. M1 used compact re-serialised JSON. M2, M5 and M6 used the raw text on the page. I reproduced both, for example M1 Figure 1 at 1,368 compact against 1,488 raw, and M5 Figure 2 at 7,199 raw against 6,624 compact.
- The page cap counts raw bytes, so future audits should quote raw bytes.

| Claim | Page | Proposer | Measured | Verdict |
|---|---|---|---|---|
| M1 F1 rows | gather tier1 / teach / desynth | 35 rows 5,247 B; 8 rows 1,232; 6 rows 904 | same | exact |
| M1 F2 rows and controls | gather tier1 / teach / desynth | 210 rows 40,833, controls 3,500; 8 / 1,587 / 611; 33 / 6,322 / 920 | same | exact |
| M1 F2 filter defect, NDN | gather tier1 | 73 runs, median 26.2 M, 50 below 30 M | 73, 26.19 M, 50 | exact |
| M1 F2 runs below cut-off by batch | gather tier1 | 10 M 5/99 vs 1/111; 20 M 44/99 vs 18/111; 30 M 80/99 vs 55/111 | same | exact |
| M1 F2 two roundings | gather desynth | lead 0.5, live 0.494 | same | exact |
| M1 F1 groups of 22 in one level only | sample sheet | batch 9, tissue 20, isolation 13, time_point 16, sex 1 | same | exact |
| M1 per-run rows | sample sheet | 37,247 B | 37,247 B, or 16,764 B as column arrays | exact |
| M1 n by contrast filters | tier1.toml and sheet | ndn_gcsf 13/13, ndn_hsct_early 13/5 | 13 and 13 donors; 5 and 13 donors, 9 and 14 runs | exact on the 2 checked |
| M2 band map, linked, strips, distance, agreement, size-genes, sex, fates, flow rows | qc tier1 | 210/20,790; 210/39,128; 37,121 to 55,371; 209/48,679; 209/48,505; 210/31,922; 210/37,245; 210/105,244; 5/671 | same | exact |
| M2 flags | flags CSV | 182 rows, 107 samples, globin 53, pca 16, spearman 9, tg_pct_passing 21 | same | exact |
| M2 tg_pct_passing scaled MAD in LDN | metrics and sheet | 0.0059 points | 0.00594 | exact |
| M2 sex: own sex disagrees | sex check and sheet | 31 of 208: 18 F with M donor, 13 M with F donor | same | exact |
| M2 teach run band map | teach CSV | 11 bad, 13 usable of 72; 3 and 4 of 48 without kept metrics | same | exact |
| M2 fastp grid | teach grid CSV | six -l 102 cells keep 0; 25 in band 7, 5 in band 6 | same | exact |
| M2 strips shared table | sheet | 23,403 B | 23,403 B | exact |
| M5 rows, all 12 figures | de tier1 | e.g. volcano 9,112 / 2,127,793, MA 9,225 / 1,827,555, controls 1,457 | same | exact |
| M5 PERMANOVA R2 | metrics CSV | 0.684, 0.089, 0.026; ndn 0.398; ss_ndn 0.126 p 0.002 | same | exact |
| M5 first p bin over flat | metrics CSV | null 1.0, adjusted 1.7, naive 7.9, unpaired9 20.1, paired 31.7 | same | exact |
| M5 MA thinning breaks pairs | de tier1 | 1,645 of 5,435 genes | same | exact |
| M5 lesson totals | de tier1 | 701, 40, 47; 4,018 vs 3,457 | same | exact |
| M5 naive vs adjusted agreement | results CSVs | 659 naive calls, 27 survive, 13 new, 11,067 neither | same | exact |
| M5 volcano caption contradiction, "every 3th" | de tier1; [scatter.py:357](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/viz/scatter.py) | as quoted | same | exact |
| M6 rows, all 5 figures | function tier1 | 14/6,222; 15/3,348; 10/1,073; 92/6,520; 14/1,415 | same | exact |
| M6 ORA zero-length bars | function tier1 | 5 of 7 down bars at 1.56e-13 | same | exact |
| M6 GSEA passing per contrast, 10 contrasts | GSEA CSVs | 36/49 to 2/48 | same; 483 rows | exact |
| M6 null split interferon | GSEA CSV | IFN alpha NES -2.65 FDR 0, gamma -2.48 FDR 0 | same | exact |
| M6 every-contrast rows | GSEA CSVs | 76,970 B | 76,970 B | exact |
| M6 background reading line against the data | function tier1 | union makes most sets less significant | 6 of 7 shown sets less significant | confirmed on the page; 44 of 48 not re-checked |
| M3/M4/M7 rows, 14 figures | desynth | e.g. 1,000/158,282; 2,500/286,778; 2,108/452,613 | same | exact |
| M3 agreement: d never below 0 | process desynth | empirical 2.5th percentile 0.0 | min 0.0, 2.5th 0.0, 0 of 1,000 below 0 | exact |
| M4 Salmon-only median reads | quant desynth | 13.3 | 13.27 | exact |
| M8 neighbour Spearman, raw/residual | ncrna desynth | all 0.717/0.356, type_a 0.529/0.212, type_b 0.846/0.497, control 0.475/0.168, case 0.624/0.418 | same | exact |
| M8 kept by min count | ncrna desynth | 1,369/344 at 10; lncRNA 1,065, 861, 569, 195, 70 at 1, 2, 5, 20, 50 | same | exact |
| M9 tau | genecard desynth | 0.429 log, 0.812 linear; without T cells 0.295 and so on | same, all 8 | exact |
| M10 power | by the Hsieh-Lavori formula | 97.5, 91.0, 73.5 %; 81 and 176 events | same by arithmetic | exact |

Not re-checked: the proposers' M2 rejudge port beyond its k = 3 match with the flags table, the M2 peer counts, the M6 ORA precision claim and lite copy, the numpy running-score check, the M10 Cramér's V values and the M7 numbers beyond rows. I did not run pixi.

## 2. Cross-cutting findings

### 2.1 The PDF and SVG cannot follow a view, and the route needs a decision first

D-051 decided that "the PDF and SVG buttons export the view as set". Every proposal says this "waits for the qc-page-fixes in-browser export", and P2S.5 names it too. Nothing exists:
- [_figure.html.j2](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/report/templates/_figure.html.j2) renders the buttons as links to files written at build time.
- No template calls `toImage` or `downloadImage`.
- `qc-page-fixes` sits at 424459f, behind main. Its worktree holds four uncommitted files, about band-map row labels and the `tg_total_reads` title.
- Plotly 4.1.1, bundled in every page, exports png, jpeg, webp, svg and full-json. It never exports PDF.

| Route | How | Exact and on style | Covers | Cost |
|---|---|---|---|---|
| A. In the browser | Plotly SVG, plus a bundled SVG-to-PDF library and an embedded font | No: it exports the page-sized Plotly render, not the journal-sized Matplotlib one, against "PDF and SVG from one render" in CONVENTIONS | every kind | a PDF library and a font file in every page, not measured |
| B. Prerendered per option | build renders each option of a select | yes | `view` selects with few options only | e.g. 10 contrasts × 2 themes × PDF and SVG at the measured 23.8 kB PDF and 18 kB SVG of the p histogram: about 0.84 MB per figure on disk |
| C. Reproduce parameters | the control writes its state into the figure's params; `ctx figure <id>` renders it | yes | every kind | needs Phase 2D `FigureSpec`; already planned for M10 KM at TODO.md:294: "the Reproduce panel rewrites its params file to the current selection" |

My recommendation:
- Use B for `view` controls with up to about 10 options, and C for `subset` and `threshold` once Phase 2D lands.
- Until then, the buttons keep the default view's files, and a line beside them says the files show the default when the view differs.
- Use A only if Rytis accepts a second rendering style.

Whatever is chosen decides how a control is designed. Under B, a free numeric input cannot export, so fixed presets beat typed values. That is one reason I turn M1's free depth input into presets below.

### 2.2 The words must follow, as a contract of the kind

- **The bar filter already breaks this.** Verified by reading `filtering()` in [_page.js.j2](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/report/templates/_page.js.j2) lines 174 to 206, which writes only `.filter-summary`, and by the NDN numbers in section 1. The same gap sits on 6 pages: the M1 overview and library-size pages of Step-1, the subset and the made-up study.
- **The fix belongs in the kind.** Each kind carries lead, caption and reading templates per state in its controls JSON, read from `content/*.toml` like every page text (P2S.3, P2S.4). It is tested under node against the Python fill, as `ctxFilter` is in `tests/test_page_script.py`.
- **"What does not follow" goes in the live line, never the caption.** 48 of 108 rendered captions on the kit's overview pages already exceed 40 words by the lint's own rule, 33 of 52 on Step-1. The templates pass only because [prose.py](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/report/prose.py) counts a placeholder as one word. Most proposals add a sentence such as "the status table stays at the run's k" to the caption. Put it in the `aria-live` summary instead.
- **Glossary marks must survive a refill.** D-048 lets a caption carry `[[term]]` marks rendered by `explain_terms`. Refilling a caption's `textContent` would drop the icons. Refill placeholder spans only. No caption in the kit carries a mark yet, so this is a design rule, not a defect today.
- **One number, one rounding.** The made-up M1 page prints one median as 0.5 in the lead and 0.494 in the live line. The lead must use `ctxFilter.number`'s rounding, the twin of `ctx.viz.bar.significant`.

### 2.3 Precompute a finite choice in Python; port to JS only a continuous knob

Several efforts are rated large only because the proposal recomputes a statistic in the page. When the options are few, Python computes every option once, through the one definition, and the page looks the result up. That needs no JS twin, keeps one definition in one place, and makes route B possible.

| Proposal | Proposed | Instead | Effort |
|---|---|---|---|
| M2 Agreement peers | none, or large | Spearman per peer choice computed by `ctx qc`, about 18 kB, estimate | medium in Python, small in the page |
| M9 F1 batch subset | live Mann-Whitney, scipy `method=auto` twin | 3 options computed by `annotate.compare_groups` | small |
| M8 F4 subset p | live t-test twin, or caption caveat | computed per option | small |
| M8 F5 targets | live KS and BH | a few context-score cut-offs computed | medium |
| M10 F4 Cramér's V | live | 2 corrections × 2 strata computed | small |
| M5 F9, F11 live rho | live Spearman over carried rows | 2 options computed | small |
| M2 flow filters | per-sample fates, 78,099 B | flows per level of one column: 14 levels × 5 links, about 9.4 kB by arithmetic from 671 B per 5 links | small |

Keep live JS twins for continuous knobs:
- the volcano and other call rules;
- the M2 judge's k;
- the M8 minimum count;
- the M10 power formula;
- the running sum, whose net bytes fall;
- GS.5's ORA and GS.7's KM.

Each twin is tested in both languages, as `ctxCalls` and `ctxSelect` are.

### 2.4 Proposals that contradict or miss what is already named or planned

| Proposal | What is already there | Consequence |
|---|---|---|
| M8 F3 volcano adds a "BH over all tested genes" basis | TODO.md:616, Rytis: the M5 volcano's controls, "biotype colours kept for called genes and BH still within each biotype"; CONVENTIONS: "Non-coding: stratify by biotype" | The pooled basis departs from his spec and the statistics rule. Drop it, or put it to Rytis as a question. The proposal's own numbers show it changes little: 81/24 both ways at the default. |
| M8 F4 neighbour: values/residuals and sample subsets | TODO.md:617, Rytis: "a drop-down of lncRNA and neighbour-gene pairs, each named by both genes", Spearman on values and residuals recomputed, the preview simulating several pairs | The proposal misses the named control and says "None planned". The pair drop-down shows both values, which makes a values/residuals toggle redundant. |
| M5 F8, F10 calls grid, 5 bases × 4 cut-offs × 5 floors | M5 contour figure "two-knob maps (DE genes across padj × \|LFC\| cut-offs, and across min-count filter × padj)", TODO.md M5 section | The same grid. Build the data once for both, and share the calls state with the volcano. |
| M6 ORA cut-off knobs | M6 "Contour figure: two-knob map of significant pathways (ORA) across DE thresholds; a click sets the thresholds and reruns the live ORA"; P2S.5; GS.5 | Build only through GS.5's engine. A second ORA port would be a fourth copy of one rule. |
| M6 GSEA | M6 knob "ranking metric" | Not proposed. It needs permutations, so it can only be a few options computed by `ctx function`. Worth asking Rytis. |
| M5 F4 volcano: contrast choice rejected on bytes | M5 knob list: "contrast selector" | The rejection is measured and right for the overview page: about 2.9 MB per contrast against 2.22 MB of headroom. Rytis needs the choice set out: per-contrast pages, or a selector after GS.3. Even at GS.3's 0.7 MB per contrast, 10 contrasts come to about 16.8 MB by arithmetic, over the cap. |
| M8 F2 per-biotype min count | M5 knob "min-count filter (precomputed grid)" | One filtering rule for both; share the definition. |
| Every proposal's export | TODO.md:294, M10: export by the Reproduce params | Route C already exists in the plan. P2S.5 plans route A. Pick one. |

### 2.5 Statistical traps

- **M5 F2, batch on the NDN PCA.** The proposal computes batch's PERMANOVA there and calls its p "valid unstratified". The p is valid for a marginal hypothesis, but that is not the question. Every pdac and ipmn sample is in the > 76 nt batch: the figure's traces "ipmn, Reads ≤ 76 nt" and "pdac, Reads ≤ 76 nt" hold 0 points. So a marginal batch R2 absorbs condition. Report batch within condition, with permutations restricted to condition levels holding both batches, or sequentially after condition, and say which.
- **M1 F1, "a group found in one level only is confounded with that level".** For tissue, 20 of 22 groups, and time point, 16 of 22, this is nesting: BM cell types come from marrow, and steady state has no time point. Confounding matters for the two groups a contrast compares. Say per contrast what differs between its two groups. I recommend a contrast preset, below.
- **M8 F4 and M9 F1 subsets.** A live Spearman over the 12 type_b samples, ρ 0.846 by my count, carries no p or interval. Precompute the p per subset, or do not offer subsets that small.
- **M5 F6 p histogram.** ndn_pdac_adjusted's first bar at 1.7 times flat against naive's 7.9 is not "more correct". Adjusting for a batch nested in condition removes power and signal alike. The words must say that.
- **M2 strata knob.** Across six strata choices the samples flagged on any metric stay between 105 and 120, by the proposer's port. Which metrics flag moves much more: globin 53 to 38. The live line should give per-metric counts, not only the total.

### 2.6 Preview pages: made-up data and 33-sample bytes

D-050 removes the made-up pages at the release, and M3.6 builds its figures with controls from the start. So the M3, M4 and M7 to M10 proposals are input to those modules' figure contracts, not changes to preview pages.

**Their teaching values are planted.** Examples:
- the 1.4-point untrimmed effect;
- "every sample called reverse";
- pooled vs within BH at 81/24;
- a tau of 0.43.

Re-measure on real data before any of these justifies a control.

**Several byte figures fail at 210 runs:**
- M4 correlation heatmap: about 7.49 MB as rows, the proposer's estimate. It needs column arrays, as GS.3 plans for M5.
- M4 detection at three thresholds: about 1.4 MB, the proposer's estimate.
- M3 per-sample counter vectors: cannot be carried, the proposer's own finding.
- M2 base quality at P2Q.8's 210 runs: about 1.24 MB of rows, by arithmetic from 121 B per row measured × the proposer's 10,218 rows.

Column arrays cut M1's per-run rows from 37,247 to 16,764 B, measured. That makes them the default encoding for every new carried table.

### 2.7 Shared state across figures, and standalone pages

Several proposals couple figures:
- the M2 judge panel driving strips, the band map and the fates band;
- one contrast and one gene set across M6;
- the M5 calls rule across the volcano and the lesson grids.

Every figure also has a standalone page that carries that figure alone. Define one page-level state, keyed by name: `calls`, `judge`, `peers`, `contrast`, `gene_set`. A figure's controls read and write that state. A standalone page shows the same controls on its one figure.

Drop couplings a standalone page cannot honour. The fates figure's band "following the page k" is one: keep that band at the run's k and say so in the live line.

### 2.8 Defects found on the way: fix at the root, before or regardless of any control

All verified on the kit unless marked.

| Where | Defect | Root-cause fix |
|---|---|---|
| M1 F2, all filter pages | the lead and caption ignore the filter | kind-level words, 2.2 |
| M1 F2, made-up | one median in two roundings | one rounding |
| M5 F4 | caption says "all drawn", then "draws 1,842 of the 5,526"; "every 3th"; an English sentence written in code at [scatter.py](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/viz/scatter.py) lines 355 to 359, which P2S.4's file list does not include | move it to content, fix the ordinal, add scatter.py to P2S.4 |
| M5 F5 | thinning keeps one of a gene's two estimates for 1,645 of 5,435 genes | thin by gene |
| M5 F7 | the fitted trend carries all 12,796 points, 481,483 B by the proposer's measure, for a smooth line | thin the trend line |
| M6 ORA | 5 of 7 down bars have zero length | show passing sets, or the sets with any hits |
| M6 GSEA | reading line 3, "the others are shown for context", points at nothing: all 15 drawn sets pass | a words fix |
| M6 background | reading line 2 and UNION_CHECK say the union makes overlaps "look rarer", while 6 of the 7 shown sets look less significant under the union | rewrite both, whatever control is built |
| M6 every contrast | hover prints fdr unformatted, e.g. 0.13016621743036838 | a hover format |
| M6 running score | hover rank `%{x:.3~g}` prints 13,456 as 1.35e+04, with no symbol | a hover format |
| M6 background caption | "the same 1681 up calls": no thousands separator, against M5's 1,683 up calls | format; the 2-gene gap is not explained by me |
| M2 size-genes | lead "look for those outside their group's dashed line", while 102 of 107 flagged samples are flagged on neither plotted metric; the 5 are the library-size flags | default the triangles to the two plotted metrics |
| M2 fastp grid | 8 equal bands over 0 to 100 %, stretched by six -l 102 cells longer than the 76 nt read, put 30 cells in the top two bands | fix the bands or the grid range; the isoline knob would be a workaround |
| [qc_metrics.toml:149](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/conf/qc_metrics.toml) | fastp rationale says the settings are ones "the page lets you move"; no control does | correct the sentence until Phase 2D's click-to-set exists |
| M2 tg_total_reads title, band-map hover "Registry metric." | flagged by the M2 group | already changed, uncommitted, in the `qc-page-fixes` worktree: title "Processed reads" |
| M3 F5, made-up | caption omits the empirical percentiles the plan's F4 gives | caption |

## 3. Shared control kinds

Four kinds cover every approved proposal, the existing one included. The groups' "filter", "group", "view", "choose" and "keyed views" are different names for the first two.

| Kind | The reader | Computed | Replaces | Figures |
|---|---|---|---|---|
| `subset` | keeps rows: one select per column with All, or a checkbox set | in the page, by a named aggregate: count, sum, median, quantiles, or a twinned statistic | today's `filter`, plus words, a `measure` series for the median, checkbox sets and figure types beyond bars | M1 F2; M2 run band map, flow; M3 F1; M8 F1; M9 F2, F3, F4; M7 F1 alternative |
| `view` | picks one of a few named alternatives: a carried column to colour, slot, order or plot, or an option computed in Python | lookup only, no statistic in the page | the groups' group, view, choose and keyed views | M1 F1; M2 fates, sex, size-genes, Agreement peers; M5 F1 to F3, F6, F8 to F12; M6 GSEA, every-contrast, background, running-score set; M3 F2 to F5; M4 F1 to F5; M7 F1, F2, F4; M8 F5; M9 F1; M10 F4 |
| `threshold` | moves numeric cut-offs, presets first, against one rule ported to JS once | in the page, by a twin tested against Python under node | today's `calls`, generalised to any rule | M5 volcano; M8 volcano; M7 F3; M1 F2 depth; M2 strips judge, Distance peers, band-map flag view; M6 ORA via GS.5; M8 F2; M10 F3 |
| `brush` | paints across linked views | in the page | unchanged; adds a detail select | M2 linked PCA |

Every kind carries:
- the words contract of 2.2;
- the page-level state of 2.7;
- one export route from 2.1.

A `view` with up to about 10 options also carries a static file per option under route B.

## 4. Verdict per proposal

Effort: proposer's, then mine where different. "Small" means reusing a built kind; "medium" means a new kind or a twin; "large" means computation the page cannot do.

### M0, M1, hub, course home

| Figure | Verdict | Kind | Effort | Why |
|---|---|---|---|---|
| M1 F1 Samples by cell type and condition | change | view, subset | medium | Keep "Colour by" for batch and isolation, plus "Count: runs or individuals". Replace the three "Show only" selects and the multi-level time point with one "Contrast" select that applies each contrast's declared filters from [tier1.toml](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/conf/datasets/tier1.toml), showing the n M5 tests with. Free filter combinations build groups no contrast tests. Drop tissue and time point as colours: they nest in cell type and condition. Offer sex only as "recorded sex, the recipient's": 31 transplanted samples express their donor's sex, by M2's check. Words per contrast, as in 2.5. |
| M1 F2 Sequencing depth per run | change | subset, threshold | medium | First make the existing filter refill the lead and caption, on all 6 pages. The cut-off: presets "none, 20, 30 M" from `total_sequences`' sourced bands, not a free 0 to 231.8 input. A sourced value can be exported. Live line per batch as proposed: 80 of 99 against 55 of 111 at 30 M. |
| M0 Setup | none agreed | | | Tables of single measurements; no plotted data. |
| Course hub | none agreed | | | Right: a movable threshold would undo the pre-registered rules. |
| Course home | none agreed | | | Navigation only. |

### M2 quality control

| Figure | Verdict | Kind | Effort | Why |
|---|---|---|---|---|
| Absolute bands, Step-1 | change | threshold, as a view of the judge state | medium; build after the strips | On Step-1 it is one row of 210 good cells, so the judge view is the figure worth having: a flag map. Make the flag map the default, with the absolute band as one row of it. Carry its 292,213 + 20,661 B on the overview only, read from the same state as the strips. |
| Samples on the PCA | approve | brush | small | "Count painted samples by" changes `cfg.detail`, which `_brush.js.j2` line 35 already reads. Defer the PC pair: PC3 is unmeasured. |
| Strips, 13 metrics | approve | threshold | medium | The best M2 control. Twin of `checkpoint.judge`; the proposer's stdlib port matches the shipped flags at k = 3 on all 13 metrics. Per-metric counts in the live line, as in 2.5. "What does not follow" in the live line. |
| Distance from its group | change | threshold, shared `peers` | medium | One Peers select for Distance and Agreement, so the two group metrics use one meaning of "group". |
| Agreement with its group | change | view, shared `peers` | large to medium | Precompute Spearman per peer choice in `ctx qc`; the "none" is avoidable. |
| Library size against genes detected | change | view | medium to small | First fix the default as in 2.8; then the select is optional. Reject the grouping option: about 300 kB per grouping by the proposer's estimate, for little over the strips. |
| Recorded sex against expression | approve | view | medium to small | The strongest single lesson in M2: 31 of 208 disagree on own sex, all sex-mismatched transplants. Label the alternative as the rule the method avoids. Build on `ctx.qc.identity` after the move. |
| Read fates | change | view | medium to small | Colour by is good. Drop the band "following the page k": the standalone page has no judge panel. |
| Where every read goes | change | view | medium to small | One "Flows for" select, flows precomputed per level, about 9.4 kB, instead of per-sample fates at 78,099 B. |
| Base quality, subset | defer | view, subset | | Admitted low value at 8 runs. Redesign at P2Q.8's 210 runs, where rows reach about 1.24 MB by arithmetic. |
| fastp grid | defer | | | (a) isoline: reject; it works around the band defect of 2.8. (b) is Phase 2D's planned click-to-set; cite it. (c) per-run grids: wall time unmeasured; P2Q.8. |
| Run band map, subset | approve | subset on rows | medium to small | "Metrics the registry keeps": 11 bad and 13 usable of 72 become 3 and 4 of 48, verified. It teaches the first FastQC misreading. |
| Run strips, subset | approve | threshold, same instance | small once built | No extra work. Its teaching point, 2 runs per stratum, ends at 210 runs. |

### M5 differential expression

| Figure | Verdict | Kind | Effort | Why |
|---|---|---|---|---|
| F1 PCA, PB | approve | view | medium to small | Colour by cell type, condition, batch, with R2 lookups: 0.684, 0.089, 0.026; p only where computed. |
| F2 PCA, NDN | change | view | medium to small | Batch PERMANOVA within condition or after it, not marginal, as in 2.5. |
| F3 PCA, steady-state NDN | defer | view | | The components select teaches only if PC3 shows something, which nobody measured. Otherwise record "none: one factor, one sample per donor". |
| F4 Volcano | none new agreed | calls | | Fix the caption defects of 2.8. Put the planned contrast selector to Rytis with the numbers in 2.4. |
| F5 MA | defer | | | Fix thinning by gene first. The top-N mark moves 4 of 50 genes on ndn_gcsf; a second contrast for it costs about 0.84 MB, the proposer's estimate. |
| F6 p-value histogram | approve | view | medium to small | About 1.8 kB, estimate, for the null, batch and pairing lessons on one figure: 1.0, 1.7 against 7.9, 20.1 against 31.7. Words per 2.5. |
| F7 Dispersion | defer | view | | Conditional by the proposer's own terms: measure whether the cloud moves first. Thin the trend line regardless. |
| F8 Calls, batch lesson | change | view | medium to small | Selects over a grid computed in Python, so every count is exact. Share the `calls` state with the volcano, and build the grid once with the planned two-knob map. 14,658 B, proposer's measure. |
| F9 naive against adjusted | change | view | medium to small | Agreement colours; ρ over called genes precomputed: 0.734 and 0.867, proposer's. Counts verified: 659, 27, 13. |
| F10 Calls, pairing lesson | approve | view | small | With F8. |
| F11 NDN against Monocyte hsct | approve | view | small | With F9. |
| F12 Named genes | change | view | medium to small | Approve the family select, precomputed: 6 against 3 of 7 pass. Reject the gene ticks: choosing the family after seeing the data is the very pitfall, and ticks need a BH twin no page has. Defer the contrast option: it touches R4 and D-046. |

### M6 function

| Figure | Verdict | Kind | Effort | Why |
|---|---|---|---|---|
| Over-representation | change | threshold via GS.5 | large | Build through P2S.5 and GS.5 only, with the cut-offs reading M5's `calls` state. "Passing sets only" fixes the zero bars now. The contrast select at default cut-offs is small: 243,719 B, proposer's. |
| Gene set enrichment | approve | view | small | Contrast and "sets shown" over the 483 rows. Share one copy with the next figure: 123,151 and 76,970 B are the same rows. |
| Interferon alpha in every contrast | approve | view | small | Gene-set select. One panel with P2S.5's FDR select on the same figure. |
| Running score | approve | view, plus a running-sum twin | medium | Set select first; net bytes fall by 82,742, proposer's. Contrast select once the GSEA contrast exists; all 10 contrasts cost 2,218,605 B, proposer's. The page's ES check against the table is a good guard. |
| The background decides the p-value | approve | view | small | Four backgrounds precomputed, 81,536 B, proposer's. The whole-genome background is the page's mistake 1 made visible. Rewrite the contradicting words regardless, as in 2.8. |

### M3, M4, M7: made-up previews, as input to their module contracts

| Figure | Verdict | Kind | Why |
|---|---|---|---|
| M3 F1 Where reads align | approve | subset, view | Needs the `measure` series: a 100 % stack's median total is always 100 %. The limit toggle is a view. |
| M3 F2 Settings | approve | view | The change view is arithmetic on 24 carried rows. |
| M3 F3 Strandedness | change | view | The forward-share view, yes. Editable thresholds are low priority: one strandedness per cohort means nothing moves until the reader passes 0.08. |
| M3 F4 Counters | change | view | Gene type, yes; needs a `gene_type` column the preview lacks. Per-sample select: defer, it cannot be carried at 210. |
| M3 F5 Agreement | approve | view | The parametric lower limit no gene can reach, verified, is a real lesson. |
| M4 F1 Mapping rate | approve | view | Tests precomputed per grouping. |
| M4 F2 Detection | change | view | Needs one absolute depth grid in the module; bytes at 210 about 1.4 MB, estimate. |
| M4 F3 Correlation | change | view | Only with column arrays: 7.49 MB as rows at 210, estimate. |
| M4 F4 Units | approve | view | The planned `counts_from_abundance` knob; high value. |
| M4 F5 Tool detection | approve | view | 72 rows, about 8 kB, estimate. |
| M7 F1 Quantified | approve | view | The CONVENTIONS 70 % rule made visible. |
| M7 F2 RNA-protein ρ | approve | view, subset | |
| M7 F3 Fold change | approve | threshold | `calls` with two layers; thin by the loosest setting. |
| M7 F4 Paired heatmap | change | view | Class select, yes. Scatter-to-heatmap link: reject, about 13.2 MB as rows, estimate. |

### M8, M9, M10

| Figure | Verdict | Kind | Effort | Why |
|---|---|---|---|---|
| M8 F1 Biotype shares | approve | subset | small | With the `measure` fix. |
| M8 F2 Expression filter | change | threshold | medium | Carry the total-count rule, 56,000 to 59,500 B by the proposer's arithmetic. Without it, every count below 10 is an upper bound the page would print as exact. Share the rule with M5's planned min-count grid. |
| M8 F3 Volcano by biotype | change | threshold | medium | Rytis's named control, TODO.md:616: biotype colours, BH within biotype. Drop the pooled basis, as in 2.4. |
| M8 F4 Neighbour | reject, replace | view | medium | Build Rytis's pair drop-down, TODO.md:617. Values and residuals both shown; no subsets of 12, as in 2.5. |
| M8 F5 Target shift | change | view | large to medium | A few context-score cut-offs precomputed with KS and BH; nothing live. |
| M9 F1 Own evidence | change | view | large to small | Three batch options computed by `compare_groups`. |
| M9 F2 Normal tissues | approve | subset | medium | GTEx per-sample sex availability is unverified; the external-data spike. |
| M9 F3 Immune cells | approve | subset, view | medium | tau 0.43 to 0.81 with the transform, 0.30 without T cells, verified. A tau twin is trivial. |
| M9 F4 Exon use | approve | subset | medium | The "qPCR pair" key number must follow, or the live line says it does not. |
| M10 F1 KM | planned | | | M10 Interactivity, GS.7. Agree: z carried exactly. |
| M10 F2 Cut-point scan | planned | | | M10 cut-point lesson, GS.7. Agree on the corrected line per range. |
| M10 F3 Power | approve | threshold | medium | Closed form, verified at 97.5, 91.0 and 73.5 %. The 0.001 preset needs the HPA source written beside it. |
| M10 F4 Cramér's V | change | view | medium to small | Precompute 4 states. Verify the Bergsma 2013 citation before it ships; the proposer flagged it as from memory. |
| M9 subcellular, blood, structure | none agreed | | | No knob in the method. |
| 12 planned figures with no page: M8 violin, miRNA volcano, miRNA-target, tri-layer; M9 isoforms, junctions, HPA panels, cross-cohort row; M10 expression by clinical variable, genome-wide screens, competing risks, swimmer and KM bands | record | | | D-051's consequence: each new figure contract states its control. Write these into the module sections of TODO.md when the figure is designed; approving them now approves nothing. |

## 5. Page budgets for the approved set

| Page | Today | Added by approved controls | After | Source |
|---|---|---|---|---|
| M1 overview, Step-1 | 4,980,019 B | per-run rows 37,247 B, or 16,764 as columns | about 5.0 MB | measured |
| M2 overview, Step-1 | 6,764,836 B | flag view 312,874; shared table 23,403; peers about 41 kB; size-genes 5,939; sex 6,147; fates 18,682; flows about 9.4 kB | about 7.18 MB | proposers' measures and estimates, mine for the flows |
| M5 overview, Step-1 | 12,779,190 B | grid 14,658; PCA R2 and p lookups a few hundred B; p histogram about 1.8 kB; trend thinning saves about 474 kB | about 12.3 MB before GS.3 | proposers', estimate for the saving |
| M6 overview, Step-1 | 5,371,213 B | GSEA rows 123,151 shared by two figures; running set select -82,742; background 81,536 | about 5.49 MB, ORA cut-offs excluded | proposers' measures |

All four stay well under 15 MB. The M6 ORA cut-offs on every contrast would add 3,406,267 B, the proposer's measure: about 8.9 MB with them.

## 6. Order of work

1. Rytis decides the export route, 2.1.
2. Fix the defects in 2.8, starting with the filter's words and the M6 background text.
3. Add the words contract and page-level state to the shared page script; extend `filter` into `subset`.
4. Add `view` and build the cheap, high-value lookups: M5 F6, F8 to F11; M6 GSEA, every contrast and background; M2 sex.
5. Add `threshold` with the `checkpoint.judge` twin: the M2 strips, Distance and the flag view; then M1's presets.
6. Then P2S.5 and GS.5 for M6 ORA, the M8 controls as Rytis named them, and the module contracts for M3, M4 and M7.

## 7. What I did not do

- I did not re-run the proposers' ports for M2 rejudge and peers, the M6 lite copy and running sum, or the M5 lesson grid. I checked their outputs only where section 1 says so.
- I did not measure PC3, any dispersion table beyond ndn_gcsf, or the 44 of 48 background claim beyond the 7 sets on the page.
- I did not run pixi, did not build a page, and did not test any control in a browser.
