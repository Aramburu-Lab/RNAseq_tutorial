# D-051 controls: the plan for approval, P2S.6

**Of 75 audit entries, 59 are built figures and 3 of those carry a control today; this plan proposes a control for 51 figures, 2 of them extensions of today's controls, defers 5, cites 2 whose controls TODO.md already plans, writes the control into the contract of 12 planned figures that have no page, and finds 5 with none: 3 pages hold no plotted data, 1 figure has its control and its remaining knobs are planned, and 1 planned figure's method has no knob.**

| Module | Entries | Control today | Proposed | Deferred | Planned elsewhere | Contract only, no page | None |
|---|---|---|---|---|---|---|---|
| M0, hub, course home | 3 | 0 | 0 | 0 | 0 | 0 | 3 |
| M1 gather | 2 | 1 | 2 | 0 | 0 | 0 | 0 |
| M2 quality control | 13 | 1 | 11 | 2 | 0 | 0 | 0 |
| M3 trim, align, count | 5 | 0 | 5 | 0 | 0 | 0 | 0 |
| M4 quantification | 5 | 0 | 5 | 0 | 0 | 0 | 0 |
| M5 differential expression | 12 | 1 | 8 | 3 | 0 | 0 | 1 |
| M6 function | 5 | 0 | 5 | 0 | 0 | 0 | 0 |
| M7 RNA versus protein | 4 | 0 | 4 | 0 | 0 | 0 | 0 |
| M8 non-coding | 9 | 0 | 5 | 0 | 0 | 4 | 0 |
| M9 gene card | 9 | 0 | 4 | 0 | 0 | 4 | 1 |
| M10 clinical | 8 | 0 | 2 | 0 | 2 | 4 | 0 |
| **Total** | **75** | **3** | **51** | **5** | **2** | **12** | **5** |

The 51 split into 26 on Step-1 pages and 25 on made-up previews, whose values are planted. Reasons for none: no plotted data on M0 setup, the course home and the hub, where a movable threshold would also undo the pre-registered rules; the M5 volcano already ships its calls control; the M9 subcellular, blood and structure panels have no knob in their method. Deferred: M5 PC3 and the M5 dispersion design wait for a measurement, the M5 MA plot and the M2 fastp grid for a defect fix, and M2 base quality for P2Q.8's 210 runs.

Sources: the six inventories and [critique.md](/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/interactivity/critique.md), which found 68 of 68 row and byte claims exact. Bytes are measured on the kit `$HOME/ctx_kit/2026-10-01_course_kit_v4/` unless marked est. or by arithmetic.

## What changed since the critique

- **The in-browser export now exists, uncommitted.** The critique says no export code exists. That no longer holds: in files dated 14:23 to 14:34 today, the `qc-page-fixes` worktree holds [export.py](/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/qcfix_wt/src/ctx/report/export.py), [_export.js.j2](/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/qcfix_wt/src/ctx/report/templates/_export.js.j2), a `vendor/` folder with jsPDF, svg2pdf.js and Liberation Sans, and `tests/test_export_browser.py`. It is the critique's route A. Once a reader moves a control, the PDF and SVG buttons redraw the page's traces at the figure files' size and type through a journal template, white theme only, with the build-time caption followed by the live line. At the defaults they keep linking the Matplotlib files.
- **It costs 1,607,890 B on every page whose figure has a control**, measured by repeating `page_assets()` on the vendored files: 507,629 B of scripts and 1,100,261 B of base64 fonts. The journal template is not counted. The worktree's `save.py` line 175 attaches export settings to every calls, brush and filter control, so the M5 overview gains it at the merge, before any new control: 14,387,080 B by arithmetic, 612,920 B under 15,000,000 B.
- **The export covers three figure contracts only.** `save.py` line 161 gives controls to `ScatterData`, `LinkedData` and `BarData`. Every control this plan puts on a strip, heatmap, density, flow, curve, box or surface needs `controls()` and the export hook on that contract first.
- **The inventory's "D-041" for the derived multimapper limit is not a decision** in [decisions.md](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/knowledge/decisions.md). The limit is recorded under M3.4, [TODO.md](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/TODO.md) line 521. Rytis's named M8 controls sit at TODO.md lines 617 and 618 in the working tree.

### Page budgets

| Overview page, Step-1 | Today, B | At the export merge, B | After this plan's controls | Basis |
|---|---|---|---|---|
| M1 gather | 4,980,019 | 6,587,909 | about 6.60 MB | per-run rows as columns, 16,764 B, replace 5,247 B |
| M2 quality control | 6,764,836 | 8,372,726 | about 8.78 MB est. | flag view 312,874 B; sample table 23,403 B; peers about 31 kB est.; sex, fates, size-genes 30,768 B; flows about 9.4 kB est. |
| M5 differential expression | 12,779,190 | 14,387,080 | about 13.93 MB est. before GS.3, about 12.1 MB est. after | trend thinning frees about 474 kB est.; lesson grid 14,658 B; GS.3's page at about 10.5 MB by TODO's arithmetic |
| M6 function | 5,371,213 | 5,371,213, no control today | about 7.34 MB est.; about 10.75 MB est. with ORA cut-offs on every contrast | GSEA rows 123,151 B; running set net -82,742 B; background 81,536 B; ORA contrasts 243,719 B; cut-offs 3,406,267 B |

Every page stays under 15 MB. M5 is the tight one, so its trend thinning lands with or before the export merge. Across the kit, if each of the up to 124 module pages that are neither setup nor a hub carried a control, the export adds about 199 MB by arithmetic, against 769,012 kB for the whole kit today by `du`.

## Shared control kinds

Four kinds cover every proposed control, the existing brush included. The groups proposed eight new kind names; group, row filter, view, choose and keyed views fold into `subset` and `view`; judge and rule fold into `threshold`; running is a `view` with a twin.

| Kind | The reader | Computed | Built from | Cost to build, once | Figures |
|---|---|---|---|---|---|
| `subset` | keeps rows: a select per column with All, a contrast preset, or a checkbox set | in the page, by a named aggregate: count, distinct count, sum, median, quantiles, tau | today's `filter`, [bar.py](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/viz/bar.py) and `_filter.js.j2` | medium: words, a `measure` series, checkbox sets, contracts beyond bars, BarData's duplicate check scoped per level | M1 F1, F2; M2 run band map; M3 F1; M7 F2; M8 F1; M9 F2, F3, F4 |
| `view` | picks one of a few named alternatives: a carried column to colour, slot, order or plot, or an option Python computed | lookup only | new | medium: per-option traces and words | M1 F1; M2 size-genes, sex, fates, flow, Distance and Agreement peers; M5 F1, F2, F6, F8 to F12; M6 all five; M3 F1 to F5; M4 F1 to F5; M7 F1, F2, F4; M8 F4, F5; M9 F1, F3; M10 F4 |
| `threshold` | moves numeric cut-offs, presets first, against one rule | in the page, by a JS twin tested against the Python definition under node | today's `calls`, [scatter.py](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/viz/scatter.py) `Calls` and `_calls.js.j2` | medium per rule: count below a line; `checkpoint.judge`; calls with a colour column or two layers; FDR cut; min count with the total-count rule; Hsieh-Lavori power | M1 F2; M2 strips, run strips, Distance, band map; M5 volcano; M6 ORA through P2S.5 and GS.5; M7 F3; M8 F2, F3; M10 F3 |
| `brush` | paints across linked views | in the page | D-033, [linked.py](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/viz/linked.py) | small: a "count painted by" select | M2 linked PCA; later P2Q.7's M5 PCA brushing |

Every kind carries:

1. **Words that follow.** Lead, caption and reading templates per state come from `content/*.toml` and are filled by the page, tested under node against the Python fill as `ctxFilter` is in `tests/test_page_script.py`. What a control cannot keep true goes in the live line, which the export appends to the caption. Refills touch placeholder spans only, so D-048 marks survive. One number has one rounding, `ctxFilter.number`, the twin of `ctx.viz.bar.significant`.
2. **One page-level state**, keyed `calls`, `judge`, `peers`, `contrast` and `gene_set`. A standalone page carries the same controls for its one figure; a coupling it cannot honour is dropped. A figure whose options are precomputed follows the shared state only when the state is one of its options, and its live line names the option it shows otherwise.
3. **Precompute a finite choice in Python** through its one definition; port to JS only a continuous knob, each twin tested in both languages. This takes six proposals from large or live statistics down to lookups: M2 Agreement peers, M9 F1 batch, M8 F4 p, M8 F5 targets, M10 F4 Cramér's V, M5 F9 and F11 ρ.
4. **The export of a changed view** through `ctxExport`, with the Matplotlib files kept for the default view.
5. **Column arrays for every new carried table**, as GS.3 plans for M5: M1's per-run rows are 37,247 B as rows and 16,764 B as columns.
6. **Defaults that reproduce today's figure**, so untouched controls draw and export what the page does now. Three defaults change to fix a defect, named in the tables: the M2 band map, M2 size against genes and the M6 over-representation bars.

## Controls by module

Effort is counted with each kind already built, its own cost being in the table above: small adds data and words to a kind; medium adds a rule, a twin or a computation in the module; large needs computation the page cannot do now.

### M0 setup, course hub, course home

| Figure | Control | Recomputes | Data added | Effort | Teaching value |
|---|---|---|---|---|---|
| No figure on any of the three pages | None: single measurements in tables, links and two counts; the hub's sanity rules were registered before the data were seen | nothing | none | none | The hub's fixed thresholds are the point of its rules. |

### M1 gather

| Figure | Control: default and range | Recomputes | Data added | Effort | Teaching value |
|---|---|---|---|---|---|
| Figure 1: Samples by cell type and condition | Contrast: all samples (default) or one of the dataset's contrasts, applying its declared filters from `tier1.toml`; Colour by: batch (default), isolation, recorded sex, the recipient's; Count: runs (default) or individuals. Tissue and time point are not offered as colours, see question 6 | stacked counts, legend, hover; lead "{groups} groups; the smallest, {name}, has {n} {runs or individuals}"; live line naming what differs between the contrast's two groups, and saying paired, null-split and donor-list contrasts are not a filter | per-run rows, 16,764 B as columns, replacing 5,247 B of aggregated rows | medium: contrast filters exported from the dataset file, a distinct-count aggregate | Individuals under a contrast's filters give the n M5 tests with, which runs overstate in 6 of the 7 contrasts a filter reproduces: Monocyte hsct is 19 runs from 8 individuals. |
| Figure 2: Sequencing depth per run | Filters cell type and condition, All by default, kept and fixed so the lead and caption follow; Depth line: none, 20 M, 30 M (default), from `total_sequences` plot_bands | bars, axis, median line, lead, caption's batch medians over the runs shown; live line "{below} of {shown} runs shown are below {cut} million: {b1} of {n1} in {batch 1}, {b2} of {n2} in {batch 2}", saying M1 removes no run | none; controls about 300 to 600 B est. | small | At 30 M, 80 of 99 runs of the 101 nt batch fall below against 55 of 111 of the 75 and 76 nt batch, so a quality rule can create the confounding it should prevent. |

Not proposed: TODO's M1 knob "minimum group size", whose default has no source; question 12.

### M2 quality control

| Figure | Control: default and range | Recomputes | Data added | Effort | Teaching value |
|---|---|---|---|---|---|
| Absolute bands, every metric that has one, sample map | Judge by: stratum median ± k MAD over every judged metric, the absolute band one row of it (default, question 7) or absolute band only; reads the `judge` state | cell colours, rows, lead "{flagged} of {cells} values lie outside their stratum's median ± {k} MAD on the bad side", 182 at k 3, equal to the flags table | every sample-level value 292,213 B and bands 20,661 B | small once the judge twin exists | On Step-1 one of 15 sample-level metrics has an absolute band, so the dataset's own rule does nearly all the judging. |
| Samples on the PCA, beside their labels and QC status | Count painted samples by: cell type (default), condition, batch, QC status | the live line's counts, through `cfg.detail` | none | small | Paint a cluster and count it by batch to see whether it is the read-length batch rather than the condition. |
| Strips, 13 metrics, tab group "Every metric, sample by sample" | One judge panel for the tab group: k 2, 2.5, 3 (default, the run's `--mad`), 3.5, 4, 5; strata cell type (default, the run's `--strata`), none, cell type and condition, batch, cell type and batch | slots, each stratum's grey bar, triangles, hover reasons, lead; live line with per-metric flag counts, saying the status table, key numbers, interpretation and confounding tests stay at the run's k and strata | shared sample table 23,403 B; per-metric config under 3 kB est.; reason templates under 1 kB est. | medium: the `checkpoint.judge` twin, whose stdlib port reproduces the shipped flags of all 13 metrics at k 3 | Samples flagged on any metric go 147, 125, 107, 92, 79 and 61 at k 2 to 5, and condition in the strata moves globin flags from 53 to 38. |
| Distance from its group | Peers: cell type and condition (default, `--peers`), cell type, cell type and batch; distances precomputed per choice by `checkpoint.pca_distance`, judged live by the panel | triangles, hover, lead; samples without a peer drawn as missing and counted | 2 extra columns of 209 values, about 13 kB est. | small with the judge | Flags at k 3 are 16, 12 and 17 by peer choice: wider groups score biology as disagreement. |
| Agreement with its group | Peers, the same state; Spearman per choice precomputed by `ctx qc` | as Distance | about 18 kB est.; the extra Spearman passes' wall time is not measured, the whole Step-1 run took 49.4 s | medium: computation in `ctx qc` | Shows whether disagreement comes from a sample or from grouping different conditions as peers. |
| Library size against genes detected | Triangles flagged on: library size or genes detected (new default, fixing a lead that points at the dashed line), any metric, one named metric. Grouping select rejected, about 300 kB per grouping est. | drawn samples, symbols, lead "{flagged} of {n} samples are flagged on {metric}"; contours unchanged | flagging metrics per sample, 5,939 B | small | Of 107 flagged samples 5 are flagged on a plotted metric, and 28 of the 49 flagged inside their group's 50 % region are globin flags. |
| Recorded sex against expression | Compare with: the donor's sex for a transplanted sample (default, the method), or each sample's own recorded sex, labelled as the rule the method avoids; shown when the sheet has `donor_sex` | colours, agreement, named points, lead "{disagree} of {n} samples disagree with their recorded sex", the caption's transplant sentence | sex and donor sex, 6,147 B | small, on `ctx.qc.identity` after the record | By the method 0 of 208 disagree, against own sex 31 of 208, every one a sex-mismatched transplant expressing its donor's sex. |
| Read fates: raw against useful fragments | Colour by: cell type (default), condition, batch, isolation | traces, legend; the median ± k MAD band stays at the run's k over all samples, said in the live line | condition, batch, isolation 18,682 B | small | Median useful fraction runs from 0.845 with CD14 beads to 0.776 after FACS, so a low bubble may describe its group's protocol. |
| Where every read goes | Flows for: all samples (default) or one level of cell type, condition or batch, links precomputed per level | links, lead "Of {total} fragments, {useful} ({share}) were assigned to a transcript", caption n | about 9.4 kB est., in place of 78,099 B of per-sample fates | small | The per-sample median decoy share is 0.260 in BM_4 against 0.150 in LDN, so where reads go depends on the cells. |
| Absolute bands, run map, teaching subset | Metrics the registry keeps: shown (default) or hidden | rows; lead from 11 bad and 13 usable of 72 cells to 3 and 4 of 48; reading line naming each kept metric's registry cause | one registry action per metric row, under 0.3 kB est. | small | Per-base sequence content fails on all 8 runs from random-hexamer priming, so a red row the registry keeps is not a sample problem. |
| Run strips, 7 metrics, teaching subset | The same judge panel; strata cell type (default), none, batch, cell type and batch | as the sample strips | batch per row under 0.5 kB est. | small with the judge | With cell type and batch every stratum holds 2 runs, under the 3 a band needs, so 0 flags show that no flag is not a pass. |
| Base quality along the read | Deferred: low value at 8 runs; redesign at P2Q.8's 210 runs, about 1.24 MB of rows by arithmetic, as column arrays | | | | |
| Reads fastp keeps, by its quality and length settings | Deferred: fix the bands first, step 0; the -q and -l click-to-set is Phase 2D's planned contour helper; per-run grids wait for P2Q.8, wall time unmeasured. Isoline knob rejected as a workaround for the band defect | | | | |

### M3 trim, align, count: made-up preview, input to M3.6's figure contracts

Values on M3, M4 and M7 to M10 are planted by the made-up data; each teaching value is what the control shows if real data behave alike, to be re-measured before it justifies a build.

| Figure | Control: default and range | Recomputes | Data added | Effort | Teaching value |
|---|---|---|---|---|---|
| Where each sample's reads align | Show only: cell type, condition, read-length batch, All by default; Multimapper limit on STAR pages: as run (default) or 1, derived by M3.4's rule | bars; live line with shown of total, median uniquely mapped % and the lowest sample, through the `measure` series because a 100 % stack's total is always 100 %; caption n and median input reads | group as 3 columns, +5,940 B; at 210 runs about 250,400 B of rows and 127,000 B of traces est. | small | Filtering by batch shows whether the unique share follows read length, which goes with condition on Step-1. |
| One sample under each alignment setting | Show: share of input reads (default) or change from as run, in percentage points; no sample picker, since the page cannot realign | bars as differences per outcome, reading line naming the unit | none | small | A 1.4-point untrimmed effect is invisible on a 0 to 100 % stack and plain as a difference. |
| Reads assigned at each strandedness setting | Show: reads at each -s (default) or forward share, the registry's two thresholds as lines; Salmon's strand bias as a third choice with `--compare-dataset`. Editable thresholds low priority: no preview call changes below 0.92 | bars coloured by call | Salmon's share, under 2 kB at 33 samples est. | small | The call reads the forward share against two thresholds, not the most reads. |
| Gene counts with and without multimappers | Genes: all (default), protein-coding, lncRNA, other. Per-sample select deferred, it cannot be carried at 210 runs | points shown; n, identical, differ, Spearman and Pearson with CIs from the precomputed `ctx.concordance` row of that type | `gene_type` +31,000 B on 1,000 rows; stats rows about 11 kB est. | medium: per-type rows in the module | A near-1 ρ over all genes can hide one gene type that disagrees, expected for lncRNAs with antisense neighbours on M3.9's real counts. |
| How far the two counts differ, gene by gene | Limits: bias ± 1.96 SD (default) or empirical 2.5 and 97.5 percentiles; Genes as the previous figure | limit lines; lead "{n_out} of {n} genes lie beyond the limits", looked up | under 1 kB est. | small | Multimappers only add reads, so the difference is never below 0: the parametric lower limit is -0.255 where the empirical one is 0.0. |

### M4 quantification: made-up preview, input to the module's contracts

| Figure | Control: default and range | Recomputes | Data added | Effort | Teaching value |
|---|---|---|---|---|---|
| Salmon mapping rate by read-length batch | Slots: read-length batch (default), cell type, condition; colour the other main factor | slot medians in the page; test looked up per grouping, Mann-Whitney U with rank-biserial r for two levels, Kruskal-Wallis with ε² for more, BH across the groupings offered. Index switch not proposed: one Salmon run per sample is planned | under 1 kB est. | small | A mapping-rate difference by condition is read against the same samples grouped by read length. |
| Genes detected against depth | Depth: steps along one absolute depth grid; Detected at: 1, 10 (default), 100 reads | medians in the page; Mann-Whitney looked up per depth and threshold; caption says how many samples reach that depth | about 146 kB at 33 samples, about 1.4 MB at 210 est.; the module must subsample every sample on one absolute grid | medium | Detection compares fairly only at a common depth, and the 10-read rule moves every count. |
| Spearman correlation between samples | Group samples by: cell type (default), batch, condition | order, within and between medians in the page, permutation p looked up | under 1 kB est.; at 210 samples column arrays first, about 7.49 MB as rows est. | medium | Which label blocks the matrix shows whether cell type or read length dominates similarity. |
| Estimated reads against TPM in one sample | Gene counts as: Salmon estimates, scaled TPM, length-scaled TPM, default the module's `--counts-from-abundance` | y values, lead's long-over-short ratio of reads per TPM | +90,008 B at 2 decimals, measured with stand-in values | small | It shows what pydeseq2 receives under each choice. |
| Genes detected by Salmon, STAR and featureCounts | Detected at a mean of: 1, 10, 100 reads, default the page's; Genes: all (default), protein-coding, lncRNA, pseudogene | six counts and the caption's shared-of-total, looked up | 72 rows, about 8 kB est. | small | The tools disagree near the threshold, where Salmon-only genes have a median of 13.3 reads. |

### M5 differential expression

| Figure | Control: default and range | Recomputes | Data added | Effort | Teaching value |
|---|---|---|---|---|---|
| Samples: tissue PB, coloured by cell type | Colour by: cell type (default), condition, batch | colours, legend; lead "The {factor} accounts for {r2} of the variation between samples"; p for cell type only, the live line saying R2 is descriptive for the others because only the first factor was permuted within individuals | R2 and p per factor, a few hundred B est. | small | Samples group by cell type, R2 0.684, against 0.089 for condition and 0.026 for batch. |
| Samples: tissue PB; cell type NDN, coloured by condition | Colour by: condition (default) or batch, the other as symbol | colours, symbols, lead and test sentence; batch tested within condition, question 9 | batch R2 and p, under 100 B est. | small | Every pdac and ipmn sample is in the > 76 nt batch, the confound CONVENTIONS requires taught. |
| Samples: tissue PB; cell type NDN; condition steady_state, coloured by batch | Deferred: a components select teaches only if PC3 shows something, unmeasured; otherwise none, one factor and one sample per donor | | PC3 for 17 samples, about 0.4 kB est. | | |
| Volcano: gcsf over steady_state | None new: the calls control ships. Caption defects in step 0; contrast choice is question 8 | | none | none | |
| Fold change against expression: gcsf over steady_state | Deferred: thin by gene first, 1,645 of 5,435 genes carry one of their two estimates; then "Mark the top N by abs log2 FC of raw or shrunk", which moves 4 of 50 genes on ndn_gcsf and needs a small contrast beside it | | top lists about 9 kB est.; a second contrast about 0.84 MB est. | | |
| p-value histogram: ndn_gcsf | Contrast: ndn_gcsf (default) or each of the other nine | bars, flat line, lead, three reading lines, caption, all looked up | about 1.8 kB est. | small | First bar over flat is 1.0 for the null split, 7.9 naive against 1.7 adjusted, 31.7 paired against 20.1 unpaired, and the words say adjusting for a batch nested in condition removes signal as well as noise. |
| Dispersion: ndn_gcsf | Deferred: a design select on the ndn_pdac pair only if the gene-wise cloud moves, unmeasured; thin the trend line regardless, it carries all 12,796 points in 481,483 B | | about 0.40 MB per design est.; thinning frees about 474 kB est. | | |
| Genes called: batch lesson | Significance: BH (default), BY, Holm, Bonferroni, raw p; below 0.1, 0.05 (default), 0.01, 0.001; fold change at least 0, 0.5, 1 (default), 1.5, 2 log2; reads the `calls` state when it is a grid setting | bars, lead "{contrast} calls the most genes: {n}", the caption's rule; lookups in a grid built by `ctx.de.stats.classify`, once, with M5's planned two-knob map | 14,658 B for both lessons | small | Under BH the naive design calls the most at all 20 settings, but on raw p with a 1.5 log2 floor the within-batch contrast leads, 441 against 376 and 327, so the naive excess sits in small fold changes. |
| ndn_pdac_naive against ndn_pdac_adjusted | Colour by: called in either (default) or agreement of the calls, 5 levels | colours, live counts, ρ over genes called in either, precomputed; the caption's ρ over all genes stays | a few bytes | small | Of the naive design's 659 calls, 27 survive adding batch and 13 appear only after it. |
| Genes called: pairing lesson | As the batch lesson | as the batch lesson | in the 14,658 B | small | ndn_gcsf_paired calls more than ndn_gcsf_unpaired9 at all 40 BH and raw-p settings, 4,018 against 3,457 at the default. |
| ndn_hsct_early against mono_hsct_early | As naive against adjusted | as naive against adjusted | a few bytes | small | 1,002 genes are called only in neutrophils, 81 only in monocytes and 28 in both the same way, ρ 0.460 over the 1,113 called in either. |
| Named genes: hsct against steady_state | p above each gene: Mann-Whitney, BH over the 7 shown (default), or DESeq2 padj, BH over all 14,668 tested. Gene ticks rejected: choosing a family after seeing the data is the pitfall. Contrast option deferred, it touches R4 and D-046 | p labels, the caption's family, Pitfalls line | none | small | 6 of 7 named genes pass 0.05 in the 7-gene family and 3 of 7 under DESeq2's padj, STAT1 0.037 against 0.12. |

### M6 function

| Figure | Control: default and range | Recomputes | Data added | Effort | Teaching value |
|---|---|---|---|---|---|
| Over-representation: gcsf over steady_state | Sets shown: passing sets only (new default, ends 5 zero-length down bars) or 7 smallest FDR per direction; Contrast: ndn_gcsf (default) or the other nine at default cut-offs, the `contrast` state shared with GSEA; FDR is P2S.5, cited; DE cut-offs padj, abs log2 FC and adjustment only through GS.5's engine, reading M5's `calls` state | bars, lead "N gene sets pass FDR < 0.05 on the up or the down calls", live call counts | contrasts 243,719 B, 290,755 B with gene lists; cut-offs 375,580 B featured, 3,406,267 B all ten | small for sets and contrast; large for cut-offs until GS.5, medium after | ORA's answer moves with a cut-off chosen before the test, 17 to 24 passing sets on ndn_gcsf over padj 0.01 to 0.1 and abs log2 FC 0 to 2. |
| Gene set enrichment: gcsf over steady_state | Contrast: ndn_gcsf (default) or nine others; Sets shown: 15 smallest FDR (default) or every tested set, 48 or 49 | bars, title, lead "N of M tested gene sets pass FDR < 0.05", caption's contrast and count; nothing statistical in the page | 483 rows, 123,151 B, one copy shared with the next figure, replacing 3,348 B | small | Passing sets run from 36 of 49 in mono_gcsf to 2 of 48 in the null split, whose 2 are the interferon sets at FDR 0. |
| HALLMARK_INTERFERON_ALPHA_RESPONSE in every contrast | Gene set: the dataset's featured set (default) or any of the 50 Hallmark sets; FDR is P2S.5, cited | bars, title, lead, caption; a set not tested in a contrast named in the live line | 0 B beyond the shared rows | small | Both interferon sets pass even in the null split while most sets do not, so a pathway that passes everywhere is weak evidence for one contrast. |
| Running score: HALLMARK_INTERFERON_ALPHA_RESPONSE in ndn_hsct_early | Gene set: the dataset's running set (default) or the 49 tested in that contrast; Contrast later, under the shared state | running sum by `lite._running`'s rule, ES, member ticks in the page; NES, p, FDR looked up; the page checks its ES against the table | +322,588 B and -405,330 B for the dropped line, net -82,742 B; all ten contrasts 2,218,605 B | medium: the running-sum twin, matching the table within 7.3e-7 on 145 sets in numpy, not yet in JS | The student sees what an enrichment score is for any pathway: an early peak, a late trough or a flat wander. |
| The background decides the p-value | Background: tested genes 12,768 (default), gene-set union 14,477, protein-coding genome 23,315, every annotated gene 77,311; Calls: up (default) or down; Bars: 7 smallest p against tested genes (default) or all 48 | bars; lead "Against {N} genes instead of {tested}, {more} of {shown} sets look more significant and {less} less", passing count, all looked up | 396 rows, 81,536 B; GS.5's live engine later needs about 1 kB est. | small | Sets passing on ndn_gcsf's up calls are 16 against tested genes, 6 against the union, 19 against protein-coding and 44 against every annotated gene, so a wrong background moves p either way. |

### M7 RNA versus protein: made-up preview, input to the module's contracts

| Figure | Control: default and range | Recomputes | Data added | Effort | Teaching value |
|---|---|---|---|---|---|
| Protein groups quantified per sample | Minimum share quantified: 50, 60, 70 (default, CONVENTIONS), 80, 90, 100 % | bars of protein groups kept, lead "{n_kept} of {n_ids} protein groups"; live line says the caption's counts do not move | about 21 kB est. | small | The 70 % rule trades proteins tested against samples per protein, and the plex sets the missingness. |
| Per-gene RNA-protein correlation, by gene class | Correlation: Spearman (default) or Pearson; Minimum share quantified, as above | ECDFs and class medians in the page; Kruskal-Wallis H, p, ε² and the lead's BH count looked up per setting | Pearson ρ and share per gene, +158,948 B, measured with stand-in floats | small | Raising the share shows whether better-measured proteins track their RNA more closely. |
| RNA against protein fold change: case vs control | abs log2 FC at least 0.5 (default), BH-adjusted p below 0.05 (default), optionally a separate protein threshold | colours, counts of both, RNA only, protein only, the reading line's rule; live line says ρ and the quadrant counts use every gene | none; controls about 1.5 kB est. | medium: calls with two layers and its twin; thin by the loosest setting at CPTAC sizes | TMT ratio compression drops protein calls first as the threshold rises. |
| Paired RNA and protein, sample by sample | Genes: changed at both (default), RNA only, protein only; Number: 6 or 12 (default, today's). Scatter-to-heatmap link rejected, about 13.2 MB as rows est. | rows swapped for the precomputed top genes, lead naming the class | about 226 kB est. | small | Protein-only genes show protein rows moving under flat RNA rows, the case the page names but does not draw. |

### M8 non-coding: made-up preview, and planned figures

| Figure | Control: default and range | Recomputes | Data added | Effort | Teaching value |
|---|---|---|---|---|---|
| Figure 1. Counted reads by biotype, per sample | Filters biotype, cell type, condition, All by default; the median line hidden while biotype is All | bars, median share of the chosen biotype over the samples shown, shown of total; live line says the lead's medians cover all samples | none; controls a few hundred B est. | small | One biotype's share per cell type or condition is readable against the axis, where a stack hides it. |
| Figure 2. Expression by biotype, against the filter | Minimum count per biotype: 10 (default, `MIN_COUNT`), presets 1, 2, 5, 20, 50; the rule shared with M5's planned min-count grid | dashed line per biotype, kept shares, the reading line and the caption's cut-off; live line says Figures 3 and 5 keep the default | the total-count rule per gene, 56,000 to 59,500 B by arithmetic, without which every count below 10 is an upper bound | medium | Lowering the lncRNA cut-off shows how many lncRNAs the default keeps out of every test: 344 kept at 10, at most 1,065 at 1. |
| Figure 3. case against control in type_a, by biotype | Rytis's control, TODO.md line 617: the M5 volcano's cut-offs and live counts, biotype colours for called genes, BH within each biotype; the pooled basis dropped | colours, fold-change lines, the lead's called counts per biotype, the reading line; a raw-p line per biotype or none | none; the page table cut to the calls columns, 235,710 B against 402,606 B today | medium: calls with a colour column | The same rule on a stratified BH shows which lncRNA calls a student would report at each cut-off. |
| Figure 4. LNC-X against its neighbour GENE-Y | Rytis's control, TODO.md line 618: a drop-down of lncRNA and neighbour pairs named by both genes, the planted pair first; no subsets of 12 samples | the pair's points, Spearman on values and on residuals in the page, p per pair looked up | about 6.9 kB per pair est., one pair measured at 6,899 B | medium: the Spearman twin | Residuals show how much of a pair's correlation comes from the shared design: 0.717 on values, 0.356 on residuals for the planted pair. |
| Figure 5. Fold changes of miR-X targets against non-targets | Targets: validated, predicted, both apart (default); context-score cut-off at a few precomputed values | curves, KS D and p, BH over the tests shown, the lead's median shifts, caption n, all looked up | about 2,720 B est. on the preview; real n_predicted × about 20 B | medium: computation in the module | A stricter prediction cut-off trades target count against shift size. |
| Planned, no page: count distributions per biotype, violin | Shares Figure 2's min-count state, or Figure 2 covers it | | same values | contract | |
| Planned, no page: miRNA volcano | M5's calls kind unchanged | | about 203 B per miRNA est. | contract | |
| Planned, no page: miRNA-target anticorrelation | Pair select from a precomputed table, ρ live, p per pair looked up | | samples × pairs × 2 values est. | contract | |
| Planned, no page: miRNA to mRNA to protein tri-layer | miRNA select among those precomputed; needs M7 | | pending M7 | contract | |

### M9 gene card: made-up preview, and planned figures

| Figure | Control: default and range | Recomputes | Data added | Effort | Teaching value |
|---|---|---|---|---|---|
| Figure 1. MADEUP1 in this dataset, by cell type | Batch: all (default), Reads ≤ 76 nt, Reads > 76 nt, three options computed by `annotate.compare_groups` | boxes, n, U, p, BH over cell types holding both groups, r, caption, lead; a cell type with one group shows points and no test | a few kB est. | small | Testing within one batch shows whether a gene's shift survives Step-1's read-length confound. |
| Figure 2. MADEUP1 across normal tissues | Donors: all (default), female, male | each tissue's median and IQR as polars' linear quantiles, the caption's samples, the lead's highest and lowest tissue | none on the preview; real GTEx per-sample sex is unverified, for the external-data spike | medium: the quantile aggregate | A sex-biased gene's tissue medians move with the donor filter, which a pooled median hides. |
| Figure 3. MADEUP1 across immune cell types | Lineages included: six checkboxes, all on; tau on: log2(nTPM + 1) (default) or nTPM | bars, tau in lead and caption, number of cell types, top cell type | none | medium: the tau aggregate | tau moves from 0.43 to 0.81 with the scale and to 0.30 without T cells, so a specificity score is read with its panel and transform. |
| Figure 4. Exon use of MADEUP1 across tissues | Tissues included: checkboxes, all on | each exon's median share, log2 ratios, levels, row CV, caption count, lead, and the qPCR pair key number, or the live line says it does not follow | none | medium | Removing the degraded liver moves the qPCR pair from E5-E6 to E6-E7, so uneven coverage can be 3' bias rather than splicing. |
| Planned, no page: isoform usage, GTEx isoform stacks | Unit: TPM or share of gene; tissue subset as Figure 4 | | pending the spike | contract | |
| Planned, no page: GTEx junction arcs | Tissue select; minimum read count hiding weak junctions | | pending the spike | contract | |
| Planned, no page: HPA tissue RNA against protein, single-cell, cell lines, cancer IHC | Antibody reliability filter; Figure 3's tau view; lineage filter; antibody select where HPA lists several | | pending the spike | contract | |
| Planned, no page: cross-cohort row, clinical row, evidence view | Cross-cohort: D-012's call thresholds per cohort, as calls does; clinical row and evidence view are GS.7 and GS.9, cited | | small est. | contract | |
| Planned, no page: subcellular schematic, blood, structure | None: categorical locations, a few blood values, a linked structure; no knob in the method | | | none | |

### M10 clinical: made-up preview, and planned figures

| Figure | Control: default and range | Recomputes | Data added | Effort | Teaching value |
|---|---|---|---|---|---|
| Figure 1. Kaplan-Meier estimate of PFI, split at the median of G00001 | Planned, cited: TODO M10 Interactivity and GS.7. Requirement from this page: z carried exactly or as ranks, so a percentile split cannot move a patient | | cohort 16,895 B at 4 decimals; about 43,800 B est. for TCGA-BRCA | large, in GS.7 | |
| Figure 2. Log-rank p of every cut-off from the 20th to the 80th percentile of G00001 | Planned, cited: the cut-point lesson and GS.7. Requirement: the corrected line recomputed per scanned range, or the live line says it holds for 20 to 80 only | | permutation matrix about 427,000 B est., or the cohort | large, in GS.7 | |
| Figure 3. Power of a Cox model by number of events and hazard ratio per SD | alpha 0.05 (default), 0.01, 0.001 with the HPA source written beside it; target power 80 % (default) or 90 %; marked HR adjusted (default) or univariable | grid, bands, isoline, mark, the lead's power, the caption's alpha | none, closed form | medium: the Hsieh-Lavori twin | The same cohort is powered at 97.5 % at 0.05 and 73.5 % at 0.001, so a portal's threshold reads in events: 81 against 176 for 80 % power. |
| Figure 4. Cramér's V between the variables of the sample sheet | Bias correction: none (default) or Bergsma 2013, citation to verify; Samples: all (default) or within each cell type; 4 states precomputed | every V and level, the lead's strongest pair, caption n and correction words; a pair constant within a stratum blank, said in the reading | about 7 kB est. | small | On 33 samples uncorrected V puts noise in the small bins: condition against batch falls from 0.149 to 0 with the correction. |
| Planned, no page: expression by clinical variable | Grouping select; each grouping's test block precomputed by `ctx.viz.annotate` | | about 98,000 B est. | contract | |
| Planned, no page: genome-wide contrasts and Cox screen | Calls kind plus a contrast select within the cap; Cox screen precomputed per D-049 | | about 203 B per gene per volcano est. | contract | |
| Planned, no page: cumulative incidence for competing risks | Event type and endpoint select over precomputed curves | | small est. | contract | |
| Planned, no page: swimmer plot, KM with CI bands, ER exercise | Cited: P2Q.4 linked view, M10 Interactivity, GS.7 | | the cohort table | contract | |

## Build order

### Step 0: defects, fixed at the root before or regardless of any control

| When | Defect | Fix |
|---|---|---|
| now | The bar filter never refills the lead or caption: with NDN the plot shows 73 of 210 runs and a median of 26.2 M while the lead keeps "135 of 210 runs have fewer than 30 million", on 6 pages | the words contract in the `filter` kind, step 3 |
| now | The made-up M1 page prints one median as 0.5 and 0.494 | one rounding |
| now | M2 size against genes points at the dashed line while 102 of 107 flags are on neither plotted metric | default the triangles to the plotted metrics |
| now | M2 fastp grid: 8 equal bands over 0 to 100 % put the 30 cells that keep reads in the top two | band edges over the cells that keep reads, -l above the read length its own "keeps none" level |
| now | [qc_metrics.toml](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/conf/qc_metrics.toml) line 149 says fastp's settings are ones "the page lets you move" | correct the sentence until Phase 2D's click-to-set exists |
| now | M3 preview caption omits the empirical percentiles its plan gives | caption, and M3.6's contract |
| with P2S.4 | M5 volcano caption says "all drawn" then "draws 1,842 of the 5,526", and prints "every 3th", from English written in [scatter.py](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/viz/scatter.py) lines 355 to 359 | add scatter.py to P2S.4's file list; move the sentence to content; fix the ordinal |
| with P2S.4 | M5 MA plot thinning splits gene pairs; M5 dispersion trend carries 12,796 points | thin by gene; thin the trend, which frees the M5 headroom |
| with P2S.4 | M6: zero-length down bars; GSEA reading line 3 points at nothing; the background reading line and UNION_CHECK say the opposite of the page's data; fdr and rank hovers unformatted; "1681" with no separator, 2 short of M5's 1,683 up calls, gap unexplained | passing sets by default; rewrite the words; hover formats; one number format |

### Steps 1 to 8

| Step | Module and figures | The shared page script gains | Waits for | Export dependency |
|---|---|---|---|---|
| 1 | qc-page-fixes merges with the export | `ctxExport` | questions 1 and 2; the M5 trend thinning | it is the export; it must read the caption the page shows |
| 2 | P2S.4, then P2S.5 as one `threshold` rule on M6's carried FDR | the FDR rule | the Step-1 record | P2S.5 is the export's first BarData consumer beyond the filter |
| 3 | M1: Figure 2's words and depth presets, then Figure 1's contrast, colour and count | words contract, page state, `subset` from `filter`, `view`, `threshold` with its first rule | step 1 | BarData only, covered |
| 4 | M2: PCA detail select; run band map; fates, flow, size against genes; the judge for sample strips, run strips and Distance; the band map's flag view; Agreement peers; the sex view after the record | brush detail select; the `checkpoint.judge` twin | step 3; the sex view waits for `ctx.qc.identity` | StripData, HeatmapData, DensityData and FlowData need `controls()` and the export hook; svg2pdf on a Sankey and a heatmap is unverified |
| 5 | M5: F6, F1, F2, F9, F11, F12, then F8 and F10 on the grid shared with the two-knob map; per-contrast volcano pages if question 8 says so | the shared `calls` state | P2S.4; the trend thinning; GS.3 for per-contrast pages | BoxData for F12 |
| 6 | M6: GSEA and every-contrast views on one row copy, background, ORA sets and contrast; the running score; ORA cut-offs | the running-sum twin; `contrast` and `gene_set` state | P2S.5; GS.5 for the cut-offs, which also turns the background into a live lookup | none beyond BarData and ScatterData |
| 7 | M3: M3.6 builds its five figures with controls from the start; M8: F1, F2, F4, F5 on the preview, F3 after the calls rework | calls with a colour column; the Spearman twin; the min-count rule shared with M5 | M3.6; for M8 F3, GS.3 and GS.4, and the qc-page-fixes merge, Rytis's condition at TODO.md line 617 | CurveData for M8 F2 and F5 |
| 8 | M9, M10, M4, M7 as each module is built | quantile and tau aggregates; the power twin; GS.7's survival twins | the external-data spike; GS.7; GS.9; each module's build | BoxData, SurfaceData, CurveData, HeatmapData |

If question 3 is answered yes, step 3 and the non-sex part of step 4 move before step 2, and P2S.5 reuses `threshold`.

## Questions for Rytis

1. **How does a changed view export?** Recommended: accept route A as built in the qc-page-fixes worktree for changed views only, keeping the Matplotlib files for each default view, on three conditions: it reads the caption the page shows, not the build's; it extends to every figure contract a control lands on; its live line, appended to the caption, says the file is the page's drawing at the figure's size. It departs from "PDF and SVG from one render" for changed views and writes the white theme only. The alternatives are prerendered files per option, exact but for selects of up to about 10 options and about 0.84 MB per figure on disk by the critique's arithmetic, or the control writing the Reproduce params for `ctx figure`, exact but after Phase 2D.
2. **Who pays its 1,607,890 B?** Recommended: subset the two Liberation Sans files to the characters the pages set, then measure, before the merge; the subset's size is not measured. Thin the M5 trend in the same change, since the M5 overview reaches 14,387,080 B at the merge by arithmetic. The alternative is to export only from overview pages, with standalone pages linking the default files; across up to 124 pages the full assets come to about 199 MB by arithmetic.
3. **Build the kinds on M1 and M2 before P2S.4 and P2S.5?** Recommended: yes, right after the merge. Neither module waits for the record except the sex view, both keep their text in `content/gather.toml` and `content/qc.toml` today, and P2S.5 then becomes one `threshold` rule rather than a fifth control. The alternative is the order TODO.md's P2S.6 sets, M1 and M2 after P2S.4 and P2S.5.
4. **Where do the words a control cannot keep true go?** Recommended: in the live line, which the export appends to the caption; the caption keeps the default view's statistics, and a `view` carries one caption per option. The alternative, refilling spans inside the caption, must protect D-048 marks and the 40-word budget, which 48 of 108 rendered captions already exceed by the critique's count.
5. **Presets or typed values?** Recommended: presets wherever a source names the values: M1 depth 20 and 30 M from `total_sequences`; M10 alpha 0.001 from HPA, with the source beside it; M7's 70 % from CONVENTIONS; M2's k steps around the run's `--mad`. A typed value stays only where a page takes one today, the volcano.
6. **M1 Figure 1: a Contrast select or free filters?** Recommended: a Contrast select applying each contrast's declared filters, with Colour by batch, isolation or recorded sex and Count runs or individuals. Tissue and time point are not offered as colours: they nest in cell type and condition, with 20 and 16 of 22 groups in one level only. The alternative is the proposer's three Show-only selects, which build groups no contrast tests.
7. **M2 sample band map: which default?** Recommended: the flag map, every judged metric against its stratum's median ± k MAD with the absolute band as one row, retitled for it; the absolute-band view becomes the other choice. On Step-1 today's figure is one row of 210 good cells.
8. **M5's planned contrast selector for the volcano?** Recommended: one standalone volcano page per contrast after GS.3, each with the calls control, the overview keeping the featured contrast. On the overview a contrast costs about 2.9 MB today and about 0.7 MB after GS.3 by derivation, so ten come to about 16.8 MB by arithmetic, over the cap. The alternative is a selector over the lesson contrasts only, measured after GS.3.
9. **M5 Figure 2: how is batch tested on the NDN PCA?** Recommended: PERMANOVA of batch within condition, permutations restricted to the condition levels holding both batches, the live line naming them. The alternative is batch after condition, sequentially. A marginal test would absorb condition, because every pdac and ipmn sample is in the > 76 nt batch.
10. **Build controls on made-up previews?** Recommended: only the two M8 controls you named; write every other M3, M4, M7 to M10 control into its module's figure contract in TODO.md, and re-measure its teaching value on real data before building it. Preview values are planted, and D-050 removes those pages at the release.
11. **M6's planned ranking-metric knob, which no group proposed?** Recommended: measure the wall time of one extra ranking at 1000 permutations on the featured contrast, then decide; if built, a precomputed view with the Wald statistic, the CONVENTIONS default, first.
12. **M1's planned minimum group size?** Recommended: no control until a source gives its default. On Step-1, groups below 2, 3 and 5 number 1, 3 and 5 by runs and 1, 3 and 7 by individuals.

## Not done

- I did not run pixi, build a page, or test any control in a browser.
- My own measurements are the export's bytes, from the vendored files with `page_assets()`'s steps repeated in Python, the kit's page count and size, and the cited lines of `save.py`, `qc_metrics.toml` and `scatter.py`. Every other number is the inventory's or the critique's, as marked.
- I did not re-run the proposers' ports: M2 rejudge and peers, the M6 lite copy and running sum, the M5 lesson grid.
- PC3, any dispersion table beyond ndn_gcsf, the subset fonts' size, and the 44-of-48 background claim beyond the 7 sets on the page remain unmeasured.
- The Distance peers row departs from both the proposer and the critique: distances are precomputed per choice rather than recomputed in the page, so no JS twin of `pca_distance` exists. Its 13 kB is my estimate at about 30 B per value.
