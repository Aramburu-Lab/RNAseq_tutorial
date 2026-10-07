# M2 background bands: proposal for approval

**27 metrics are in scope. 11 get a verified sourced cut-off, 4 get a proposed course convention, and 6 keep the dataset-only check with no band. For the other 6, good, usable and bad do not apply, so they get a reference line or call regions instead. Of the 15 metrics with a band, 10 would only shade and 5 would judge. Four of the five are FastQC verdicts that already judge today; `salmon_percent_mapped` is the one new judging band. Under every recommended row, 0 Step-1 flags change.**

Shading will not show on any M2 strip until the code in decisions 1 to 4 exists. Today `strip()` draws only `good` and `usable`, and the registry's validator refuses 7 of the proposed pairs.

Every count below was recounted on 2026-10-01 by [recount.py](/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/bands/proposal_check/recount.py). [apply_check.py](/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/bands/proposal_check/apply_check.py) runs the registry's own validator on [proposal.json](/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/bands/proposal.json) applied to an in-memory copy of `conf/qc_metrics.toml`. The counts come from two files:

- **M** is [2026-09-30_qc_tier1_metrics_v1.csv](/Users/rytis/ctx_tier1/outputs/02_qc/2026-09-30_qc_tier1_metrics_v1.csv). Its unit is the sample: 210 samples, one run each.
- **T** is [2026-09-30_qc_tier1_teach_read_metrics_v1.csv](/Users/rytis/ctx_tier1/outputs/02_qc/2026-09-30_qc_tier1_teach_read_metrics_v1.csv). Its unit is the run: the 8 teaching runs.
- `2026-09-30_qc_tier1_teach_metrics_v1.csv` holds only `sample`, so the subset's Trim Galore, Salmon and count values are read from its 8 rows in M.

S1 to S16 refer to the quoted sources at the end.

| Metric | Title | Unit | Good / usable, direction | Mode | Source, source type | Step-1 good / usable / bad / missing | Flags that change |
|---|---|---|---|---|---|---|---|
| total_sequences | Reads | reads | 30 M / 20 M aligned reads, higher; unchanged plot_bands | shade, as today | ENCODE June 2017 page and audit code, consortium_standard, S1 | T `fastqc-total_sequences`: 0 / 0 / 8 / 0. The course cut every run to 1,000,000 reads. M has no column | 0. Judging would add 8 on T, flagging all 8 runs instead of 4, every one because of subsampling |
| read_length | Read length | nt | does not apply; one line at 50 | shade, line only | ENCODE minimum read length, consortium_standard, S2 | T `fastqc-avg_sequence_length`: 8 above, 0 below, 0 missing. M has no column; the Trim Galore JSONs give 210 above | 0 |
| percent_gc | GC content | % | does not apply; one line at 49.8 | shade, line only | course convention: the mean of MultiQC 1.35's hg38 transcriptome GC curve, derived, S3 | T `fastqc-percent_gc`: 7 below, 1 above, 0 missing. M has no column | 0; the MAD flag on ERR7425270 stays |
| fastp_surviving | Reads passing fastp | % | 95 / 90, higher | shade, at default settings only | course convention: at its defaults, fastp drops a read only when more than 40 % of bases are below Q15, it has more than 5 N, or it is under 15 nt. 95 is about 5 times the worst loss among the 8 runs, 1.03 %. S4 | T `fastp-pct_surviving_reads`: 8 / 0 / 0 / 0. M has no column | 0 |
| fastqc_per_base_quality | Per-base quality | status | no change: the verdict uses lower quartile 10 / 5 and median 25 / 20; the shading uses 28 / 20, higher | judge by FastQC status and shade, both as today | FastQC 0.12.1 limits.txt and QualityBoxPlot.java, tool_code, S5 | T `per_base_sequence_quality`: 8 / 0 / 0 / 0 | 0 |
| fastqc_adapter | Adapter content | status | 5 / 10, lower | judge by status as today; add shading | FastQC 0.12.1 limits.txt, tool_code, S6 | T `adapter_content`: 5 / 1 / 2 / 0 | 0; ERR7425288 and ERR7425289 are already flagged |
| fastqc_n_content | N bases | status | 5 / 20, lower | judge by status as today; add shading | FastQC 0.12.1 limits.txt, tool_code, S6 | T `per_base_n_content`: 8 / 0 / 0 / 0 | 0 |
| fastqc_sequence_content | Per-base sequence content | status | 10 / 20, lower; positions 1 to 13 grey | shade only; action stays keep | FastQC 0.12.1 limits.txt, tool_code, S6; the 13 positions from Hansen et al. 2010, S7 | T `per_base_sequence_content`: 0 / 0 / 8 / 0 by status. Counting from position 14: 7 / 1 / 0 / 0 | 0. Setting the action to flag would add 8 |
| fastqc_gc_distribution | GC distribution | status | 15 / 30, lower | judge by status as today; add shading | FastQC 0.12.1 limits.txt, tool_code, S6 | T `per_sequence_gc_content`: 4 / 3 / 1 / 0 | 0; ERR7425271 is already flagged |
| fastqc_overrepresented | Overrepresented sequences | status | 0.1 / 1, lower | shade only; action stays keep | FastQC 0.12.1 limits.txt, tool_code, S6 | T `overrepresented_sequences`: 6 / 2 / 0 / 0 | 0 |
| tg_pct_with_adapter | Reads with adapter | % | does not apply as computed; one line at 33.3, the share expected by chance | shade, line only; the MAD check stays | course convention, derived from the per-length expectation in Trim Galore 2.3.0's report.rs, S8 | M `trim_galore-tg_pct_with_adapter`: 2 below the line, 208 above, 0 missing. FastQC's 5 / 10 on this column would make 210 of 210 bad | 0. Decision 6 offers a measure that can carry a band |
| tg_pct_passing | Reads kept by Trim Galore | % | 99.9 / 99, higher | shade | course convention: the example report in Trim Galore 2.3.0's length guide removes 552 of 1,000,000 reads, 0.1 %. S9 | M `trim_galore-tg_pct_passing`: 210 / 0 / 0 / 0. The subset: 8 good | 0; the 21 MAD flags stay |
| tg_pct_quality_trimmed | Bases trimmed for quality | % | 2.3 / 20, lower | shade | course convention. Good is Trim Galore's 2.3 % in Liao and Shi 2020, Table S1, S10. Usable is the 20 % of bases the q30_rate line allows below Q30, an analogy | M `trim_galore-tg_pct_quality_trimmed`: 210 / 0 / 0 / 0. The subset: 8 good | 0; the 1 MAD flag stays |
| gene_strand_bias | Strand of the gene-assigned reads | fraction | does not apply; call regions: reverse at or below 0.2, unstranded where the shares differ by at most 0.1, forward at or above 0.8 | shade regions; judge by the call as today | nf-core 3.27.0 thresholds, transferred; course_convention per D-039, S11 | not measurable: there is no `outputs/03_process`, and Step-1 has no genome alignment. Missing 210 of 210 and 8 of 8 | cannot be counted |
| salmon_percent_mapped | Reads Salmon assigned | % | 30 / 10, higher | judge | Salmon 2.7.0 diagnostics.rs, tool_code, S12 | M `salmon-percent_mapped`: 209 / 1 / 0 / 0; the usable one is ERX6994414 at 29.29. T has no Salmon column, so 8 are missing | 0; the 4 MAD flags stay |
| salmon_strand_bias | Strand of the reads | fraction | does not apply; nf-core call regions, plus Salmon's 0.49 to 0.51 strip for runs quantified as U | shade; judge by the call as today | nf-core 3.27.0, pipeline_default, S11; Salmon 2.7.0 libtype.rs, tool_code, S13 | M `salmon-strand_mapping_bias`: 210 unstranded. 209 inside the strip, 1 outside: ERX6994473 at 0.5138 | 0. Judging on the strip would add 1 flag row and flag 0 new samples |
| salmon_compatible_ratio | Reads compatible with the library type | %; the registry says fraction | one line at 95, higher | shade, line only | Salmon 2.7.0 libtype.rs, tool_code, S13 | M `salmon-compatible_fragment_ratio`: 210 at 100, 0 below, 0 missing. A single-end U run reads 100 by construction | 0 |
| library_size | Library size | reads | 30 M / 20 M aligned reads, higher | shade | ENCODE June 2017 page and audit code, consortium_standard, S1 | M `library_size`: 29 / 81 / 100 / 0. The subset's 8 rows: 2 / 4 / 2 / 0 | 0. Judging would add 95 flag rows and move 45 samples from pass to flag |
| genes_detected | Genes detected | genes | none | MAD only | none: ENCODE's 8,000 / 4,000 is for long-read RNA-seq | M `genes_detected`: no band; median 22,418 | 0 |
| top_gene_share | Share of the top gene | fraction | none; fix the definition first | MAD only | none sourced | M `top_gene_share`: no band; median 5.34 %. The 18S-antisense lncRNA ENSG00000280614 is the top gene in 102 of 210 samples | 0 |
| mt_fraction | Mitochondrial transcripts | fraction | none, no line | MAD only | none: the only number found, 5 %, is a per-cell single-cell filter | M `mt_fraction`: no band. A 5 % line would put 20 in bad, 18 of them monocytes | 0. Judging at 5 % would add 15 rows and flag 5 more monocytes |
| mt_rrna_fraction | Mitochondrial rRNA | fraction | one line at 0.1, lower | shade, line only | RNA-SeQC 2.4.2 report.py and GTF.cpp, tool_code, S14 | M `mt_rrna_fraction`: 204 below, 6 above, 0 missing; all 6 are already MAD-flagged | 0; judging would also add 0 |
| rrna_fraction | Nuclear rRNA | fraction | none; fix the definition first | MAD only | none: every rRNA cut-off found counts 18S and 28S reads, which this column cannot see | M `rrna_fraction`: no band; max 6.8e-5 | 0. Decision 10 would remove 11 |
| globin_fraction | Globin | fraction | 0.001 / 0.01, lower | shade | course convention placed on Step-1: every peripheral-blood sample is at or below 0.062 %. Shin et al. 2014 gives the scale only, S15 | M `globin_fraction`: 207 / 3 / 0 / 0; the usable 3 are ERX6994368, ERX6994432 and ERX6994441, all BM_4. The subset: 8 good | 0; the 53 MAD flags stay |
| replicate_spearman | Agreement with its group | rho | none on this column | MAD only | ENCODE's 0.9 / 0.8 is verified, S16, but ENCODE computes the statistic another way | M `replicate_spearman`: at 0.8, 67 / 0 / 142 / 1 | 0. Judging at 0.8 would add 133 rows and move 58 samples from pass to flag |
| pca_distance | Distance from its group | VST units | none | MAD only | none: the distance has no absolute scale | M `pca_distance`: no band; 1 missing | 0 |
| sex_check | Recorded sex against expression | category | does not apply; shade the two half-planes of the XIST against chrY scatter | shade the scatter only | course_convention, rule R8 | M `sex_agrees`: 208 agree, 2 not recorded, 0 disagree | 0 |

## What Rytis must decide

Each item gives the recommended answer first.

1. **Let `strip()` draw `plot_bands`, and a new single line, when a metric has no `good` or `usable`. Recommended: yes.** [figures.py](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/qc/figures.py) `strip()` draws only `Metric.spans()`, and those come from `good` and `usable`. Today only the per-base quality figure and the M1 depth figure read `plot_bands`. Without this change, 11 shade rows draw nothing on M2.
2. **Check the order of `plot_bands` against `direction`, and read the side from the pair itself where direction is `none`. Recommended: yes.** [registry.py](/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/src/ctx/qc/registry.py) lines 200 to 204 refuse any pair not sorted high to low. apply_check.py shows 7 rows refused today: the 5 FastQC pairs, tg_pct_quality_trimmed and globin_fraction. The alternative is a `plot_direction` field. It costs one more field, and in return it catches a swapped pair.
3. **Add a `plot_line` field for a single sourced reference value, sourced by `plot_source`. Recommended: yes.** It would carry five lines: read_length 50, percent_gc 49.8, tg_pct_with_adapter 33.3, salmon_compatible_ratio 95 and mt_rrna_fraction 0.1. The loader ignores unknown keys, so `plot_line`, `plot_expected_bias_positions` and `plot_unstranded_tolerance` load without error and draw nothing until code reads them. A check for unknown keys belongs in the same change. For mt_rrna_fraction and salmon_compatible_ratio there is an alternative: set `good` and `usable` both to the line. Today's code already draws and judges that, and it changes 0 Step-1 flags. But it turns RNA-SeQC's display line and Salmon's warning into a fail.
4. **Shade the FastQC status tabs behind per-position curves, as the per-base quality figure does, and leave the verdict alone. Recommended: yes.** Adapter, N and sequence content each need a curve by position. GC distribution needs a deviation strip, with the histogram drawn over MultiQC's hg38 transcriptome curve. Overrepresented sequences needs a bar of the top sequence's share. The grey positions 1 to 13 go in the new field `plot_expected_bias_positions = 13`, because QC thresholds live only in the registry.
5. **Shade total_sequences on the teaching tab, with a caption saying the course cut every run to 1,000,000 reads. Recommended.** All 8 runs land in ENCODE's bad region by construction. The alternative is to leave this one strip unshaded.
6. **For tg_pct_with_adapter, copy nf-core's FastQC on the raw reads for all 210 samples. Recommended.**
   - The Step-1 run executed 210 FastQC 0.12.1 tasks, but `fastqc/raw` was not copied into `data/raw/tier1_rnaseq`.
   - With it, FastQC's own adapter verdict judges the measurement its limits were written for. total_sequences, read_length and percent_gc also gain full-cohort columns.
   - Every FastQC status metric would then judge the full cohort too, so full-cohort flags would change. I cannot count them until the files are here.
   - `ctx gather` would need to read `fastqc/raw` through MultiQC.
   - Fallback: a new metric, the share of reads with at least 12 bp of adapter removed, read from MultiQC's `trim_galore_adapter_length_plot.txt`. It would take plot_bands [5, 10], lower, course_convention "FastQC 0.12.1 adapter warn 5, error 10, transferred". On M it gives 117 / 48 / 45 / 0. Its MAD check alone flags 16 samples, all among the 45; judging would flag 45, 29 more than that.

   The command to run on the cluster; it should print 210:

   ```
   find /proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/tier1/nfcore_rnaseq/results -path '*fastqc/raw/*_raw_fastqc.zip' | wc -l
   ```
7. **Keep tg_pct_quality_trimmed's usable line at 20 %, labelled an analogy to q30_rate's 20 % of bases. Recommended.** The research's reasoning for it is refuted; see the last section. Nothing on Step-1 comes near the line: the maximum is 1.05 %.
8. **For salmon_strand_bias, draw Salmon's 0.49 to 0.51 strip for runs quantified as U, display only, through `plot_unstranded_tolerance = 0.01`. Recommended.** Leave ctx's `<` as it is. nf-core's docs and schema say "less than" while its code uses `<=`, and no Step-1 call moves either way.
9. **For salmon_compatible_ratio, change the registry unit to %. Recommended.** The column holds 100.0. Fix the "fraction" in the metrics `_README.md` where the README is written, not in the shipped file.
10. **Fix the definitions of top_gene_share and rrna_fraction before giving either a band. Recommended.**
    - top_gene_share: leave out the Mt_rRNA genes and the 18S-antisense lncRNAs ENSG00000280614 and ENSG00000281383. The median falls from 5.34 % to 3.22 %; the maximum is 17.9 %, DEFA1B in BM_1.
    - rrna_fraction: count 18S reads, either through those two genes in unstranded libraries or against an rRNA reference. GSE100925 is reverse-stranded, so the antisense route works only for unstranded cohorts.
    - Until then, set rrna_fraction's action to keep. Its 11 MAD flags sit on shares of 1.25e-5 to 6.8e-5. Keep removes those 11 flag rows and unflags 4 samples.
11. **Leave the 53 globin MAD flags as they are. Recommended.** `checkpoint.judge` only adds band flags, never removes one. A rule letting `good` veto a MAD flag would remove 50 flag rows and unflag 28 samples, and it would change what every band in the registry means.
12. **Label derived numbers course_convention, for percent_gc's 49.8 % as well. Recommended.** The skeptic ruled that 33.3 % is derived from Trim Galore's code, not printed by it. 49.8 % is likewise a mean computed from MultiQC's curve file, so I label it course_convention, where the research had tool_code. No skeptic asked for this; I made the change so the two lines are labelled the same way.
13. **No band on replicate_spearman now; replace "UNVERIFIED" in its rationale with the confirmed ENCODE quote. Recommended.** If the column is ever recomputed ENCODE's way, show 0.8 as display only, and either drop same-donor pairs or hold them to 0.9. 39 of 210 samples share a donor with a peer, in 7 of 22 groups.
14. **No line on mt_fraction. Recommended.** The 5 % line comes from a per-cell single-cell filter. Mercer et al. 2011 give about 5 % as the bulk lower bound for white blood cells, the opposite sense.
15. **Shade only the half-planes on sex_check. Recommended.** The caption names ERX6994414 and ERX6994457, the two HSCT samples with the strongest cross-sex signal. A "both expressed" region would be unsourced.
16. **Fix `_show` in checkpoint.py as a separate task. Recommended.** It formats with `.3g`, so 21 tg_pct_passing reasons read "100 is below 100". The fix belongs in `_show`, for every value near 100.

## Refuted citations and corrected counts

In each case below, the proposal uses the skeptic's correction.

| Metric | What the research said | What the source or the recount says |
|---|---|---|
| tg_pct_quality_trimmed | Good at 5, rounded up from Liao and Shi's 2.3 to 4.6 % | That range spans two tools. Supplementary Table S1 gives 4.6 for Trimmomatic with sliding-window trimming and 2.3 for Trim Galore, on SEQC UHRR and HBRR, 100 bp paired-end, Trim Galore 0.6.2, quality and adapter trimming together. Good is now 2.3 |
| tg_pct_quality_trimmed | A run past 20 % would also fail the course's Q30 line | That does not follow. BWA trimming removes every base after the cut point, Q30 or not, and q30_rate exists only for the 8 teaching runs |
| tg_pct_quality_trimmed | Trim Galore guide: "running sum of `(cutoff - Q)`", given as v2.3.0 | That wording is on the live site and the dev branch. The v2.3.0 tag reads `(Q - cutoff)`, and the sentence begins "Trimming uses the BWA algorithm:" |
| tg_pct_with_adapter | 33.3 % as tool_code; median excess over it 5.2 points | Trim Galore prints one expectation per length and never this sum, so the line is derived and labelled course_convention. The median excess is 5.13 points |
| tg_pct_with_adapter | The 12 bp share read from the Trim Galore JSON | D-028 makes MultiQC the only parser of tool logs. MultiQC's `trim_galore_adapter_length_plot.txt` gives the same values for 210 of 210 samples |
| tg_pct_with_adapter | Shading a new 12 bp metric changes no flag | Every new registry metric gets the MAD check, which would flag 16 samples |
| tg_pct_passing | 99.9 / 99 by analogy to FastQC's overrepresented limit | Adapter dimers are at most 0.00026 %, and most too-short reads are not one sequence. The source is now Trim Galore's own example report, 552 of 1,000,000 |
| salmon_compatible_ratio | Judge at good = usable = 95 | Salmon only warns here. Judging at 95 turns that warning into a fail, so the proposal is now a display line |
| salmon_percent_mapped | Quote begins "Below 30%" | The comment opens "The thresholds are judgement calls, not statistics", and both codes are severity warning; the source text now carries both |
| fastqc_per_base_quality | Mode judge with good 28, usable 20 | 28 / 20 is FastQC's plot background, not its verdict. The verdict uses lower quartile 10 / 5 and median 25 / 20 |
| read_length | plot_bands [75, 50] from ENCODE and Salmon | Salmon's 75 bp is advice on choosing the index k, not on read quality, and 13 samples sit exactly on it. Only ENCODE's 50 nt line is kept |
| total_sequences | Judging adds 62 flags and 28 samples on the full cohort | Hypothetical: total_sequences has no full-cohort column. The 62 and 28 are what judging tg_total_reads at ENCODE's lines would do |
| fastp_surviving | Dropped at 40 % of bases below Q15; worst Step-1 loss 1.03 % | Dropped at more than 40 %. The 1.03 % is the worst of the 8 teaching runs; the full cohort has no fastp data |
| library_size | GTEx's 10 M as a single-end line; the column is "reads counted on genes" | GTEx counted mapped reads in paired-end libraries, and whether it counted mates is unverified, so that line is withdrawn. The column equals Salmon's transcriptome-mapped reads: the ratio runs from 0.99996 to 1.0001 |
| library_size | RNAseqQC rule at 0.1 times the dataset median | 0.1 is the researcher's reading of "same order of magnitude", not a number in RNAseqQC |
| mt_rrna_fraction | 15 samples exceed 0.1 | 6 do: ERX6994412, ERX6994416, ERX6994417, ERX6994429, ERX6994434 and ERX6994497 |
| mt_rrna_fraction | RNA-SeQC's 0.1 belongs to a combined rRNA rate, not here | RNA-SeQC's GTF.cpp matches transcript_type "rRNA". On GENCODE that counts Mt_rRNA, rRNA and rRNA_pseudogene, which is this column plus rrna_fraction to within 0.007 points |
| rrna_fraction | Every rRNA cut-off found measures total rRNA | Not true of RNA-SeQC on GENCODE, which does not count 18S either |
| globin_fraction | Peripheral-blood maximum 0.097 %, a monocyte | ERX6994457 is a bone-marrow monocyte. The peripheral-blood maximum is 0.062 %, an LDN sample |
| genes_detected | Median 22,435 | The median is 22,418 |
| top_gene_share | Median 5.35 % | The median is 5.34 % |
| mt_rrna_fraction | Median 1.65 % | The median is 1.64 % |
| mt_fraction | A 5 % judging line flags mostly monocytes; its authors reject it | All 5 newly flagged samples are monocytes. Osorio and Cai say the 5 % threshold should be reconsidered, not that they reject it |
| gene_strand_bias | Relabel source_type as pipeline_default | Rejected: D-039 made this transfer a course convention, and the quantity differs from both of nf-core's uses. The registry's claim that featureCounts leaves opposite-strand reads out holds only at -s 0; a 4-read toy gives 0.667 against RSeQC's 1.0 |
| salmon_strand_bias | Change ctx's `<` to nf-core's `<=` as a fix | Optional, not a fix: nf-core's own docs and schema say "less than" |
| salmon_percent_mapped | GTEx "Read mapping rate >0.2" | Relayed by PanKbase only; the GTEx primary was never opened. Rejected |
| replicate_spearman | Step-1 peers are anisogenic | 39 of 210 samples share a donor with a peer, and those pairs are isogenic by ENCODE's definition |
| fastqc_gc_distribution | Sorted cells dominated by few transcripts can warn | The researcher's inference. FastQC says only that its model assumes a random, diverse library |

## Sources, quoted

All were retrieved on 2026-10-01.

- **S1** ENCODE, Bulk RNA-seq Data Standards, updated June 2017, https://www.encodeproject.org/data-standards/rna-seq/long-rnas/: "Each replicate should have 30 million aligned reads". ENCODE's audit at https://github.com/ENCODE-DCC/encoded/blob/dev/src/encoded/audit/analysis.py, commit 4638915d from 2022-11-15, uses `upper_threshold=30000000, middle_threshold=20000000`. It gives level `WARNING` for 20 to 30 million and `NOT_COMPLIANT` below 20 million.
- **S2** The same ENCODE page: "The read length should be a minimum of 50 base pairs." The audit is at https://github.com/ENCODE-DCC/encoded/blob/dev/src/encoded/audit/experiment.py: `check_file_read_length_rna(f, 50,`, with "insufficient read length" at `NOT_COMPLIANT`.
- **S3** MultiQC 1.35, https://github.com/MultiQC/MultiQC/blob/v1.35/multiqc/modules/fastqc/fastqc_theoretical_gc/fastqc_theoretical_gc_hg38_txome.txt, headed "Human Transcriptome (UCSC hg38)". Its mean, computed this session, is 49.81.
- **S4** fastp 1.1.0 `--help`, the installed binary: -q "Default 15 means phred quality >=Q15 is qualified."; -u "Default 40 means 40%".
- **S5** FastQC 0.12.1, https://github.com/s-andrews/FastQC/blob/v0.12.1/uk/ac/babraham/FastQC/Graphs/QualityBoxPlot.java: `getY(28), baseWidth, getY(20)-getY(28)`.
- **S6** FastQC 0.12.1, https://github.com/s-andrews/FastQC/blob/v0.12.1/Configuration/limits.txt: "adapter warn 5", "adapter error 10", "n_content warn 5", "n_content error 20", "sequence warn 10", "sequence error 20", "gc_sequence warn 15", "gc_sequence error 30", "overrepresented warn 0.1", "overrepresented error 1".
- **S7** FastQC 0.12.1 Help, Per Base Sequence Content: "Nearly all RNA-Seq libraries will fail this module because of this bias". Hansen, Brenner and Dudoit 2010, Nucleic Acids Res 38:e131, doi:10.1093/nar/gkq224: "the first 13 positions".
- **S8** Trim Galore 2.3.0, https://github.com/FelixKrueger/TrimGalore/blob/v2.3.0/src/report.rs, line 638: `let expect = stats.total_reads as f64 / 4f64.powi(length as i32);`.
- **S9** Trim Galore 2.3.0 length guide, https://www.trimgalore.com/guide/length/, example report: "shorter than the length cutoff of 20 bp: 552 (0.1%)".
- **S10** Liao and Shi 2020, NAR Genom Bioinform 2(3):lqaa068, doi:10.1093/nargab/lqaa068. Supplementary Table S1, row TrimGalore, Trimmed (%): 2.3 for UHRR and 2.3 for HBRR.
- **S11** nf-core/rnaseq 3.27.0. `nextflow.config` sets `stranded_threshold` 0.8 and `unstranded_threshold` 0.1. https://github.com/nf-core/rnaseq/blob/3.27.0/docs/usage.md: "These thresholds apply to both the strandedness inferred from Salmon outputs".
- **S12** Salmon 2.7.0, https://github.com/COMBINE-lab/salmon/blob/v2.7.0/crates/salmon-core/src/diagnostics.rs: "The thresholds are judgement calls, not statistics"; "Below 30% something is usually wrong; below 10% it is almost certainly the wrong reference."
- **S13** Salmon 2.7.0, https://github.com/COMBINE-lab/salmon/blob/v2.7.0/crates/salmon-core/src/libtype.rs: `1.0 - compatible_fragment_ratio > 0.05`, "greater than 5% of the fragments disagreed with the expected library type"; `(self.strand_mapping_bias - 0.5).abs() > 0.01`.
- **S14** RNA-SeQC 2.4.2, https://github.com/getzlab/rnaseqc/blob/v2.4.2/python/rnaseqc/report.py: `'rRNA Rate': 0.1` and `'rRNA Rate': 'gt'`. src/GTF.cpp: `const boost::regex ribosomalPattern("rRNA");`.
- **S15** Shin et al. 2014, PLoS One 9(3):e91041, doi:10.1371/journal.pone.0091041: "mean 81% vs. 17% of reads mapping to globin transcripts". Whole blood; used as scale only.
- **S16** ENCODE, the S1 page: ">0.9 between isogenic replicates and >0.8 between anisogenic replicates". ENCODE's audit only warns, with "low replicate concordance" at level `WARNING`.

## Not done

- **Figure sizes.** The request relayed with this run was that static PDF and SVG figures are too small and should be publication-ready on a transparent background. This proposal does not address it, and it is still open.
- **Full-cohort FastQC and fastp.** None is on the laptop, so every FastQC count here covers the 8 teaching runs only.
- **gene_strand_bias.** It has no Step-1 values, and how STAR GeneCounts treats reads on opposite-strand genes was not tested.
- **Unopened primaries.** The GTEx primary for PanKbase's 0.2 mapping rate and the Subread Users Guide were not opened. Whether RNA-SeQC v1 counted both mates in GTEx v6p was not verified.
- **E-MTAB-11190 protocol.** It was not opened. The claim that its short inserts come from Nextera chemistry, rather than from sorted blood, is an inference.
- **The registry.** Nothing in `conf/qc_metrics.toml` or `src/` was edited.
