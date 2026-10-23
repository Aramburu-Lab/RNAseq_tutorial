# NF.1: M2's read stage and M3 on nf-core/rnaseq 3.27.0's outputs (D-113)

Draft for Rytis, 2026-10-23. Nothing here is built until he has read it and answered the four
questions at the end. Decided already: every read-level result comes from nf-core/rnaseq 3.27.0
(D-113), and M3's strand call is RSeQC `infer_experiment` (Rytis 2026-10-23).

Sources: the inventory of today's M2 and M3 (an agent's read of `src/ctx/qc`, `src/ctx/process`,
`conf/qc_metrics.toml` and the tests, 2026-10-23), and nf-core/rnaseq 3.27.0's output and
parameter pages (`nf-co.re/rnaseq/3.27.0/docs/output`, `/parameters`, read 2026-10-23). The
paths below are the output page's; the first finished run (NF.2) is checked against them before
any reader is written, because our launcher overrides one publish path (FastQC after trimming).

## What nf-core publishes, per sample, that the pages can read

| Result | Path under `results/` | Read today |
|---|---|---|
| Trim Galore report | `trimgalore/*_trimming_report.{txt,json}` (JSON by our `site.config`, D-025) | yes, M2 fates |
| FastQC raw | `fastqc/` per the page; `fates.has_raw_fastqc` expects `fastqc/raw/` | yes, if copied |
| FastQC trimmed | `trimgalore/fastqc/`, or `fastqc/trim/` under our override | no |
| STAR log | `star_salmon/log/*.Log.final.out` | yes, MultiQC `star` |
| Salmon on the alignment | `star_salmon/<sample>/quant.sf`; `star_salmon/salmon.merged.gene_counts_length_scaled.tsv` | no |
| Salmon alone | `salmon/...` as today | yes, M1, M4, M5 |
| samtools | `star_salmon/samtools_stats/<sample>.{sorted,markdup.sorted}.bam.{flagstat,idxstats,stats}` | flagstat only, via MultiQC |
| RSeQC strand | `star_salmon/rseqc/infer_experiment/*.infer_experiment.txt` | no |
| RSeQC other | `read_distribution/*.txt`, `junction_annotation/`, `inner_distance/` (paired only), `read_duplication/`, `bam_stat/` | no |
| Qualimap | `star_salmon/qualimap/<sample>/rnaseq_qc_results.txt` | no |
| dupRadar, Preseq | `star_salmon/dupradar/`, `star_salmon/preseq/*.lc_extrap.txt` | no |
| featureCounts | `star_salmon/featurecounts/*.featureCounts.txt.summary`, biotype only | no |
| Resources per task | `pipeline_info/execution_trace_*.txt` | copied, not shown |
| Genome BAM | `star_salmon/<sample>.markdup.sorted.bam`, published by default | never copied to a laptop |

Not produced, so gone from M3: STAR's per-gene counts (ReadsPerGene), featureCounts gene counts
at `-s 0/1/2` and at `-O`, `-O --fraction`, `-M`, ctx's forward-share strand call, and the
featured sample's re-alignments (limit 1, half genome, untrimmed, unspliced).

## M2's read stage, figure by figure

| Today | Input today | After NF.1 |
|---|---|---|
| Read fates bubble, read flow Sankey | Trim Galore, Salmon | unchanged; Salmon on the alignment added as a fate where present |
| FastQC curves, overrepresented sequences | FastQC raw via MultiQC | unchanged, read from nf-core's raw FastQC; the copy list gains it |
| Run band map, run strips, read tables | FastQC, fastp, seqkit per run | FastQC and Trim Galore per run; `q30_rate` (seqkit), `fastp_surviving` and `fastp_passed_reads` leave the registry's read stage; Trim Galore's `tg_*` already cover survival |
| fastp parameter grid | 36 fastp runs by ctx | Question 1 |
| Sample-level views, status, confounding | counts, trim and quant metrics, M3's table | unchanged; the align stage reads NF.1's M3 table |

## M3, figure by figure

| Today | After NF.1 | From |
|---|---|---|
| F1 mapping outcome bars | unchanged | STAR `Log.final.out` |
| Alignment table | unchanged columns; reads in and out of trimming from Trim Galore, not fastp | STAR, Trim Galore, samtools flagstat |
| F2 strand: assigned reads at `-s 0/1/2` | per sample, RSeQC's three fractions (sense, antisense, undetermined) as stacked bars, the call drawn at nf-core's thresholds | `infer_experiment.txt` |
| F2b forward share, strandedness table | the call per sample: stranded if one direction is at least `stranded_threshold` 0.8, unstranded if both are within `unstranded_threshold` 0.1 of each other, else undetermined; Salmon's `expected_format` beside it (D-039's second opinion kept) | RSeQC, `salmon/*/lib_format_counts.json`; both thresholds are nf-core 3.27.0's defaults (`nextflow.config`), registered in `conf/qc_metrics.toml` with that source |
| F3, F4: five counter comparisons (STAR vs featureCounts, s0 vs called, `-O`, `--fraction`, `-M`) | one comparison: Salmon on STAR's alignment against Salmon alone, per gene, Spearman with CI and Bland and Altman bias and limits (`concordance.per_unit`, kept) | `star_salmon/` and `salmon/` matrices; this keeps learning outcome 3, "what an aligner, a pseudo-aligner and a read counter each produce, and why counts differ" |
| F1b featured-sample settings, incomplete-reference table | Question 2 | |
| `--compare-dataset` context | unchanged | gather counts, qc metrics |
| Reference table (`reference.json`, start codons) | dropped: ctx no longer prepares the genome; the page names the GENCODE release and nf-core's FASTA and GTF from the run's params | `pipeline_info/params_*.json` |
| Folders table | dropped | |
| Cluster steps table | from nf-core's execution trace: per process, tasks, wall time, peak RSS | `pipeline_info/execution_trace_*.txt` |
| New figures from nf-core QC | Question 3 | |

## Code

1. `ctx gather`: `COPIED_FILES` gains the files above that the pages read, never a BAM; `AGREE`
   gains `aligner`; `skip_alignment` false is accepted beside true. One copy command for the Mac
   (an rsync include list) in the plan's Mac block, its byte count measured on the first run.
2. One reader for an nf-core outdir's alignment, `ctx.process.nfcore` (STAR logs, samtools,
   RSeQC, the two Salmon matrices, the trace), through MultiQC where MultiQC has the module
   (`star`, `samtools`, `rseqc`, `salmon`, `featurecounts`) and directly where it does not.
3. `ctx process --rnaseq-dir`, as `ctx qc` already takes it, replaces `--align-dir`; the planner
   (`pipeline.plan`) wires it the same way.
4. Removed once NF.1's pages pass: `process/cluster.py`, `cluster_route.py`, `fc_compact.py`,
   `qc/cluster.py`, the featureCounts settings matrices in `process/count.py` and their tests,
   `--align-dir` on both commands and its cards. The ED.6 cluster plan is marked superseded.
   Question 2 decides whether ctx's own HISAT2 route on chr22 stays as a laptop exercise.
5. Registry (`conf/qc_metrics.toml`): `rseqc_sense`, `rseqc_antisense`, `rseqc_undetermined` and
   the strand call, sourced to nf-core 3.27.0; seqkit and fastp read-stage metrics removed.

## Fixture and tests

- A committed nf-core 3.27.0 outdir with STAR: nf-core's own `-profile test,docker` (yeast) on
  this Mac with `aligner: star_salmon`, cut to the files the pages read by
  `scripts/python/make_rnaseq_fixture.py`, as P2A.2 made today's Salmon-only fixture. Its
  README states the run, the versions and the SHA-256 of every file.
- Per reader: a test on the fixture and a planted mutation (a fraction column swapped, a
  threshold read from anywhere but the registry, the two Salmon matrices swapped).
- End to end: `ctx qc --rnaseq-dir` and `ctx process --rnaseq-dir` on the fixture under full and
  core; the lite pages in Chromium at 1,280 and 375 px; page sizes under 15 MB.
- On real data, once NF.2 ends: tier1's M2 and M3 rebuilt from its nf-core run, every figure
  counted against the run's samples (210), the strand call compared with nf-core's own MultiQC
  strand check table, and the ceiling stated: samples called against samples in the samplesheet.

## Order

1. NF.2 ends; the first run's files checked against the table above.
2. The fixture.
3. gather's copy list and checks.
4. M3's reader, figures and tables; then M2's read stage.
5. The removals in Code 4, the registry, the cards, CONVENTIONS and a D-entry.
6. The six projects' pages rebuilt; the core review set rebuilt.

## Questions for Rytis

1. **The fastp parameter grid (M2).** nf-core trims with Trim Galore, and the grid ran fastp 36
   times on one run. Drop it, or keep it as a laptop-only exercise on the teaching subset?
2. **The featured sample's settings and the incomplete reference (M3).** Limit 1, half genome,
   untrimmed and unspliced need ctx's own aligner, which a laptop runs only as HISAT2 on chr22.
   Drop them, or keep ctx's HISAT2 route as an optional laptop exercise (a second route to
   maintain)?
3. **New figures from nf-core's QC (M3).** Candidates, each per sample: RSeQC read distribution
   (CDS, UTR, intron, intergenic), RSeQC junction annotation (known, partly novel, novel
   splices), Qualimap's 5' to 3' bias, dupRadar's duplication against expression, Preseq's
   library complexity, featureCounts' biotypes. Which to add; the lite edition shows the same.
4. **What M3 stops teaching.** featureCounts' overlap and multimapper settings (`-O`,
   `--fraction`, `-M`) have no nf-core counterpart. Accept that, with multimappers taught from
   STAR's own multimapped share in F1?
