# M4 Quantification on real data: plan (2026-10-07)

**Outcome.** The featured M4 page can show all 210 Step-1 runs from Salmon outputs nf-core/rnaseq already wrote on Kebnekaise: one copy, no new quantification. Library type, decoys, mapping, sample correlation and length scaling's invariant work today; transcript lessons need copy Q1; the cross-method comparison needs STAR on the same reads (M3.11 or M3.13). Numbers are measured today unless marked UNVERIFIED.

## 1. Data

**Featured dataset: Step-1 (`tier1`).** D-050 puts the full cohort on read-level pages ([decisions.md:98](knowledge/decisions.md)); M3 and M5 use Step-1 too. The 8-run subset stays the laptop exercise, out of the kit.

**On this Mac**

| Path | Holds | Size |
|---|---|---|
| `~/ctx_tier1/data/raw/tier1_rnaseq` | Length-scaled gene matrix (78,470 × 210); `meta_info.json` and `lib_format_counts.json` for each sample; Trim Galore JSON reports; `pipeline_info` | 243 MB |
| The breast outdirs | Same layout (22, 94 and 103 samples), no dataset file yet | 14, 49 and 65 MB |
| `resources/gencode/` | GENCODE 50 GTF | 124,650,284 B |
| M3's `trim_b72d1253088d/reads` | The 8 trimmed runs | 465,927,339 B |

None of them holds a `quant.sf` file, a transcript table, a transcript-to-gene table or unscaled gene counts; `COPIED_FILES` leaves those on the cluster ([rnaseq.py:29-35](src/ctx/gather/rnaseq.py)).

**What these files already show**

| Fact | Value |
|---|---|
| Index | One for all 4 datasets: Salmon 2.7.0, 644,914 targets, 194 decoys |
| Bias correction | Off in 429/429 runs; optimiser `vb` |
| Library type | Step-1: U in 210/210 runs, strand bias 0.494–0.514. breast_tissue: SR in 103/103. Paired-end breast sets: IU |
| Fragment length | Single-end runs use the prior 250 ± 25. Paired-end runs are observed: 111–239 nt in breast_mono |
| Step-1 mapping | 29.3–90.0 %, median 78.6 %. Decoys take a median 17.7 % of processed reads (5.5–52.0 %) |
| Length-scaled column sums | Equal to `num_mapped` within 3.3e-8 relative, in 429/429 runs |
| Reads processed − mapped − decoy − filtered | At least 11,770, in 429/429 runs |
| GENCODE 50 genes without a row | 471 of 78,941. 112 are on chrY, 62 of them in the PARs; why is UNVERIFIED |
| Genes with ≥ 10 length-scaled counts in any sample | 42,539. Protein-coding: 17,879 of 20,131 annotated |
| Cluster trace | Salmon 64–1,249 s at 10.2 GB per sample; R tximport, all 210: 10.8 GB |

**Cluster block Q0** (Rytis; read only, about 40 lines):

```
# Q0. Kebnekaise login node, read only. Paste back all of it.
P=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/tier1/nfcore_rnaseq/results
if test -s "$P/salmon/salmon.merged.gene_counts_length_scaled.tsv"; then
    find "$P/salmon" -maxdepth 1 -type f -printf '%s %f\n' | sort -k2
    for f in quant.sf quant.genes.sf aux_info/ambig_info.tsv libParams/flenDist.txt cmd_info.json; do
        printf '%s: ' "$f"; find "$P/salmon" -mindepth 2 -path "*/$f" -printf '%s\n' | awk '{s+=$1} END {print NR, "files,", s, "bytes"}'
    done
    cat "$(find "$P/salmon" -mindepth 2 -maxdepth 2 -name cmd_info.json | sort | head -1)"
    find "$P/genome" -maxdepth 3 -type f -printf '%s %P\n' | sort -k2
else
    echo "NO MATRIX under $P/salmon: nothing checked"
fi
```

Q1 and Q2 are written in full from Q0's output, since both depend on what it prints ([lessons.md:47](knowledge/lessons.md)).

- **Q1** is step 7's rsync with 11 more includes: four gene tables, `transcript_*.tsv` and `tx2gene.tsv`; for each sample, `quant.sf`, `ambig_info.tsv`, `flenDist.txt` and `cmd_info.json`; and `duplicate_clusters.tsv`.
  - Check: the matrix's SHA-256 must stay the hash gather recorded, 630e5514…eaa5fd.
  - Size: UNVERIFIED. R's tximport task read 5.9 GB and wrote 1.8 GB, which bounds `quant.sf` and the tables.
- **Q2** comes after M3.11 uploads the 8 trimmed runs. It runs nf-core's own Salmon image and index on them twice: once with nf-core's options plus `--dumpEq`, once with `--gcBias --seqBias`.
  - Request: 4 cores, 16 GB, 1 h. The basis is nf-core's own subsampled strandedness runs, 15–37 s each at 10.0–10.3 GB.

**Downloads that need your permission**

| What | Source | Bytes | Licence |
|---|---|---|---|
| Q1 and Q2 | Your cluster | Q0 measures | E-MTAB-11190: the IDF names no licence ([TODO.md:2442](TODO.md)) |
| Yeast test data, to rebuild `tests/fixtures/rnaseq` | nf-core/test-datasets@626c8fab, per 3.27.0's `conf/test.config` | 24,675,021 for its 14 input files (HEAD) | MIT |
| nf-core/rnaseq 3.27.0 code and any uncached image | github.com/nf-core/rnaseq | UNVERIFIED | MIT |
| pytximport 0.13.0 (option 3b only) | bioconda | 51,475, plus 5 dependencies | GPL-3.0-only |

## 2. Analyses and figures

Figures reuse the preview's contracts ([preview_quant.py:168-386](scripts/python/preview_quant.py)). M2 already judges Salmon's metrics ([qc_metrics.toml:418](conf/qc_metrics.toml)) and read fates ([fates.py:77-90](src/ctx/qc/fates.py)); M4 explains them, never judging twice.

| # | Figure, contract | Data, when | Statistics | Control (D-051) |
|---|---|---|---|---|
| F1 | Salmon's fragment fates: mapped, decoy, below score, no hit. `BarData` | MultiQC table, **now** | Mann-Whitney by batch within cell type, rank-biserial r, BH. Unadjusted today: NDN p 0.069, r 0.25; monocytes 0.026, 0.31 | subset |
| F2 | Genes detected against depth. `CurveData` | Binomial thinning of `gene_counts`, seed 42, **Q1** | Mann-Whitney at the common depth | depth; reads per gene |
| F3 | Spearman between samples. `HeatmapData` | Length-scaled counts, **now** | Within − between cell types, labels permuted within donor ([annotate.py:241](src/ctx/viz/annotate.py)). The preview permutes freely ([preview_quant.py:227-241](scripts/python/preview_quant.py)) | group by |
| F4 | Reads against TPM by effective length. `ScatterData` | Gene TPM, counts and lengths, **Q1** | descriptive | count table |
| F5 | Genes detected by Salmon, STAR and featureCounts. `BarData` | **M3.11 + Q2** (8 runs); M3.13 (210) | counts; `concordance.agreement` | threshold; gene type |
| F6 | `ndn_gcsf` fold change, estimates against length-scaled counts | **Q1** | `concordance.agreement` | gene type |
| F7 | Isoforms expressed against isoforms annotated | Transcript TPM, GTF, **Q1** | Mann-Whitney, coding against lncRNA | gene type |
| F8 | Reads Salmon splits between genes | Q2's `--dumpEq`, 8 runs | descriptive | gene type |

**Laptop exercise** (`tier1_teach`, chr22; 4 threads, load average 4–8):

| Step | Measured |
|---|---|
| Index | 1.72 s at 392,904,704 B; with chr22 decoys, 8.48 s at 1,693,171,712 B |
| Quantification | 1.07–3.33 s per run, ≤ 239,599,616 B; all U; 2.55–3.20 % mapped (1.92–2.69 % with decoys) |
| Bias correction, 2 runs | 2.66–3.94 s; gene totals within ±0.31 log2 (95 %), transcripts up to 405 reads |

**Sanity checks**

- **Ceilings.** Genes counted ≤ 78,470 rows ≤ 78,941 GENCODE 50 genes; isoforms expressed ≤ annotated (protein-coding median 15).
- **Invariants, before any page.** Length-scaled column sums equal `num_mapped`; every read has a fate; ctx's gene tables equal nf-core's R tables; single-end Length − EffectiveLength is 248.7–250 on transcripts ≥ 300 nt (exactly 250 for 14,453 of 14,646 on chr22).
- **Planted truth (`chr22_se`).** SR ×4 and U ×1, as planted. In se_a, se_b and se_d every gene without a same-strand neighbour equals its planted reads, and the overlapping trio its planted total (458, 425, 436).
- **Published context, not a ceiling.** Against qPCR on 13,045 genes, Salmon reached R 0.845 and STAR-HTSeq 0.821; on fold changes, 0.929 and 0.933 ([Everaert 2017](https://doi.org/10.1038/s41598-017-01617-3); stranded, paired-end reads).
- **What chr22 cannot teach.** Genes that only chr22 Salmon detects get 4.1–19.3× the reads the full index gives them. A laptop Venn would teach an incomplete reference, not multimapping.

## 3. The build

**`src/ctx/quant/`.** Each file stays under the 600-line ratchet.

- `tximport.py`: the one gene summary, in tximport's order ([summarizeToGene.R](https://github.com/thelovelab/tximport/blob/devel/R/summarizeToGene.R)).
- `salmon.py`: transcriptome with M3's ATG check; keyed index and quantification folders under `outputs/04_quant/tools/`; `salmon_probe`, as D-036 does for STAR.
- `outdir.py`: Polars readers by column name; logs only through `multiqc_io`.
- `checks.py`, `tables.py`, `figures.py`, `page.py`, `run.py`, as in M3.

**CLI:** `ctx quant --dataset --dataset-file [--rnaseq-dir | --reads] [--genome] [--counts-from-abundance] [--libtype] [--gc-bias] [--seq-bias] [--decoys] [--compare-dataset] [--registry]`.

- The laptop path finds M3's trimmed reads through M3's newest `folders` table, or refuses.
- M4 writes only under `outputs/04_quant/`.

**Tests, each with the mutations it catches**

- `test_tximport_matches_r`: ctx against R tximport 1.30.0, run in the cached nf-core image on the chr22_se runs and a planted 4-sample set. My 30-line port matched R within 3.5e-15 relative. Mutations: arithmetic for geometric mean; row means before the fill; scaling to TPM; count-weighted lengths; Length for EffectiveLength. The planted set catches all 5; chr22_se alone misses the first 2, so the order test needs it.
- `test_salmon_probe`: a Salmon that maps nothing is refused.
- `test_planted_truth`: the library types and gene totals above.
- `test_fates`, `test_sums`: a dropped decoy term; scaling to TPM.
- `test_restricted_permutation`: a planted donor effect that free permutation calls significant.
- `test_door`: every input problem listed before Salmon runs.
- End to end on both fixtures; `-m budget` on Step-1.

**Budgets.**

- `runtime.py:26` makes quant Track A; the laptop path's measured steps fit Track B too.
- The cohort path must stream: R needed 10.8 GB. Polars summed a synthetic 644,914-row `quant.sf` in 0.032 s, 661 MB for 10 files.
- The page's size against the 15 MB cap is unmeasured; a 210 × 210 heatmap is the risk.

**From the preview into `src/`.** The figure builders and `rho_levels` move to `figures.py`, the page to `page.py`; `preview_quant.py` keeps its generators and calls `ctx.quant.page`, as `preview_process.py` does. `quant.toml` is rewritten from measurements, its planned cards move to `cards.toml`, and `glossary.yaml` gains decoys.

## 4. Dependencies and order

1. Now, in parallel: Q0 and M4.2 to M4.4.
2. Then decision 2 and Q1. `COPIED_FILES` and every `_rnaseq/` rsync line change together ([test_repo.py:262-283](tests/test_repo.py)), and the yeast fixture is rebuilt ([test_rnaseq_fixture.py:52](tests/test_rnaseq_fixture.py)).
3. After Q1: F2, F4, F6, F7. After M3.11 and Q2: F5. Last: controls.
4. `cli.py`, `cards.toml` and `glossary.yaml` are shared with other planners; merge in turn. M8 can reuse Q1's transcript tables.

## 5. Decisions for Rytis

1. **Featured dataset.** Recommended: **Step-1**. Cost: single-end effective lengths rest on a prior. The breast datasets have no dataset file yet.
2. **What Q1 copies.** Recommended: **all of it**; Q0 measures the size. Merged tables alone lose per-sample ambiguity; nothing leaves only F1 and F3.
3. **Gene summaries.** Recommended: **ctx's own port, checked against R**: about 40 lines and a D-entry, since CONVENTIONS.md:33 names pytximport. pytximport brings GPL-3.0, 5 packages and pandas (barred by CONVENTIONS.md:29); its memory is unmeasured.
4. **Cross-method comparison.** Recommended: **M3.11 adds featureCounts; Q2 quantifies the same reads**. The cohort's Salmon instead mixes in trimmer differences.
5. **Budget.** Recommended: **Track A, plus a Track B test on the cohort path**.
6. **Salmon defaults.** Recommended: **nf-core's, with bias correction off and offered as a knob**. The preview's cards say it is on ([quant.toml:67-79](src/ctx/report/content/quant.toml)).
7. **F4 has no test.** Recommended: **extend M3's R6 exception to it**.

## 6. Steps

- [ ] M4.0 This plan, a D-entry, and the corrections below (S)
- [ ] M4.1 Q0 (S, Rytis)
- [ ] M4.2 The fixture builder: Salmon, then R tximport (M)
- [ ] M4.3 `tximport.py` and its tests (M)
- [ ] M4.4 `salmon.py` and its probe (M)
- [ ] M4.5 Q1; `COPIED_FILES`; the rsync lines; the yeast fixture (M)
- [ ] M4.6 Readers, checks and the input refusals (M)
- [ ] M4.7 Tables (L)
- [ ] M4.8 Figures, content, cards and glossary (L)
- [ ] M4.9 `ctx quant`, page, end-to-end tests (M)
- [ ] M4.10 Step-1 page; budgets; browser at 1280 and 375 px, both themes (M)
- [ ] M4.11 The laptop exercise (M)
- [ ] M4.12 Q2 (S)
- [ ] M4.13 F5 (M)
- [ ] M4.14 Controls (M)
- [ ] M4.15 The preview through `page.py`, review, lessons (S)

**Corrections found**

- [TODO.md:2530](TODO.md) calls Step-1 stranded; all 210 runs are U.
- [quant.toml:9](src/ctx/report/content/quant.toml) claims bias correction and EM; the runs used neither (VB).
- [quant.toml:16](src/ctx/report/content/quant.toml)'s paralog claim is untestable on chr22.
- [TODO.md:42](TODO.md) puts full-transcriptome Salmon in Track A; the cohort's full index needed 10.2 GB.

**Not done.** No pixi, cluster command, page build or repository edit. Unmeasured: cluster file sizes, pytximport's memory, page size. Fates are assumed disjoint from non-negative remainders, not documentation. Others were editing TODO.md and decisions.md; line numbers are as of 15:45.
