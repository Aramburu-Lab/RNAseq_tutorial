#!/bin/bash
# F.1 M8: a scratch root for tcga_brca. outputs/ is a real folder holding copies of the newest
# gather and M5 tables M8 reads, with their READMEs; data and resources link the real root's.
set -euo pipefail
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a84b7900e11d976d1
real=/Users/rytis/ctx_breast
root=$WT/outputs/scratch_roots/breast_m8
mkdir -p "$root/outputs/01_gather" "$root/outputs/05_de"
ln -sfn "$real/data" "$root/data"
ln -sfn "$real/resources" "$root/resources"
for what in samplesheet counts genes mirna_samplesheet mirna_mature_counts mirna_hairpin_counts \
            mirna_hairpins; do
  for f in "$real"/outputs/01_gather/2026-10-16_gather_tcga_brca_${what}_v1*; do
    cp "$f" "$root/outputs/01_gather/"
  done
done
for f in "$real"/outputs/05_de/2026-10-16_de_tcga_brca_tumour_normal_results_v1*; do
  cp "$f" "$root/outputs/05_de/"
done
ls -l "$root/outputs/01_gather" "$root/outputs/05_de" | awk '{print $5, $9}'
