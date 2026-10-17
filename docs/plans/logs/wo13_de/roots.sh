#!/bin/bash
# Package wo13-de: scratch roots inside this worktree. TCGA roots link the real roots' newest
# gather tables file by file; Step-1's links data/ and resources/ and gathers afresh. Every
# write lands under the scratch root's own outputs/.
#   bash roots.sh <label>
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ab82f97e3ba379d33
label=$1
for pair in lung:tcga_lung:2026-10-16 breast:tcga_brca:2026-10-16; do
  IFS=: read -r which dataset day <<< "$pair"
  real=/Users/rytis/ctx_$which
  root=$WT/outputs/scratch_roots/${which}_$label
  mkdir -p "$root/outputs/01_gather"
  [ -e "$root/resources" ] || ln -s "$real/resources" "$root/resources"
  for what in samplesheet counts genes; do
    f=${day}_gather_${dataset}_${what}_v1.parquet
    [ -e "$root/outputs/01_gather/$f" ] || ln -s "$real/outputs/01_gather/$f" "$root/outputs/01_gather/$f"
  done
done
root=$WT/outputs/scratch_roots/tier1_$label
mkdir -p "$root/outputs"
[ -e "$root/data" ] || ln -s /Users/rytis/ctx_tier1/data "$root/data"
[ -e "$root/resources" ] || ln -s /Users/rytis/ctx_tier1/resources "$root/resources"
ls -l "$WT"/outputs/scratch_roots/*_"$label"/outputs/01_gather
