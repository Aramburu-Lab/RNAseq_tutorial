#!/bin/bash
# Package wo13-de: the S.5 measurements, one at a time so no run presses on another's memory.
# "after" runs this worktree's src, "before" main's at 628ad3a (outputs/scratch_roots/src_main).
#   bash measure.sh <label>
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ab82f97e3ba379d33
here=$WT/docs/plans/logs/wo13_de
label=$1
before=$WT/outputs/scratch_roots/src_main/src
bash "$here/roots.sh" "$label" > /dev/null
bash "$here/roots.sh" "${label}_before" > /dev/null
bash "$here/step.sh" "breast_$label" tcga_brca de
bash "$here/step.sh" "lung_$label" tcga_lung de
bash "$here/step.sh" "lung_$label" tcga_lung qc
SRC=$before bash "$here/step.sh" "lung_${label}_before" tcga_lung qc
for which in breast lung; do
  echo "== $which: tables of main's code (base) against this run's" >> "$here/${label}_compare.log"
  bash "$here/env.sh" python docs/plans/logs/wo13_de/compare_tables.py \
    "outputs/scratch_roots/${which}_base/outputs/05_de" \
    "outputs/scratch_roots/${which}_$label/outputs/05_de" >> "$here/${label}_compare.log"
done
echo "== lung qc: main's code against this run's" >> "$here/${label}_compare.log"
bash "$here/env.sh" python docs/plans/logs/wo13_de/compare_tables.py \
  "outputs/scratch_roots/lung_${label}_before/outputs/02_qc" \
  "outputs/scratch_roots/lung_$label/outputs/02_qc" >> "$here/${label}_compare.log"
cat "$here/${label}_compare.log"
