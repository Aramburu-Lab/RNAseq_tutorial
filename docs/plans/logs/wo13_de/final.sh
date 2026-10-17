#!/bin/bash
# Package wo13-de: the final code measured. Step-1's 92 metrics and every table against main's
# run (tier1_main); then lung de and qc three times and breast de once, each table compared
# with main's run on the same inputs (lung_rep1_before, breast_base), one run at a time.
#   bash final.sh <label>
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ab82f97e3ba379d33
here=$WT/docs/plans/logs/wo13_de
label=$1
roots=$WT/outputs/scratch_roots
compare() { bash "$here/env.sh" python docs/plans/logs/wo13_de/compare_tables.py "$1" "$2"; }
bash "$here/step1.sh" "$label"
{
  echo "== Step-1 de: main's run against $label"
  compare "$roots/tier1_main/outputs/05_de" "$roots/tier1_$label/outputs/05_de"
  echo "== Step-1 function: main's run against $label"
  compare "$roots/tier1_main/outputs/06_function" "$roots/tier1_$label/outputs/06_function"
} > "$here/${label}_compare.log"
for i in 1 2 3; do
  bash "$here/roots.sh" "${label}$i" > /dev/null
  bash "$here/step.sh" "lung_${label}$i" tcga_lung de > /dev/null
  bash "$here/step.sh" "lung_${label}$i" tcga_lung qc --strata project --peers project > /dev/null
done
bash "$here/step.sh" "breast_${label}1" tcga_brca de > /dev/null
{
  echo "== lung de: main's run against ${label}1"
  compare "$roots/lung_rep1_before/outputs/05_de" "$roots/lung_${label}1/outputs/05_de"
  echo "== lung qc: main's run against ${label}1"
  compare "$roots/lung_rep1_before/outputs/02_qc" "$roots/lung_${label}1/outputs/02_qc"
  echo "== breast de: main's run against ${label}1"
  compare "$roots/breast_base/outputs/05_de" "$roots/breast_${label}1/outputs/05_de"
} >> "$here/${label}_compare.log"
for f in "$here"/lung_"${label}"?_{de,qc}.log "$here/breast_${label}1_de.log"; do
  echo "$(basename "$f" .log) $(grep -E 'real|maximum resident|peak memory|^exit' "$f" | awk '{print $1 $2}' | tr '\n' ' ')"
done | tee "$here/${label}_summary.log"
cat "$here/${label}_compare.log"
