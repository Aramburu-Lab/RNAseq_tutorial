#!/bin/bash
# Package wo13-de: Step-1 gather, de and function on a fresh scratch root, then its 92
# recorded metrics compared value by value (D-046).
#   [SRC=<another src>] bash step1.sh <label>
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ab82f97e3ba379d33
here=$WT/docs/plans/logs/wo13_de
label=$1
raw=/Users/rytis/ctx_tier1/data/raw
bash "$here/roots.sh" "$label" > /dev/null
bash "$here/step.sh" "tier1_$label" tier1 gather --fetchngs-dir "$raw/tier1_meta" --rnaseq-dir "$raw/tier1_rnaseq"
bash "$here/step.sh" "tier1_$label" tier1 de
bash "$here/step.sh" "tier1_$label" tier1 function
out=$WT/outputs/scratch_roots/tier1_$label/outputs
bash "$here/env.sh" python docs/plans/logs/wo13_p/compare_metrics.py \
  deliverables/2026-09-30_de_tier1_metrics_v1.csv \
  "$(ls "$out"/05_de/*_tier1_metrics_v1.csv | tail -1)" \
  "$(ls "$out"/06_function/*_tier1_metrics_v1.csv | tail -1)" | tee "$here/tier1_${label}_metrics.log"
