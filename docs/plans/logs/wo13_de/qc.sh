#!/bin/bash
# Package wo13-de: ctx qc on TCGA lung, as the real root ran it (--strata project --peers
# project), with main's code and with this worktree's, one after the other; then compared.
#   bash qc.sh <label>
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ab82f97e3ba379d33
here=$WT/docs/plans/logs/wo13_de
label=$1
bash "$here/roots.sh" "qc$label" > /dev/null
bash "$here/roots.sh" "qc${label}_before" > /dev/null
bash "$here/step.sh" "lung_qc$label" tcga_lung qc --strata project --peers project
SRC=$WT/outputs/scratch_roots/src_main/src bash "$here/step.sh" "lung_qc${label}_before" \
  tcga_lung qc --strata project --peers project
bash "$here/env.sh" python docs/plans/logs/wo13_de/compare_tables.py \
  "outputs/scratch_roots/lung_qc${label}_before/outputs/02_qc" \
  "outputs/scratch_roots/lung_qc$label/outputs/02_qc" | tee "$here/qc${label}_compare.log"
