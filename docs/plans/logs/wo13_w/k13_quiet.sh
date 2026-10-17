#!/bin/bash
# Package W, K.13: k13_stress.py with no load, for the baseline of each traced step.
#   bash k13_quiet.sh <label> <drivers> <deadline s> <rounds>
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aea85c67a862f2e82
SCRATCH=/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad
label=$1
log=$WT/docs/plans/logs/wo13_w/k13_stress_${label}.log
echo "== $label: drivers $2, deadline $3 s, rounds $4, no load, $(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$log"
bash "$WT/docs/plans/logs/wo13_w/py.sh" "$WT/docs/plans/logs/wo13_w/k13_stress.py" \
  "$SCRATCH/w_k13_$label" "$2" "$3" "$4" >> "$log" 2>&1
echo "== exit $?" >> "$log"
tail -20 "$log"
