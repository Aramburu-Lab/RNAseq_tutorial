#!/bin/bash
# Package W, K.13: the heatmap-cycle test N times under full CPU load, each with its own
# basetemp, so no other session's pytest garbage is deleted inside the timing.
#   bash k13_repeat.sh <label> <runs>
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aea85c67a862f2e82
SCRATCH=/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad
label=$1 runs=$2
log=$WT/docs/plans/logs/wo13_w/k13_${label}.log
for i in $(seq "$runs"); do
  echo "== run $i of $runs ($label)" >> "$log"
  bash "$WT/docs/plans/logs/wo13_w/loaded.sh" "$log" -q -p no:cacheprovider \
    --basetemp="$SCRATCH/w_k13_bt_$label" --durations=2 tests/test_heatmap_cycle.py
done
grep -E "^== |passed|failed|call " "$log"
