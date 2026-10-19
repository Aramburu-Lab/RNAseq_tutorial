#!/bin/bash
# wo13-f1: M5's and M6's new kinds redrawn from Step-1's pages in clean folders (redraw.py).
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a00e5034a0fa89898
PY=$WT/docs/plans/logs/wo13_f1_py.sh
out=$WT/outputs/scratch_roots/tier1_f1/outputs
clean=$WT/outputs/scratch_roots/redraw_step1
mkdir -p "$clean"
for pair in "05_de blood_m5_dispersion" "05_de blood_m5_pca" "05_de blood_m5_box" \
            "06_function blood_m6_running" "06_function blood_m6_background"; do
  set -- $pair
  page=$(ls "$out/$1"/*_tier1_overview_v1.html)
  echo "== $2"
  "$PY" python "$WT/docs/plans/logs/wo13_f1/redraw.py" "$page" "$2" "$clean/$2" 2>&1 | grep -v "^$"
done
