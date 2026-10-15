#!/bin/bash
# Package P: one set of builds for a label from one src: tier1 de, function, report; breast de.
#   bash run_set.sh <label> <src dir> [gather]
# With "gather", the base roots' gathers run first (once).
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a4bc57162a0d35df8
B=$WT/docs/plans/logs/wo13_p/build.sh
label=$1 src=$2
if [ "$3" = gather ]; then
  bash "$B" base "$src" tier1 gather &
  bash "$B" base "$src" breast gather &
  wait
fi
bash "$B" "$label" "$src" tier1 de function report &
bash "$B" "$label" "$src" breast de &
wait
echo "set $label done"
