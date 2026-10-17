#!/bin/bash
# Package wo13-de: lung de and qc, main's code and this worktree's in turn, three times each,
# for the spread of repeated runs (lessons.md: state the spread, not one run).
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ab82f97e3ba379d33
here=$WT/docs/plans/logs/wo13_de
before=$WT/outputs/scratch_roots/src_main/src
for i in 1 2 3; do
  for side in before after; do
    label=rep${i}_$side
    bash "$here/roots.sh" "$label" > /dev/null
    src=$WT/src
    [ "$side" = before ] && src=$before
    SRC=$src bash "$here/step.sh" "lung_$label" tcga_lung de > /dev/null
    SRC=$src bash "$here/step.sh" "lung_$label" tcga_lung qc --strata project --peers project > /dev/null
  done
done
for f in "$here"/lung_rep*_{de,qc}.log; do
  echo "$(basename "$f" .log) $(grep -E 'real|maximum resident|peak memory|^exit' "$f" | awk '{print $1}' | tr '\n' ' ')"
done | tee "$here/repeat_summary.log"
