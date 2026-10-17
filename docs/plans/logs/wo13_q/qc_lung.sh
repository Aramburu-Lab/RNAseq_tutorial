#!/bin/bash
# Package Q: ctx qc on tcga_lung in the scratch root, timed by /usr/bin/time -l.
# Usage: qc_lung.sh <label>; the log goes to docs/plans/logs/wo13_q/qc_lung_<label>.log.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a7158d7cd646891a9
ROOT=$WT/.scratch_q/root
LOG=$WT/docs/plans/logs/wo13_q/qc_lung_$1.log
"$WT/docs/plans/logs/wo13_q/py.sh" /usr/bin/time -l python -m ctx --root "$ROOT" qc \
  --dataset tcga_lung --dataset-file "$WT/conf/datasets/tcga_lung.toml" \
  --strata project --peers project > "$LOG" 2>&1
echo "exit $?" >> "$LOG"
tail -25 "$LOG"
