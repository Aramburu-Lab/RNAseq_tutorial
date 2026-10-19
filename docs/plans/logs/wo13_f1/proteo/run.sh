#!/bin/bash
# wo13-f1 M7: ctx proteo with this worktree's src on one scratch root, timed.
# Usage: run.sh <scratch name> <dataset> <dataset file> [ctx proteo options]
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a0cad45c60690dabd
name=$1 dataset=$2 spec=$3
shift 3
uptime
/usr/bin/time -l "$WT/docs/plans/logs/wo13_f1_proteo_py.sh" python -m ctx \
  --root "$WT/outputs/scratch_roots/$name" proteo --dataset "$dataset" --dataset-file "$spec" "$@"
echo "exit $?"
