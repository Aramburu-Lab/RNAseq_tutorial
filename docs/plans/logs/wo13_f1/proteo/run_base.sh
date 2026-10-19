#!/bin/bash
# wo13-f1 M7: ctx proteo with BASE's src (25baea5, git archive into the scratchpad) on one
# scratch root of this worktree, timed: the run before this change, for a byte comparison.
# Usage: run_base.sh <base src parent> <scratch name> <dataset> <dataset file> [options]
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a0cad45c60690dabd
base=$1 name=$2 dataset=$3 spec=$4
shift 4
export PATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH
export PYTHONPATH=$base/src
export PYTHONDONTWRITEBYTECODE=1
cd "$base" || exit 1
uptime
/usr/bin/time -l python -m ctx --root "$WT/outputs/scratch_roots/$name" proteo \
  --dataset "$dataset" --dataset-file "$spec" "$@"
echo "exit $?"
