#!/bin/bash
# Package wo13-de: one timed ctx step on a scratch root, its log in this folder.
#   [SRC=<another src>] bash step.sh <root label e.g. lung_base> <dataset> <step> [extra ctx args]
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ab82f97e3ba379d33
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
export PATH=$BIN:$PATH PYTHONPATH=${SRC:-$WT/src} PYTHONDONTWRITEBYTECODE=1
label=$1
dataset=$2
step=$3
shift 3
root=$WT/outputs/scratch_roots/$label
log=$WT/docs/plans/logs/wo13_de/${label}_${step}.log
echo "== $step $(date -u +%Y-%m-%dT%H:%M:%SZ) src=$PYTHONPATH" >> "$log"
/usr/bin/time -l python -m ctx --root "$root" "$step" --dataset "$dataset" \
  --dataset-file "$WT/conf/datasets/$dataset.toml" "$@" >> "$log" 2>&1
echo "exit $?" >> "$log"
grep -E "^== |^exit |real|maximum resident|peak memory" "$log" | tail -6
