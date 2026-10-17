#!/bin/bash
# Package wo13-e2a: Step-1 gather and qc on a fresh scratch root inside this worktree, data/
# and resources/ linked read-only from ~/ctx_tier1; every write lands in the scratch root.
#   [SRC=<another src>] bash step1.sh <label> [steps, default "gather qc"]
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a68d1b8f18deeebdd
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
export PATH=$BIN:/usr/bin:/bin PYTHONPATH=${SRC:-$WT/src} PYTHONDONTWRITEBYTECODE=1
label=$1
steps=${2:-gather qc}
raw=/Users/rytis/ctx_tier1/data/raw
root=$WT/outputs/scratch_roots/tier1_$label
mkdir -p "$root/outputs"
[ -e "$root/data" ] || ln -s /Users/rytis/ctx_tier1/data "$root/data"
[ -e "$root/resources" ] || ln -s /Users/rytis/ctx_tier1/resources "$root/resources"
log=$WT/docs/plans/logs/wo13_e2a/tier1_${label}.log
for step in $steps; do
  extra=()
  [ "$step" = gather ] && extra=(--fetchngs-dir "$raw/tier1_meta" --rnaseq-dir "$raw/tier1_rnaseq")
  [ "$step" = qc ] && extra=(--rnaseq-dir "$raw/tier1_rnaseq")
  echo "== $step $(date -u +%Y-%m-%dT%H:%M:%SZ) src=$PYTHONPATH" >> "$log"
  /usr/bin/time -l python -m ctx --root "$root" "$step" --dataset tier1 \
    --dataset-file "$WT/conf/datasets/tier1.toml" "${extra[@]}" >> "$log" 2>&1
  echo "exit $?" >> "$log"
done
grep -E "^== |^exit |real|peak memory" "$log" | tail -8
