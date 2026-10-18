#!/usr/bin/env bash
# wo13-v: build M10 for one cohort in a scratch root of this worktree. The real root's newest
# gather tables are linked file by file into a fresh outputs/01_gather, so ctx writes only here.
# Usage: build.sh <label> <real root> <slug>
set -euo pipefail
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ad2999f5c1b3365b1
label=$1
real=$2
slug=$3
root=$WT/.wo13v_roots/$label/$slug
rm -rf "$root"
mkdir -p "$root/outputs/01_gather"
for kind in samplesheet counts genes cdr cases molecular_tests treatments annotations; do
  newest=$(find "$real/outputs/01_gather" -name "*_gather_${slug}_${kind}_v*.parquet" | sort | tail -1)
  ln -s "$newest" "$root/outputs/01_gather/"
done
/usr/bin/time -l bash "$WT/docs/plans/logs/wo13_v/env.sh" ctx --root "$root" clinical \
  --dataset "$slug" --dataset-file "$WT/conf/datasets/$slug.toml" 2>&1 | tail -22
