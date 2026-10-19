#!/usr/bin/env bash
# wo13-s2: build M10 for one cohort in a scratch root of this worktree. The real root's newest
# gather tables are linked file by file into a fresh outputs/01_gather, so ctx writes only here.
# A root that exists is refused, never removed. Usage: build_m10.sh <label> <real root> <slug>
set -euo pipefail
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aa11ca42a7ace51a9
label=$1
real=$2
slug=$3
root=$WT/.wo13s2_roots/$label/$slug
if [ -e "$root" ]; then echo "exists: $root"; exit 1; fi
mkdir -p "$root/outputs/01_gather"
for kind in samplesheet counts genes cdr cases molecular_tests treatments annotations; do
  newest=$(find "$real/outputs/01_gather" -name "*_gather_${slug}_${kind}_v*.parquet" | sort | tail -1)
  ln -s "$newest" "$root/outputs/01_gather/"
done
/usr/bin/time -l bash "$WT/docs/plans/logs/wo13_s2_py.sh" ctx --root "$root" clinical \
  --dataset "$slug" --dataset-file "$WT/conf/datasets/$slug.toml" 2>&1 | tail -22
