#!/bin/bash
# wo13-f1: Step-1 gather, de and function on a fresh scratch root with this worktree's src, then
# its 92 recorded metrics compared value by value (wo13_de/step1.sh's recipe, D-046). The real
# root's data/ and resources/ are linked, never written; every write lands in the scratch root.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a00e5034a0fa89898
PY=$WT/docs/plans/logs/wo13_f1_py.sh
root=$WT/outputs/scratch_roots/tier1_f1
raw=/Users/rytis/ctx_tier1/data/raw
if [ -e "$root" ]; then echo "$root exists: move it aside first"; exit 1; fi
mkdir -p "$root/outputs"
ln -s /Users/rytis/ctx_tier1/data "$root/data"
ln -s /Users/rytis/ctx_tier1/resources "$root/resources"
spec=$WT/conf/datasets/tier1.toml
for step in "gather --fetchngs-dir $raw/tier1_meta --rnaseq-dir $raw/tier1_rnaseq" de function; do
  # shellcheck disable=SC2086
  "$PY" python -m ctx --root "$root" $step --dataset tier1 --dataset-file "$spec" > /dev/null 2>&1
  echo "${step%% *}: exit $?"
done
out=$root/outputs
"$PY" python docs/plans/logs/wo13_p/compare_metrics.py \
  deliverables/2026-09-30_de_tier1_metrics_v1.csv \
  "$(ls "$out"/05_de/*_tier1_metrics_v1.csv | tail -1)" \
  "$(ls "$out"/06_function/*_tier1_metrics_v1.csv | tail -1)"
