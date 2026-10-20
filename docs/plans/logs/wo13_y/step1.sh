#!/bin/bash
# wo13-y: Step-1 gather, de and function on a fresh scratch root with this worktree's src, its
# 92 recorded metrics compared value by value (wo13_f1/step1.sh's recipe, D-046); then qc and
# quant, as the 2026-10-16 run of ~/ctx_tier1 took them, for the M2 and M4 page sizes. The real
# root's data/ and resources/ are linked, never written; every write lands in the scratch root.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a93c44dfd6e2a237b
PY=$WT/docs/plans/logs/wo13_y_py.sh
root=$WT/outputs/scratch_roots/tier1_y
raw=/Users/rytis/ctx_tier1/data/raw
if [ -e "$root" ]; then echo "$root exists: move it aside first"; exit 1; fi
mkdir -p "$root/outputs"
ln -s /Users/rytis/ctx_tier1/data "$root/data"
ln -s /Users/rytis/ctx_tier1/resources "$root/resources"
spec=$WT/conf/datasets/tier1.toml
for step in "gather --fetchngs-dir $raw/tier1_meta --rnaseq-dir $raw/tier1_rnaseq" de function \
    "qc --rnaseq-dir $raw/tier1_rnaseq" \
    "quant --rnaseq-dir $raw/tier1_rnaseq --salmon-dir $raw/tier1_salmon"; do
  # shellcheck disable=SC2086
  "$PY" python -m ctx --root "$root" $step --dataset tier1 --dataset-file "$spec" \
    > "$root/${step%% *}.log" 2>&1
  echo "${step%% *}: exit $?"
done
out=$root/outputs
"$PY" python docs/plans/logs/wo13_p/compare_metrics.py \
  deliverables/2026-09-30_de_tier1_metrics_v1.csv \
  "$(ls "$out"/05_de/*_tier1_metrics_v1.csv | tail -1)" \
  "$(ls "$out"/06_function/*_tier1_metrics_v1.csv | tail -1)"
for page in "$out"/02_qc/*_tier1_overview_v1.html "$out"/04_quant/*_tier1_overview_v1.html \
    /Users/rytis/ctx_tier1/outputs/02_qc/2026-10-16_qc_tier1_overview_v1.html \
    /Users/rytis/ctx_tier1/outputs/04_quant/2026-10-16_quant_tier1_overview_v1.html; do
  echo "$(stat -f %z "$page") B $page"
done
