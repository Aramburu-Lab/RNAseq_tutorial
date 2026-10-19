#!/bin/bash
# wo13-f1: Step-1's M4 in the scratch root step1.sh built, with the outdir and copy Q1 the
# 2026-10-16 run of ~/ctx_tier1 read; then both figures with a spec redrawn from their panels.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a00e5034a0fa89898
PY=$WT/docs/plans/logs/wo13_f1_py.sh
root=$WT/outputs/scratch_roots/tier1_f1
clean=$WT/outputs/scratch_roots/redraw_m4
raw=/Users/rytis/ctx_tier1/data/raw
"$PY" python -m ctx --root "$root" quant --dataset tier1 --dataset-file "$WT/conf/datasets/tier1.toml" \
  --rnaseq-dir "$raw/tier1_rnaseq" --salmon-dir "$raw/tier1_salmon" > /dev/null 2>&1
echo "quant: exit $?"
page=$(ls "$root"/outputs/04_quant/*_tier1_overview_v1.html)
ls -l "$page" /Users/rytis/ctx_tier1/outputs/04_quant/2026-10-16_quant_tier1_overview_v1.html | awk '{print $5, $9}'
mkdir -p "$clean"
for id in blood_m4_fates blood_m4_estimates; do
  echo "== $id"
  "$PY" python "$WT/docs/plans/logs/wo13_f1/redraw.py" "$page" "$id" "$clean/$id" 2>&1 | grep "panels\|exit\|byte-equal\|DIFFERS\|Error"
done
