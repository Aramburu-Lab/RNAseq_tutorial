#!/bin/bash
# wo13-f1: Step-1's M2 in the scratch root step1.sh built (gather there), with --rnaseq-dir as
# the 2026-10-16 run of ~/ctx_tier1 took it; then three figures redrawn from their panels.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a00e5034a0fa89898
PY=$WT/docs/plans/logs/wo13_f1_py.sh
root=$WT/outputs/scratch_roots/tier1_f1
clean=$WT/outputs/scratch_roots/redraw_m2
"$PY" python -m ctx --root "$root" qc --dataset tier1 --dataset-file "$WT/conf/datasets/tier1.toml" \
  --rnaseq-dir /Users/rytis/ctx_tier1/data/raw/tier1_rnaseq > /dev/null 2>&1
echo "qc: exit $?"
page=$(ls "$root"/outputs/02_qc/*_tier1_overview_v1.html)
ls -l "$page" /Users/rytis/ctx_tier1/outputs/02_qc/2026-10-16_qc_tier1_overview_v1.html | awk '{print $5, $9}'
mkdir -p "$clean"
for id in blood_m2_strip blood_m2_fastqc blood_m2_linked_pca; do
  echo "== $id"
  "$PY" python "$WT/docs/plans/logs/wo13_f1/redraw.py" "$page" "$id" "$clean/$id" 2>&1 | grep "panels\|exit\|byte-equal\|DIFFERS\|Error"
done
