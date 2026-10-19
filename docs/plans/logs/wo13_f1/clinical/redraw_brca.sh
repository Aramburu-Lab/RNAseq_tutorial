#!/bin/bash
# wo13-f1 M10: every M10 figure kind of the scratch BRCA page redrawn from its Reproduce panel
# (redraw.py takes the first panel of an ID, so ER is its OS curve); then the km_score page's
# live panel compared with the module's.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a43de516483292309
PY=$WT/docs/plans/logs/wo13_f1_clinical_py.sh
page=$(ls "$WT"/outputs/scratch_roots/brca_f1_head/outputs/12_clinical/*_tcga_brca_overview_v1.html)
clean=$WT/outputs/scratch_roots/redraw_m10
for id in breast_m10_er breast_m10_expression breast_m10_scan_score breast_m10_competing \
          breast_m10_power breast_m10_confounding; do
  echo "== $id"
  "$PY" python "$WT/docs/plans/logs/wo13_f1/redraw.py" "$page" "$id" "$clean/$id" 2>&1 \
    | grep "panels\|exit\|byte-equal\|DIFFERS\|Error"
done
echo "== live panel of km_score"
"$PY" python "$WT/docs/plans/logs/wo13_f1/clinical/live_alike.py" \
  "$(dirname "$page")"/*_km_score_v1.html "$clean"/breast_m10_km_score/repro_*/*_km_score_v1.html
