#!/bin/bash
# wo13-f1 M9: genecard once more on tier1_m9 after the family table was renamed own_family
# (box_family put "_box_" in a table name, which test_genecard_own reads as a box figure).
# The earlier run is moved aside; then the junctions panel is redrawn from the new page.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aafb5beef9062b91f
PY=$WT/docs/plans/logs/wo13_f1_genecard_py.sh
root=$WT/outputs/scratch_roots/tier1_m9
aside=$WT/outputs/scratch_roots/tier1_m9_aside
mv "$root/outputs/11_genecard" "$aside/11_genecard_box_family"
"$PY" /usr/bin/time -l python -m ctx --root "$root" genecard --dataset tier1 \
  --dataset-file "$WT/conf/datasets/tier1.toml" --offline \
  --genes ISG15,IFI6,STAT1,OAS2,CD274,AIM2,GBP5 > /dev/null
echo "genecard: exit $?"
ls -l "$root"/outputs/11_genecard/*_overview_v1.html
page=$(ls "$root"/outputs/11_genecard/*_overview_v1.html)
"$PY" python docs/plans/logs/wo13_f1/redraw.py "$page" blood_m9_junctions \
  "$WT/outputs/scratch_roots/redraw_tier1_junctions"
