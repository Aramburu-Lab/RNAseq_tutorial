#!/bin/bash
# wo13-f1 M9: tier1.sh's genecard stopped in F1's views (an OwnCaption where export.settings
# takes a string, on BASE's code path). Its first attempt is moved aside, the scratch de's size
# factors moved out of 05_de, and genecard run again: the own panels then say why they are
# absent, as on the real root's 2026-10-16 run, which had no size factors either.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aafb5beef9062b91f
PY=$WT/docs/plans/logs/wo13_f1_genecard_py.sh
root=$WT/outputs/scratch_roots/tier1_m9
aside=$WT/outputs/scratch_roots/tier1_m9_aside
mkdir -p "$aside"
[ -d "$root/outputs/11_genecard" ] && mv "$root/outputs/11_genecard" "$aside/11_genecard_first"
mv "$root"/outputs/05_de/*_tier1_size_factors_v1* "$aside/"
spec=$WT/conf/datasets/tier1.toml
"$PY" /usr/bin/time -l python -m ctx --root "$root" genecard --dataset tier1 \
  --dataset-file "$spec" --offline --genes ISG15,IFI6,STAT1,OAS2,CD274,AIM2,GBP5
echo "genecard: exit $?"
