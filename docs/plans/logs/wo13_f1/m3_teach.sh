#!/bin/bash
# wo13-f1: Step-1's teaching M3 (tier1_teach, HISAT2, as ~/ctx_tier1's 2026-10-07 run) in the
# scratch root step1.sh built, gather's tier1_teach sheet copied in from the real root; then
# two figures redrawn from their panels.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a00e5034a0fa89898
PY=$WT/docs/plans/logs/wo13_f1_py.sh
root=$WT/outputs/scratch_roots/tier1_f1
clean=$WT/outputs/scratch_roots/redraw_m3
raw=/Users/rytis/ctx_tier1/data/raw
for f in /Users/rytis/ctx_tier1/outputs/01_gather/2026-10-07_gather_tier1_teach_samplesheet_v1*; do
  cp "$f" "$root/outputs/01_gather/"
done
/usr/bin/time -l "$PY" python -m ctx --root "$root" process --dataset tier1_teach \
  --dataset-file "$WT/conf/datasets/tier1_teach.toml" --genome "$raw/reference/GRCh38.primary_assembly.chr22.fa.gz" \
  --aligner hisat2 --compare-dataset tier1 --fetchngs-dir "$raw/tier1_teach" > "$clean.process.log" 2>&1
echo "process: exit $?"
tail -20 "$clean.process.log" | grep "real\|maximum resident\|peak memory\|Error"
page=$(ls "$root"/outputs/03_process/*_tier1_teach_overview_v1.html)
ls -l "$page" /Users/rytis/ctx_tier1/outputs/03_process/2026-10-07_process_tier1_teach_overview_v1.html | awk '{print $5, $9}'
mkdir -p "$clean"
for id in blood_m3_mapping blood_m3_counters; do
  echo "== $id"
  "$PY" python "$WT/docs/plans/logs/wo13_f1/redraw.py" "$page" "$id" "$clean/$id" 2>&1 | grep "panels\|exit\|byte-equal\|DIFFERS\|Error"
done
