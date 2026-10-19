#!/bin/bash
# F.1 M8: tcga_brca's M8 in the scratch root breast_root.sh built, with the command the
# 2026-10-16 page of ~/ctx_breast prints (its dataset file this worktree's copy); then one
# figure of each half redrawn from its panel in a clean folder, and the page's panel bytes.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a84b7900e11d976d1
PY=$WT/docs/plans/logs/wo13_f1_ncrna_py.sh
root=$WT/outputs/scratch_roots/breast_m8
clean=$WT/outputs/scratch_roots/redraw_m8
/usr/bin/time -l "$PY" python -m ctx --root "$root" ncrna --dataset tcga_brca \
  --dataset-file "$WT/conf/datasets/tcga_brca.toml" > "$WT/docs/plans/logs/wo13_f1/ncrna/breast_run.out" 2> "$WT/docs/plans/logs/wo13_f1/ncrna/breast_run.err"
echo "ncrna: exit $?"
grep "real\|maximum resident\|peak memory" "$WT/docs/plans/logs/wo13_f1/ncrna/breast_run.err"
page=$(ls "$root"/outputs/08_ncrna/*_tcga_brca_overview_v1.html)
"$PY" python "$WT/docs/plans/logs/wo13_f1/panel_bytes.py" "$page" \
  /Users/rytis/ctx_breast/outputs/08_ncrna/2026-10-16_ncrna_tcga_brca_overview_v1.html
ls -l "$root"/outputs/08_ncrna/*_figure_options_v1* "$root"/outputs/08_ncrna/*_gene_biotypes_v1* | awk '{print $5, $9}'
mkdir -p "$clean"
for id in "$@"; do
  echo "== $id"
  "$PY" python "$WT/docs/plans/logs/wo13_f1/redraw.py" "$page" "$id" "$clean/$id" 2>&1 | grep "panels\|exit\|byte-equal\|DIFFERS\|Error"
done
