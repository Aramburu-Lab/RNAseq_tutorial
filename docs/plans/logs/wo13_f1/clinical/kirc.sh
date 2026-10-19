#!/bin/bash
# wo13-f1 M10: TCGA-KIRC's M10 in a scratch root, gather's 2026-10-16 tables copied in, run with
# this worktree's src (as given) and with BASE's src (argument "base"); page sizes printed.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a43de516483292309
PY=$WT/docs/plans/logs/wo13_f1_clinical_py.sh
which=${1:-head}
root=$WT/outputs/scratch_roots/kirc_f1_$which
real=/Users/rytis/ctx_kidney
mkdir -p "$root/outputs/01_gather"
[ -e "$root/data" ] || ln -s "$real/data" "$root/data"
[ -e "$root/resources" ] || ln -s "$real/resources" "$root/resources"
for t in samplesheet counts genes cdr cases molecular_tests treatments annotations; do
  cp "$real/outputs/01_gather/2026-10-16_gather_tcga_kirc_${t}_v1".* "$root/outputs/01_gather/"
  cp "$real/outputs/01_gather/2026-10-16_gather_tcga_kirc_${t}_v1_README.md" "$root/outputs/01_gather/"
done
if [ "$which" = base ]; then
  env PYTHONPATH=/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/base_src/src \
    PATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH PYTHONDONTWRITEBYTECODE=1 \
    /usr/bin/time -l python -m ctx --root "$root" clinical --dataset tcga_kirc \
    --dataset-file "$WT/conf/datasets/tcga_kirc.toml" > "$root/run.log" 2>&1
else
  "$PY" /usr/bin/time -l python -m ctx --root "$root" clinical --dataset tcga_kirc \
    --dataset-file "$WT/conf/datasets/tcga_kirc.toml" > "$root/run.log" 2>&1
fi
echo "clinical ($which): exit $?"
grep "real\|maximum resident\|peak memory" "$root/run.log"
ls -l "$root"/outputs/12_clinical/*_overview_v1.html "$real"/outputs/12_clinical/2026-10-16_clinical_tcga_kirc_overview_v1.html | awk '{print $5, $9}'
ls -l "$root"/outputs/12_clinical/*_{figure_options,panel_values}_v1* 2>/dev/null | awk '{print $5, $9}'
