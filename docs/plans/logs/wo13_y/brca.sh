#!/bin/bash
# wo13-y: TCGA-BRCA's M10 in a scratch root, gather's 2026-10-16 tables copied in, run with this
# worktree's src (wo13_f1/clinical/brca.sh's recipe); its page size beside the real root's.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a93c44dfd6e2a237b
PY=$WT/docs/plans/logs/wo13_y_py.sh
root=$WT/outputs/scratch_roots/brca_y
real=/Users/rytis/ctx_breast
if [ -e "$root" ]; then echo "$root exists: move it aside first"; exit 1; fi
mkdir -p "$root/outputs/01_gather"
ln -s "$real/data" "$root/data"
ln -s "$real/resources" "$root/resources"
for t in samplesheet counts genes cdr cases molecular_tests treatments annotations; do
  cp "$real/outputs/01_gather/2026-10-16_gather_tcga_brca_${t}_v1".* "$root/outputs/01_gather/"
  cp "$real/outputs/01_gather/2026-10-16_gather_tcga_brca_${t}_v1_README.md" "$root/outputs/01_gather/"
done
"$PY" /usr/bin/time -l python -m ctx --root "$root" clinical --dataset tcga_brca \
  --dataset-file "$WT/conf/datasets/tcga_brca.toml" > "$root/run.log" 2>&1
echo "clinical: exit $?"
for page in "$root"/outputs/12_clinical/*_overview_v1.html \
    "$real"/outputs/12_clinical/2026-10-16_clinical_tcga_brca_overview_v1.html; do
  echo "$(stat -f %z "$page") B $page"
done
