#!/bin/bash
# wo13-y: Step-1's M2 and TCGA-BRCA's M10 with main's code (an export of c33b4ae's src, whose
# qc and clinical code is main's), each in its own scratch root, to tell whether the pages over
# the 15 MB cap are this package's or main's. tier1_y's gather and brca_y's copied gather
# tables are linked in, read only.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a93c44dfd6e2a237b
SRC=$WT/outputs/scratch_y/main_src/src
export PATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH
export PYTHONPATH=$SRC PYTHONDONTWRITEBYTECODE=1 CTX_CHROMIUM=/Applications/Chromium.app/Contents/MacOS/Chromium
python -c "import ctx; print('ctx from', ctx.__file__)"
raw=/Users/rytis/ctx_tier1/data/raw
for which in tier1 brca; do
  root=$WT/outputs/scratch_roots/${which}_main
  if [ -e "$root" ]; then echo "$root exists: move it aside first"; exit 1; fi
  mkdir -p "$root/outputs"
done
root=$WT/outputs/scratch_roots/tier1_main
ln -s /Users/rytis/ctx_tier1/data "$root/data"
ln -s /Users/rytis/ctx_tier1/resources "$root/resources"
ln -s "$WT/outputs/scratch_roots/tier1_y/outputs/01_gather" "$root/outputs/01_gather"
python -m ctx --root "$root" qc --dataset tier1 --dataset-file "$WT/conf/datasets/tier1.toml" \
  --rnaseq-dir "$raw/tier1_rnaseq" > "$root/qc.log" 2>&1
echo "qc (main): exit $?"; grep -o "[0-9]* bytes, over the" "$root/qc.log"
root=$WT/outputs/scratch_roots/brca_main
ln -s /Users/rytis/ctx_breast/data "$root/data"
ln -s /Users/rytis/ctx_breast/resources "$root/resources"
ln -s "$WT/outputs/scratch_roots/brca_y/outputs/01_gather" "$root/outputs/01_gather"
python -m ctx --root "$root" clinical --dataset tcga_brca \
  --dataset-file "$WT/conf/datasets/tcga_brca.toml" > "$root/run.log" 2>&1
echo "clinical (main): exit $?"; grep -o "[0-9]* bytes, over the" "$root/run.log"
ls -l "$root"/outputs/12_clinical/*_overview_v1.html 2>/dev/null | awk '{print $5, $9}'
