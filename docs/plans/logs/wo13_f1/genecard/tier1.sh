#!/bin/bash
# wo13-f1 M9: Step-1 gather, de and genecard on a fresh scratch root with this worktree's src.
# data/ is linked read only; resources/ is a real folder holding a copy of the response cache,
# which genecard writes, and links to the read-only GTF, dataset and gene-set folders. The
# genes are the 2026-10-16 run's, read --offline from the copied cache.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aafb5beef9062b91f
PY=$WT/docs/plans/logs/wo13_f1_genecard_py.sh
real=/Users/rytis/ctx_tier1
root=$WT/outputs/scratch_roots/tier1_m9
raw=$real/data/raw
if [ -e "$root" ]; then echo "$root exists: move it aside first"; exit 1; fi
mkdir -p "$root/outputs" "$root/resources"
ln -s "$real/data" "$root/data"
cp -R "$real/resources/cache" "$root/resources/cache"
cp "$real/resources/manifest.json" "$root/resources/manifest.json"
for f in gencode datasets msigdb; do ln -s "$real/resources/$f" "$root/resources/$f"; done
spec=$WT/conf/datasets/tier1.toml
for step in "gather --fetchngs-dir $raw/tier1_meta --rnaseq-dir $raw/tier1_rnaseq" de; do
  # shellcheck disable=SC2086
  "$PY" python -m ctx --root "$root" $step --dataset tier1 --dataset-file "$spec" > /dev/null 2>&1
  echo "${step%% *}: exit $?"
done
"$PY" /usr/bin/time -l python -m ctx --root "$root" genecard --dataset tier1 \
  --dataset-file "$spec" --offline --genes ISG15,IFI6,STAT1,OAS2,CD274,AIM2,GBP5
echo "genecard: exit $?"
