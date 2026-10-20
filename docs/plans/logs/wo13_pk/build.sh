#!/bin/bash
# Package PK: build pages in a scratch root of this worktree from this worktree's src.
#   bash build.sh <label> <tier1|breast|ccrcc> <steps: gather de genecard ...>
# The root links the real root's data and each entry of its resources (read only by every
# step run here) and holds its own outputs. Gather's folder is linked from the real root
# unless gather is a step; with LINK_DE=1, so is M5's, for a genecard run on the real calls.
# genecard runs --offline: the response cache is read, never written.
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aa7a952e67debab5e
label=$1 which=$2
shift 2
export PATH=$BIN:$PATH PYTHONPATH=${SRC:-$WT/src} PYTHONDONTWRITEBYTECODE=1
case $which in
  tier1) real=/Users/rytis/ctx_tier1 dataset=tier1
         cards=(--genes ISG15,IFI6,STAT1,OAS2,CD274,AIM2,GBP5) gather=() ;;
  breast) real=/Users/rytis/ctx_breast dataset=tcga_brca cards=() gather=() ;;
  ccrcc) real=/Users/rytis/ctx_kidney dataset=cptac_ccrcc cards=()
         gather=(--gdc-dir $real/data/raw/cptac_ccrcc_star) ;;
  *) echo "tier1, breast or ccrcc"; exit 2 ;;
esac
root=$WT/outputs/scratch_roots/${which}_$label
mkdir -p "$root/outputs" "$root/resources"
[ -e "$root/data" ] || ln -s "$real/data" "$root/data"
for entry in "$real"/resources/*; do
  [ -e "$root/resources/${entry##*/}" ] || ln -s "$entry" "$root/resources/${entry##*/}"
done
if [[ " $* " != *" gather "* ]] && [ ! -e "$root/outputs/01_gather" ]; then
  ln -s "$real/outputs/01_gather" "$root/outputs/01_gather"
fi
if [ "${LINK_DE:-0}" = 1 ] && [ ! -e "$root/outputs/05_de" ]; then
  ln -s "$real/outputs/05_de" "$root/outputs/05_de"
fi
file=$WT/conf/datasets/$dataset.toml
log=$WT/docs/plans/logs/wo13_pk/build_${which}_$label.log
for step in "$@"; do
  extra=()
  [ "$step" = gather ] && extra=("${gather[@]}")
  [ "$step" = genecard ] && extra=(--offline "${cards[@]}")
  echo "== $step $(date -u +%Y-%m-%dT%H:%M:%SZ) $(git -C "$WT" rev-parse --short HEAD)" >> "$log"
  /usr/bin/time -l python -m ctx --root "$root" "$step" --dataset "$dataset" \
    --dataset-file "$file" "${extra[@]}" >> "$log" 2>&1
  echo "exit $?" >> "$log"
done
grep -E "^== |exit |real|maximum resident" "$log" | tail -12
