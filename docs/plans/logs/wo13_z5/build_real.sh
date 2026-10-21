#!/bin/bash
# wo13-z5: real roots' pages rebuilt into scratch roots, before (f287e07's src) and after (this
# worktree's), nothing written into the real root. Usage:
#   bash build_real.sh <before|after> <real root> <dataset> <scratch parent> <src> <steps...>
# The scratch root links the real root's data, each entry of its resources and its gather
# folder, read only by the steps run here; each other step's folder is linked when a step
# needs it (DE_LINK=1 links 05_de for function and later).
which=$1 real=$2 dataset=$3 parent=$4 src=$5
shift 5
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a6b39a2692f007e36
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
export PATH=$BIN:/usr/bin:/bin:/usr/sbin:/sbin PYTHONPATH=$src PYTHONDONTWRITEBYTECODE=1
root=$parent/${dataset}_$which
mkdir -p "$root/outputs" "$root/resources"
[ -e "$root/data" ] || ln -s "$real/data" "$root/data"
for entry in "$real"/resources/*; do
  [ -e "$root/resources/${entry##*/}" ] || ln -s "$entry" "$root/resources/${entry##*/}"
done
[ -e "$root/outputs/01_gather" ] || ln -s "$real/outputs/01_gather" "$root/outputs/01_gather"
if [ -n "${DE_LINK:-}" ] && [ ! -e "$root/outputs/05_de" ]; then
  ln -s "$real/outputs/05_de" "$root/outputs/05_de"
fi
file=$WT/conf/datasets/$dataset.toml
raw=$real/data/raw
log=$parent/build_${dataset}_$which.log
cd "$root" || exit 1
for step in "$@"; do
  extra=()
  [ "$step" = genecard ] && extra=(--offline)
  [ "$step" = qc ] && [ "$dataset" = tier1 ] && extra=(--rnaseq-dir "$raw/tier1_rnaseq")
  [ "$step" = quant ] && extra=(--rnaseq-dir "$raw/tier1_rnaseq" --salmon-dir "$raw/tier1_salmon")
  echo "== $step $which $(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$log"
  /usr/bin/time -l python -m ctx --root "$root" "$step" --dataset "$dataset" \
    --dataset-file "$file" "${extra[@]}" >> "$log" 2>&1
  echo "exit $?" >> "$log"
done
grep -E "^== |^exit |real|maximum resident" "$log" | tail -12
