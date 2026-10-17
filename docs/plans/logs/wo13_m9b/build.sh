#!/bin/bash
# Package M9b: build pages in a scratch root of this worktree from this worktree's src.
#   bash build.sh <label> <tier1|breast> <steps: de genecard ...>
# The root links the real root's data, each entry of its resources and its gather folder, all
# read only by the steps run here, and holds its own outputs; with DE_FROM=<label>, M5's folder
# is linked from that scratch root. genecard runs --offline: the cache is read, never written.
. "$(dirname "$0")/env.sh"
label=$1 which=$2
shift 2
case $which in
  tier1) real=/Users/rytis/ctx_tier1 dataset=tier1 ;;
  breast) real=/Users/rytis/ctx_breast dataset=tcga_brca ;;
  *) echo "tier1 or breast"; exit 2 ;;
esac
root=$WT/outputs/scratch_roots/m9b_${which}_$label
mkdir -p "$root/outputs" "$root/resources"
[ -e "$root/data" ] || ln -s "$real/data" "$root/data"
for entry in "$real"/resources/*; do
  [ -e "$root/resources/${entry##*/}" ] || ln -s "$entry" "$root/resources/${entry##*/}"
done
[ -e "$root/outputs/01_gather" ] || ln -s "$real/outputs/01_gather" "$root/outputs/01_gather"
if [ -n "${DE_FROM:-}" ] && [ ! -e "$root/outputs/05_de" ]; then
  ln -s "$WT/outputs/scratch_roots/m9b_${which}_$DE_FROM/outputs/05_de" "$root/outputs/05_de"
fi
file=$WT/conf/datasets/$dataset.toml
log=$WT/docs/plans/logs/wo13_m9b/build_${which}_$label.log
for step in "$@"; do
  extra=()
  [ "$step" = genecard ] && extra=(--offline ${GENES:+--genes "$GENES"})
  echo "== $step $(date -u +%Y-%m-%dT%H:%M:%SZ) $(git -C "$WT" rev-parse --short HEAD) $(uptime | sed 's/.*load/load/')" >> "$log"
  /usr/bin/time -l python -m ctx --root "$root" "$step" --dataset "$dataset" \
    --dataset-file "$file" "${extra[@]}" >> "$log" 2>&1
  echo "exit $?" >> "$log"
done
grep -E "^== |exit |real|maximum resident|peak memory" "$log" | tail -12
