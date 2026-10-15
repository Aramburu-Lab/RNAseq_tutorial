#!/bin/bash
# Package P: build Step-1 (tier1) and TCGA-BRCA pages in a scratch root of this worktree.
#   bash build.sh <label> <src dir> <tier1|breast> <steps: gather de function report ...>
# The root links the real data and resources (read only by every step run here) and holds its
# own outputs; a label's 01_gather can be linked from the "base" root to skip gather.
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a4bc57162a0d35df8
label=$1 src=$2 which=$3
shift 3
export PATH=$BIN:$PATH PYTHONPATH=$src PYTHONDONTWRITEBYTECODE=1
case $which in
  tier1) real=/Users/rytis/ctx_tier1 dataset=tier1
         gather=(--fetchngs-dir $real/data/raw/tier1_meta --rnaseq-dir $real/data/raw/tier1_rnaseq) ;;
  breast) real=/Users/rytis/ctx_breast dataset=tcga_brca
          gather=(--gdc-dir $real/data/raw/tcga_brca --mirna-dir $real/data/raw/tcga_brca_mirna
                  --isoform-dir $real/data/raw/tcga_brca_isoform) ;;
  lusc) real=/Users/rytis/ctx_lung dataset=tcga_lusc
        gather=(--gdc-dir $real/data/raw/tcga_lusc --mirna-dir $real/data/raw/tcga_lusc_mirna
                --isoform-dir $real/data/raw/tcga_lusc_isoform) ;;
  ccrcc) real=/Users/rytis/ctx_kidney dataset=cptac_ccrcc
         gather=(--gdc-dir $real/data/raw/cptac_ccrcc_star) ;;
  *) echo "tier1, breast, lusc or ccrcc"; exit 2 ;;
esac
root=$WT/outputs/scratch_roots/${which}_$label
mkdir -p "$root/outputs"
[ -e "$root/data" ] || ln -s "$real/data" "$root/data"
[ -e "$root/resources" ] || ln -s "$real/resources" "$root/resources"
base=$WT/outputs/scratch_roots/${which}_base/outputs/01_gather
if [ "$label" != base ] && [ ! -e "$root/outputs/01_gather" ] && [ -d "$base" ]; then
  ln -s "$base" "$root/outputs/01_gather"
fi
file=$WT/conf/datasets/$dataset.toml
log=$WT/docs/plans/logs/wo13_p/build_${which}_$label.log
for step in "$@"; do
  extra=()
  [ "$step" = gather ] && extra=("${gather[@]}")
  echo "== $step $(date -u +%Y-%m-%dT%H:%M:%SZ) src=$src" >> "$log"
  /usr/bin/time -l python -m ctx --root "$root" "$step" --dataset "$dataset" \
    --dataset-file "$file" "${extra[@]}" >> "$log" 2>&1
  echo "exit $?" >> "$log"
done
grep -E "^== |exit |real|maximum resident" "$log" | tail -12
