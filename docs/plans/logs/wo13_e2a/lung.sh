#!/bin/bash
# Package wo13-e2a: ctx qc on TCGA lung as the real root ran it (--strata project --peers
# project), on a scratch root that links the real root's 2026-10-16 gather tables file by file.
#   [SRC=<another src>] bash lung.sh <label>
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a68d1b8f18deeebdd
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
export PATH=$BIN:/usr/bin:/bin PYTHONPATH=${SRC:-$WT/src} PYTHONDONTWRITEBYTECODE=1
label=$1
real=/Users/rytis/ctx_lung
root=$WT/outputs/scratch_roots/lung_$label
mkdir -p "$root/outputs/01_gather"
[ -e "$root/resources" ] || ln -s "$real/resources" "$root/resources"
for what in samplesheet counts genes; do
  f=2026-10-16_gather_tcga_lung_${what}_v1.parquet
  [ -e "$root/outputs/01_gather/$f" ] || ln -s "$real/outputs/01_gather/$f" "$root/outputs/01_gather/$f"
done
log=$WT/docs/plans/logs/wo13_e2a/lung_${label}.log
echo "== qc $(date -u +%Y-%m-%dT%H:%M:%SZ) src=$PYTHONPATH" >> "$log"
/usr/bin/time -l python -m ctx --root "$root" qc --dataset tcga_lung \
  --dataset-file "$WT/conf/datasets/tcga_lung.toml" --strata project --peers project >> "$log" 2>&1
echo "exit $?" >> "$log"
grep -E "^== |^exit |real|peak memory" "$log" | tail -4
ls -l "$root"/outputs/02_qc/*overview*
