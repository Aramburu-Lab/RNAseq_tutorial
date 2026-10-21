#!/bin/bash
# wo13-z1: one cohort's M10 into a scratch root, run with a given src; nothing written into the
# real root. Usage: build.sh <real root> <slug> <src dir> <scratch root>. The real root's data/
# and resources/ are linked; gather's newest table of each kind for the slug is copied in.
real=$1 slug=$2 src=$3 root=$4
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-af769670337776f4c
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
if [ -e "$root" ]; then echo "$root exists: move it aside first"; exit 1; fi
mkdir -p "$root/outputs/01_gather"
ln -s "$real/data" "$root/data"
ln -s "$real/resources" "$root/resources"
for t in samplesheet counts genes cdr cases molecular_tests treatments annotations; do
  newest=$(ls "$real/outputs/01_gather/" | grep -E "^[0-9-]+_gather_${slug}_${t}_v1\.(parquet|csv)$" | sort | tail -1)
  day=${newest%%_gather_*}
  cp "$real/outputs/01_gather/${day}_gather_${slug}_${t}_v1".* "$root/outputs/01_gather/"
  cp "$real/outputs/01_gather/${day}_gather_${slug}_${t}_v1_README.md" "$root/outputs/01_gather/" 2>/dev/null
done
export PYTHONPATH=$src PYTHONDONTWRITEBYTECODE=1 PATH=$BIN:$PATH
cd "$root" || exit 1
/usr/bin/time -l "$BIN/python" -m ctx --root "$root" clinical --dataset "$slug" \
  --dataset-file "$WT/conf/datasets/$slug.toml" > "$root/run.log" 2>&1
echo "clinical $slug: exit $?"
for page in "$root"/outputs/12_clinical/*_overview_v1.html; do
  echo "$(stat -f %z "$page") B $page"
done
