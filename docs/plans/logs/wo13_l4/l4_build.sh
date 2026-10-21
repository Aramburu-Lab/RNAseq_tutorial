#!/bin/bash
# L.4: rebuild only M7 of ccRCC into a scratch root whose inputs link the real root.
#   bash l4_build.sh <label>
S=/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aba0cd8c07b2db4ce
real=/Users/rytis/ctx_kidney
root=$S/l4_ccrcc_$1
mkdir -p "$root/outputs"
[ -e "$root/data" ] || ln -s "$real/data" "$root/data"
[ -e "$root/resources" ] || ln -s "$real/resources" "$root/resources"
[ -e "$root/outputs/01_gather" ] || ln -s "$real/outputs/01_gather" "$root/outputs/01_gather"
[ -e "$root/outputs/05_de" ] || ln -s "$real/outputs/05_de" "$root/outputs/05_de"
touch "$S/l4_stamp_$1"
export PATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH
export PYTHONPATH=$WT/src PYTHONDONTWRITEBYTECODE=1
cd "$S" || exit 1
/usr/bin/time -l python -m ctx --root "$root" proteo --dataset cptac_ccrcc \
  --dataset-file "$WT/conf/datasets/cptac_ccrcc.toml" > "$S/l4_$1.log" 2>&1
echo "exit $?"
tail -4 "$S/l4_$1.log"
echo "written under the real root since the start:"
find "$real" -newer "$S/l4_stamp_$1" | head
