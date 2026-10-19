#!/usr/bin/env bash
# wo13-s2: ctx run --check on the real KIRC folders (read only), on tcga_kirc.toml as it is and
# on a copy without [clinical.published.table3_use], both in this worktree's scratch root.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aa11ca42a7ace51a9
R=/Users/rytis/ctx_kidney/data/raw
S=$WT/.wo13s2_roots/check_kirc
mkdir -p "$S/data/raw" "$S/resources"
# The resources --check looks for, linked from the real root, which ctx only reads.
for d in gdc_gencode_v36 gtex_v11 mirna_targets; do
  [ -e "$S/data/raw/$d" ] || ln -s "$R/$d" "$S/data/raw/$d"
done
[ -e "$S/resources/msigdb" ] || ln -s /Users/rytis/ctx_kidney/resources/msigdb "$S/resources/msigdb"
cp "$WT/conf/datasets/tcga_kirc.toml" "$S/tcga_kirc_as_is.toml"
python3 - "$WT/conf/datasets/tcga_kirc.toml" "$S/tcga_kirc_without.toml" <<'PY'
import sys
text = open(sys.argv[1]).read()
start = text.index("[clinical.published.table3_use]")
open(sys.argv[2], "w").write(text[:start] + text[text.index("\n\n", start) + 2:])
PY
for name in as_is without; do
  echo "=== tcga_kirc_$name.toml"
  bash "$WT/docs/plans/logs/wo13_s2_py.sh" ctx --root "$S" run --check --dataset tcga_kirc \
    --dataset-file "$S/tcga_kirc_$name.toml" --gdc-dir "$R/tcga_kirc" \
    --mirna-dir "$R/tcga_kirc_mirna" --isoform-dir "$R/tcga_kirc_isoform"
  echo "exit $?"
done
