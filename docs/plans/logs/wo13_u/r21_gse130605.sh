#!/bin/bash
# wo13-u R2.1: GSE130605 again in a fresh root, R2's dataset file plus the hgnc_complete
# resource of conf/dataset_template_counts.toml; R2's root (~/ctx_r2_gse130605) is read, never
# written. Its resources and inputs are copied, so ctx setup downloads only HGNC's file
# (16,973,125 B, stated before the run).
set -euo pipefail
PY=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aa1d340947c025c9d/docs/plans/logs/wo13_u_py.sh
TEMPLATE=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aa1d340947c025c9d/conf/dataset_template_counts.toml
OLD=~/ctx_r2_gse130605
NEW=~/ctx_u_gse130605
SLUG=r2_gse130605
if [ -e "$NEW" ]; then echo "$NEW exists: move it aside first"; exit 1; fi
mkdir -p "$NEW/data/raw/$SLUG" "$NEW/resources"
cp -R "$OLD/resources/." "$NEW/resources/"
cp "$OLD/data/raw/$SLUG/"* "$NEW/data/raw/$SLUG/"
cp "$OLD/$SLUG.csv" "$NEW/"
# The template's hgnc_complete block, inserted before [sdrf].
python3 - "$TEMPLATE" "$OLD/$SLUG.toml" "$NEW/$SLUG.toml" <<'EOF'
import re, sys
template, old, new = (open(p).read() if i < 2 else p for i, p in enumerate(sys.argv[1:]))
block = re.search(r"\n# HGNC's complete set.*?\ndest = \"resources/hgnc/[^\n]*\n", template, re.S)[0]
assert old.count("\n[sdrf]") == 1
open(new, "w").write(old.replace("\n[sdrf]", block + "\n[sdrf]"))
EOF
MATRIX="$NEW/data/raw/$SLUG/GSE130605_raw_counts_GRCh38.p13_NCBI.tsv.gz"
"$PY" python -m ctx --root "$NEW" setup --dataset-file "$NEW/$SLUG.toml" | tail -3
"$PY" python -m ctx --root "$NEW" run --dataset "$SLUG" --dataset-file "$NEW/$SLUG.toml" --counts "$MATRIX" --check | tail -2
/usr/bin/time -l "$PY" python -m ctx --root "$NEW" run --dataset "$SLUG" --dataset-file "$NEW/$SLUG.toml" --counts "$MATRIX" 2>&1 | grep -v "^/Users" | tail -25
