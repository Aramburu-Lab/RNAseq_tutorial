#!/usr/bin/env bash
# wo13-v: build M10 for the five cohorts in scratch roots under one label.
# Usage: build_all.sh <label>
L=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ad2999f5c1b3365b1/docs/plans/logs/wo13_v
while read -r real slug; do
  echo "=== $slug"
  bash "$L/build.sh" "$1" "$real" "$slug"
done <<'EOF'
/Users/rytis/ctx_breast tcga_brca
/Users/rytis/ctx_kidney tcga_kirc
/Users/rytis/ctx_lung tcga_luad
/Users/rytis/ctx_lung tcga_lusc
/Users/rytis/ctx_tcga_hnsc tcga_hnsc
EOF
