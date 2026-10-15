#!/bin/bash
# Which commits hold each GDC dataset file's bytes, against what the gather READMEs recorded.
cd "$(dirname "$0")/../../../.." || exit 1
for n in tcga_kirc cptac_ccrcc tcga_brca cptac_brca; do
  echo "HEAD $n $(git show HEAD:conf/datasets/$n.toml | shasum -a 256 | cut -c1-64)"
  for c in $(git log --format=%h -- conf/datasets/$n.toml); do
    echo "  $c $(git show "$c:conf/datasets/$n.toml" 2>/dev/null | shasum -a 256 | cut -c1-64)"
  done
done
for f in ~/ctx_kidney/outputs/01_gather/*_counts_v1_README.md ~/ctx_breast/outputs/01_gather/*_counts_v1_README.md; do
  echo "$(basename "$f") $(grep -o 'sha256:[a-z_]*\.toml | [0-9a-f]*' "$f" | head -1)"
done
