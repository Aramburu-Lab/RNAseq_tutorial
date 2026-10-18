#!/usr/bin/env bash
# R3: M7 on each recorded cohort in a scratch root that links the real root's inputs, then its
# metrics against the 2026-10-16 recorded run's, value by value. Nothing is written to the real
# roots. Usage: recorded_equal.sh <scratch folder>
set -euo pipefail
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ac5a02af0b00495c4
S="$1"
for pair in cptac_ccrcc:ctx_kidney cptac_brca:ctx_breast tier1:ctx_tier1; do
    slug=${pair%%:*}; real="$HOME/${pair##*:}"; root="$S/$slug"
    mkdir -p "$root/outputs"
    ln -sfn "$real/data" "$root/data"
    ln -sfn "$real/resources" "$root/resources"
    for d in 01_gather 05_de; do ln -sfn "$real/outputs/$d" "$root/outputs/$d"; done
    "$WT/docs/plans/logs/r3/env.sh" ctx --root "$root" proteo --dataset "$slug" \
        --dataset-file "$WT/conf/datasets/$slug.toml" > "$S/$slug.log" 2>&1
    "$WT/docs/plans/logs/r3/env.sh" python -c "
import sys, glob, polars as pl
new = pl.read_parquet(glob.glob('$root/outputs/07_proteo/*_proteo_${slug}_metrics_v1.parquet')[0])
old = pl.read_parquet('$real/outputs/07_proteo/2026-10-16_proteo_${slug}_metrics_v1.parquet')
j = old.join(new, on='metric', how='full', suffix='_new')
same = j.filter((pl.col('value') == pl.col('value_new')) | (pl.col('value').is_nan() & pl.col('value_new').is_nan()))
print('$slug', 'metrics', old.height, new.height, 'equal', same.height)
v = sorted(set(pl.read_parquet(glob.glob('$root/outputs/07_proteo/*_sanity_v1.parquet')[0])['rule_verdict']))
print('$slug', 'verdicts', v)
"
done
