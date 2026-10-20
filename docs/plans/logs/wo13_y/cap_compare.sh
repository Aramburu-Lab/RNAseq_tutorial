#!/bin/bash
# wo13-y: Step-1's M2 and TCGA-BRCA's M10 with the cap lifted (uncapped.py), from two copies of
# src outside git, main's (c33b4ae's export) and this package's (a copy of the worktree's), so
# the Reproduce panels of both name no commit and the pages differ by the code alone.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a93c44dfd6e2a237b
export PATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH
export PYTHONDONTWRITEBYTECODE=1 CTX_QC_REGISTRY=$WT/conf/qc_metrics.toml
rm -rf "$WT/outputs/scratch_y/y_src" && mkdir -p "$WT/outputs/scratch_y/y_src" && cp -R "$WT/src" "$WT/outputs/scratch_y/y_src/src"
find "$WT/outputs/scratch_y/y_src" -name __pycache__ -prune -exec rm -rf {} +
for side in main y; do
  export PYTHONPATH=$WT/outputs/scratch_y/${side}_src/src
  for which in tier1 brca; do
    root=$WT/outputs/scratch_roots/${which}_cap_$side
    if [ -e "$root" ]; then echo "$root exists: move it aside first"; exit 1; fi
    mkdir -p "$root/outputs"
    case $which in
      tier1) real=/Users/rytis/ctx_tier1 gather=$WT/outputs/scratch_roots/tier1_y/outputs/01_gather
             args=(qc --dataset tier1 --dataset-file "$WT/conf/datasets/tier1.toml"
                   --rnaseq-dir "$real/data/raw/tier1_rnaseq") ;;
      brca) real=/Users/rytis/ctx_breast gather=$WT/outputs/scratch_roots/brca_y/outputs/01_gather
            args=(clinical --dataset tcga_brca --dataset-file "$WT/conf/datasets/tcga_brca.toml") ;;
    esac
    ln -s "$real/data" "$root/data"
    ln -s "$real/resources" "$root/resources"
    ln -s "$gather" "$root/outputs/01_gather"
    python "$WT/docs/plans/logs/wo13_y/uncapped.py" --root "$root" "${args[@]}" > "$root/run.log" 2>&1
    echo "$side $which: exit $?; $(head -1 "$root/run.log")"
    for page in "$root"/outputs/*/*_overview_v1.html; do echo "  $(stat -f %z "$page") B $page"; done
  done
done
