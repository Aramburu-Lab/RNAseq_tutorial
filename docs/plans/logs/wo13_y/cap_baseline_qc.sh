#!/bin/bash
# wo13-y: cap_baseline.sh's qc run again with the QC registry named: an export of main's src has
# no conf/ beside it, so qc found no registry there. Same scratch root, same code.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a93c44dfd6e2a237b
export PATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH
export PYTHONPATH=$WT/outputs/scratch_y/main_src/src PYTHONDONTWRITEBYTECODE=1
export CTX_QC_REGISTRY=$WT/conf/qc_metrics.toml
root=$WT/outputs/scratch_roots/tier1_main
python -c "import ctx; print('ctx from', ctx.__file__)"
python -m ctx --root "$root" qc --dataset tier1 --dataset-file "$WT/conf/datasets/tier1.toml" \
  --rnaseq-dir /Users/rytis/ctx_tier1/data/raw/tier1_rnaseq > "$root/qc2.log" 2>&1
echo "qc (main): exit $?"; grep -A1 "bytes, over" "$root/qc2.log"
ls -l "$root"/outputs/02_qc/*_overview_v1.html 2>/dev/null | awk '{print $5, $9}'
