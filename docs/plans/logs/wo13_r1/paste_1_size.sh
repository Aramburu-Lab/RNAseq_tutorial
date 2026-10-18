export CTX_REPO=/Users/rytis/ctx_r1_clone RAW="$HOME/ctx_tcga_hnsc/data/raw"
git clone /Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a00f69bc2d903ff45 "$CTX_REPO"
pixi install --manifest-path "$CTX_REPO/pixi.toml" -e fetch --frozen
CTX_SITE=laptop DATA=0 bash "$CTX_REPO/scripts/bash/fetch_gdc.sbatch" TCGA-HNSC "$RAW/tcga_hnsc"
