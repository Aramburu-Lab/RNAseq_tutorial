#!/bin/bash
# Package M9, H.2 (approved by Rytis 2026-10-15): the lung cards online, each run timed, the
# response cache's manifest and size read before and after each.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a209f728a6cb9c172
ROOT=$HOME/ctx_lung
LOGS=$WT/docs/plans/logs
cache() { python3 -c "import json,sys; print(len(json.load(open(sys.argv[1]))))" "$ROOT/resources/cache/manifest.json"; du -sk "$ROOT/resources/cache" | cut -f1; }
echo "before: $(cache | paste -sd' ' -) (entries, KiB)" > "$LOGS/M9_h2_cache.log"
bash "$LOGS/M9_ctx.sh" "$ROOT" M9_h2_luad_online genecard --dataset tcga_luad \
  --dataset-file "$WT/conf/datasets/tcga_luad.toml" --genes NKX2-1,SFTPC,EGFR,MKI67 --online
echo "after luad: $(cache | paste -sd' ' -) (entries, KiB)" >> "$LOGS/M9_h2_cache.log"
bash "$LOGS/M9_ctx.sh" "$ROOT" M9_h2_lusc_online genecard --dataset tcga_lusc \
  --dataset-file "$WT/conf/datasets/tcga_lusc.toml" --genes SOX2,TP63,KRT5,MKI67 --online
echo "after lusc: $(cache | paste -sd' ' -) (entries, KiB)" >> "$LOGS/M9_h2_cache.log"
cat "$LOGS/M9_h2_cache.log"
