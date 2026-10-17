#!/bin/bash
# Package W, E.1a: m9_browser.py on a scratch root's newest M9 page and its M5 page.
#   bash m9_check.sh <root folder name, as tier1_m9>
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aea85c67a862f2e82
SCRATCH=/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad
out=$WT/outputs/scratch_roots/$1/outputs
m9=$(ls -t "$out"/11_genecard/*_overview_v1.html | head -1)
m5=$(ls -t "$out"/05_de/*_overview_v1.html | head -1)
log=$WT/docs/plans/logs/wo13_w/m9_browser_$1.log
echo "== $(date -u +%Y-%m-%dT%H:%M:%SZ) ${m9##*/} ${m5##*/}" > "$log"
bash "$WT/docs/plans/logs/wo13_w/py.sh" "$WT/docs/plans/logs/wo13_w/m9_browser.py" "$m9" "$m5" \
  "$SCRATCH/w_m9_$1_raw.json" >> "$log" 2>&1
echo "exit $?" >> "$log"
tail -80 "$log"
