#!/bin/bash
# Package W, P.1: GS.6 (Open beside) and the picker in Microsoft Edge, through the same CDP
# driver: tests/test_pick_browser.py with CTX_CHROMIUM at Edge's binary.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aea85c67a862f2e82
SCRATCH=/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad
log=$WT/docs/plans/logs/wo13_w/gs6_edge.log
edge="/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge"
echo "== $(date -u +%Y-%m-%dT%H:%M:%SZ) $("$edge" --version)" > "$log"
CTX_CHROMIUM="$edge" bash "$WT/docs/plans/logs/wo13_w/pytest.sh" -v -p no:cacheprovider \
  --basetemp="$SCRATCH/w_gs6_edge_bt" tests/test_pick_browser.py tests/test_ora_browser.py \
  tests/test_genecard_browser.py >> "$log" 2>&1
echo "exit $?" >> "$log"
grep -E "PASSED|FAILED|SKIPPED|ERROR|passed|failed|^exit|^==" "$log"
