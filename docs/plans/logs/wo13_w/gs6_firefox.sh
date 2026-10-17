#!/bin/bash
# Package W, P.1: GS.6 in Firefox (gs6_firefox.mjs) on the Step-1 scratch root's M5 and M9.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aea85c67a862f2e82
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
log=$WT/docs/plans/logs/wo13_w/gs6_firefox.log
m5=$(ls -t "$WT"/outputs/scratch_roots/tier1_m9/outputs/05_de/*_overview_v1.html | head -1)
echo "== $(date -u +%Y-%m-%dT%H:%M:%SZ) $(/Applications/Firefox.app/Contents/MacOS/firefox --version)" > "$log"
GS6_PROFILE_DIR=$WT/outputs/scratch_roots/ff "$BIN/node" "$WT/docs/plans/logs/wo13_w/gs6_firefox.mjs" /Applications/Firefox.app/Contents/MacOS/firefox \
  "file://$m5" >> "$log" 2>&1
echo "exit $?" >> "$log"
cat "$log"
