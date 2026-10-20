#!/bin/bash
# Package PK: the hand-off tests in four runs at once, each building its own root, so each runs
# under the others' load, as the full suite's xdist workers ran it when the reveal first failed.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aa7a952e67debab5e
SCRATCH=/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad
label=${1:-load}
for i in 1 2 3 4; do
  rm -rf "$SCRATCH/pk_${label}_$i"
  bash "$WT/docs/plans/logs/wo13_pk/pytest.sh" -q -p no:cacheprovider \
    --basetemp="$SCRATCH/pk_${label}_$i" tests/test_handoff_browser.py \
    > "$WT/docs/plans/logs/wo13_pk/${label}_$i.log" 2>&1 &
done
wait
uptime
tail -qn 1 "$WT"/docs/plans/logs/wo13_pk/"${label}"_*.log
grep -h "AssertionError" "$WT"/docs/plans/logs/wo13_pk/"${label}"_*.log
