#!/bin/bash
# Package wo13-de: the proofs on the commit, in turn: final.sh, then every planted mutation.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ab82f97e3ba379d33
here=$WT/docs/plans/logs/wo13_de
bash "$here/final.sh" "$1"
rm -f "$here/mutations.log"
bash "$here/env.sh" python docs/plans/logs/wo13_de/mutate.py docs/plans/logs/wo13_de/mutations.json \
  docs/plans/logs/wo13_de/mutations.log
