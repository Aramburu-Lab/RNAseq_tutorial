#!/bin/bash
# Package W: main's environment, this worktree's src (SRC overrides). bash py.sh <python args>
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aea85c67a862f2e82
export PATH=$BIN:$PATH PYTHONPATH=${SRC:-$WT/src} PYTHONDONTWRITEBYTECODE=1
cd "$WT" || exit 1
exec python "$@"
