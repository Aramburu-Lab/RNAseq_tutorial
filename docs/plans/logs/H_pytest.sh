#!/bin/bash
# Package H: pytest with main's environment and this worktree's src; arguments pass through.
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a01b45b7d9228a7fd
export PATH=$BIN:$PATH PYTHONPATH=$WT/src
cd "$WT" || exit 1
exec python -m pytest "$@"
