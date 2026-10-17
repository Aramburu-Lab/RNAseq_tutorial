#!/bin/bash
# Package wo13-e2b: run one command with main's environment and this worktree's src.
#   bash env.sh <command ...>
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a6353f37f240ddb6a
export PATH=$BIN:$PATH PYTHONPATH=${SRC:-$WT/src} PYTHONDONTWRITEBYTECODE=1
cd "$WT" || exit 1
exec "$@"
