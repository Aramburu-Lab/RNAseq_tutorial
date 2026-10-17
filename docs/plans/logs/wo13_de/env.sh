#!/bin/bash
# Package wo13-de: run one command with main's environment and this worktree's src.
#   [SRC=<another src>] bash env.sh <command ...>
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ab82f97e3ba379d33
export PATH=$BIN:$PATH PYTHONPATH=${SRC:-$WT/src} PYTHONDONTWRITEBYTECODE=1
cd "$WT" || exit 1
exec "$@"
