#!/bin/bash
# Package P: main's environment, this worktree's src. bash py.sh <python args>
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a4bc57162a0d35df8
export PATH=$BIN:$PATH PYTHONPATH=${SRC:-$WT/src} PYTHONDONTWRITEBYTECODE=1
cd "$WT" || exit 1
exec python "$@"
