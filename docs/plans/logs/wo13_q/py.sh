#!/bin/bash
# Package Q: run a command with main's environment and this worktree's src.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a7158d7cd646891a9
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
export PATH=$BIN:$PATH
export PYTHONPATH=$WT/src
cd "$WT" || exit 1
exec "$@"
