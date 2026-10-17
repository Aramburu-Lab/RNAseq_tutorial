#!/bin/bash
# Package W: pytest with main's environment and this worktree's src; arguments pass.
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aea85c67a862f2e82
export PATH=$BIN:$PATH PYTHONPATH=$WT/src PYTHONDONTWRITEBYTECODE=1
cd "$WT" || exit 1
exec python -m pytest "$@"
