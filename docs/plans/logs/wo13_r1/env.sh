#!/usr/bin/env bash
# wo13-r1: main's environment binaries with this worktree's src first (no pixi in a worktree).
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a00f69bc2d903ff45
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
export PATH=$BIN:$PATH PYTHONPATH=$WT/src PYTHONDONTWRITEBYTECODE=1
cd "$WT" || exit 1
"$@"
