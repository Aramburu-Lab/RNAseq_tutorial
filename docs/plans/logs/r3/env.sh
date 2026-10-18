#!/usr/bin/env bash
# R3: run a command on main's default environment with this worktree's src first (lessons.md:
# never run pixi in a worktree; PYTHONPATH absolute so test subprocesses follow it).
set -euo pipefail
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ac5a02af0b00495c4
export PATH="/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH"
export PYTHONPATH="$WT/src"
export PYTHONDONTWRITEBYTECODE=1
cd "$WT"
exec "$@"
