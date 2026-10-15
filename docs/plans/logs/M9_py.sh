#!/bin/bash
# Package M9 (wo13-m9): run a Python script with main's environment and this worktree's src.
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
export PATH=$BIN:$PATH
export PYTHONPATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a209f728a6cb9c172/src
exec python "$@"
