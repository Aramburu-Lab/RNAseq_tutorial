#!/usr/bin/env bash
# R3: the full suite once, as the brief asks, with its skip reasons.
set -uo pipefail
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ac5a02af0b00495c4
export PATH="/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH"
export PYTHONPATH="$WT/src"
cd "$WT"
pytest -q -n auto --dist loadscope -rs -p no:cacheprovider
echo "suite exit $?"
