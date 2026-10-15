#!/bin/bash
# Package M9 (wo13-m9): ruff check, ruff format --check and mypy on the files named.
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a209f728a6cb9c172
export PATH=$BIN:$PATH PYTHONPATH=$WT/src
cd "$WT" || exit 1
ruff check "$@"
ruff format --check "$@"
mypy src/ctx/genecard
interrogate -q src/ctx/genecard
