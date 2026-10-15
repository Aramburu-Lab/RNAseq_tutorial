#!/bin/bash
# Package M9 (wo13-m9): the static gates by hand, main's environment, this worktree's src.
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a209f728a6cb9c172
export PATH=$BIN:$PATH PYTHONPATH=$WT/src
cd "$WT" || exit 1
step() { echo "== $*"; "$@"; echo "exit $?"; }
step ruff check src tests scripts docs
step ruff format --check src tests scripts docs
step mypy src
step interrogate -q src/ctx
step labdocs lint
step labdocs render --check
step python -m ctx.report.prose
step python -m ctx.report.glossary_lint
