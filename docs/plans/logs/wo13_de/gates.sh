#!/bin/bash
# Package wo13-de: the static gates by hand, main's environment, this worktree's src.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ab82f97e3ba379d33
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
export PATH=$BIN:$PATH PYTHONPATH=$WT/src PYTHONDONTWRITEBYTECODE=1
cd "$WT" || exit 1
step() { echo "== $*"; "$@"; echo "exit $?"; }
step ruff check .
step ruff format --check .
step mypy src
step interrogate -q src/ctx
step labdocs lint
step labdocs render --check
step python -m ctx.report.prose
step python -m ctx.report.glossary_lint
