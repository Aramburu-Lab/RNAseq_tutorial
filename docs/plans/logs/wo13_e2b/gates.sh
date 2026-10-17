#!/bin/bash
# Package wo13-e2b: the static gates by hand, main's environment, this worktree's src.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a6353f37f240ddb6a
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
step python -m pytest -q tests/test_templates.py
