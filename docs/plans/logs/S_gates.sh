#!/bin/bash
# Package S: the static gates by hand, main's environment, this worktree's src.
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a8ab53e1e7f33bad4
export PATH=$BIN:$PATH PYTHONPATH=$WT/src
cd "$WT" || exit 1
step() { echo "== $*"; "$@"; echo "exit $?"; }
step ruff check src tests scripts docs/plans/logs
step ruff format --check src tests scripts docs/plans/logs
step mypy src
step interrogate -q src/ctx
step labdocs lint
step labdocs render --check
step python -m ctx.report.prose
step python -m ctx.report.glossary_lint
