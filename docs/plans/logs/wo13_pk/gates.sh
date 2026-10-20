#!/bin/bash
# Package PK: the static gates by hand, main's environment, this worktree's src.
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aa7a952e67debab5e
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
step python -m pytest -q -p no:cacheprovider tests/test_templates.py
wc -l src/ctx/report/templates/*.j2 | sort -n | tail -3
