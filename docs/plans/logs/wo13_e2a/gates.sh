#!/bin/bash
# Package wo13-e2a: the static gates by hand, main's environment, this worktree's src.
#   bash gates.sh [pytest args...]   with args, also runs pytest on them
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a68d1b8f18deeebdd
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
export PATH=$BIN:/opt/homebrew/bin:/usr/bin:/bin PYTHONPATH=$WT/src PYTHONDONTWRITEBYTECODE=1
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
wc -l src/ctx/report/templates/*.j2 | sort -n | tail -3
if [ $# -gt 0 ]; then step python -m pytest -q -p no:cacheprovider "$@"; fi
