#!/bin/bash
# Package M9b: the static gates by hand, main's environment, this worktree's src.
. "$(dirname "$0")/env.sh"
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
