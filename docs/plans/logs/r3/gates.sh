#!/usr/bin/env bash
# R3's gates on every file, as preflight runs them, apart from the full suite (suite.sh).
set -uo pipefail
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ac5a02af0b00495c4
export PATH="/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH"
export PYTHONPATH="$WT/src"
export PYTHONDONTWRITEBYTECODE=1
cd "$WT"
status=0
step() { echo "== $*"; "$@" || { echo "FAILED: $*"; status=1; }; }
step ruff check .
step ruff format --check .
step mypy src
step python -m interrogate -c pyproject.toml src/ctx --fail-under 100 -q
step labdocs lint
step labdocs render --check
step python -m ctx.report.prose
step python -m ctx.report.glossary_lint
exit $status
