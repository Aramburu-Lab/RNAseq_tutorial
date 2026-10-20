#!/bin/bash
# wo13-y: every gate as CI runs it, against this worktree; exits non-zero on the first failure.
W=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a93c44dfd6e2a237b
export PATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH
export PYTHONPATH=$W/src PYTHONDONTWRITEBYTECODE=1
cd "$W" || exit 1
set -e
echo "== ruff check"; ruff check .
echo "== ruff format"; ruff format --check .
echo "== mypy"; mypy src
echo "== interrogate"; interrogate src
echo "== labdocs lint"; labdocs lint
echo "== labdocs render --check"; labdocs render --check
echo "== prose"; python -m ctx.report.prose
echo "== glossary"; python -m ctx.report.glossary_lint
echo "== lengths over 600 among changed files"
over=$(git diff --name-only main -- '*.py' '*.j2' '*.mjs' '*.js' '*.sh' | while read -r p; do
  if [ -f "$p" ]; then n=$(wc -l < "$p"); [ "$n" -gt 600 ] && echo "$p $n"; fi
done || true)
[ -z "$over" ] || { echo "$over"; exit 1; }
echo "== all gates passed"
