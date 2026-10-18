#!/usr/bin/env bash
# wo13-v: every gate but the full suite, on this worktree, with main's environment binaries.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ad2999f5c1b3365b1
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
export PATH=$BIN:$PATH PYTHONPATH=$WT/src PYTHONDONTWRITEBYTECODE=1
cd "$WT" || exit 1
echo "== ruff check"; ruff check . | tail -20
echo "== ruff format"; ruff format --check . | tail -5
echo "== mypy"; mypy src | tail -3
echo "== interrogate"; python -m interrogate -vv src/ctx | grep -E "MISSED|RESULT"
echo "== prose"; python -m ctx.report.prose | tail -5
echo "== glossary"; python -m ctx.report.glossary_lint | grep -v "^warning" | tail -5
echo "== labdocs lint"; labdocs lint | tail -10
echo "== labdocs render --check"; labdocs render --check | tail -10
