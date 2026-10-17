#!/bin/bash
# wo13-t: every gate on the package's files, as CI runs them, against this worktree.
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
export PATH=$BIN:$PATH PYTHONPATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ae58c42b269c90dca/src PYTHONDONTWRITEBYTECODE=1
cd /Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ae58c42b269c90dca
echo "== ruff check"; ruff check .
echo "== ruff format"; ruff format --check .
echo "== mypy"; mypy src
echo "== interrogate"; interrogate src
echo "== labdocs lint"; labdocs lint
echo "== labdocs render --check"; labdocs render --check
echo "== prose"; python -m ctx.report.prose
echo "== glossary"; python -m ctx.report.glossary_lint
