#!/bin/bash
# wo13-r2: every gate on the package's files, as CI runs them, against this worktree.
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
export PATH=$BIN:$PATH PYTHONPATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a525b7995fa273222/src PYTHONDONTWRITEBYTECODE=1
cd /Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a525b7995fa273222 || exit 1
echo "== ruff check"; ruff check .; echo "exit $?"
echo "== ruff format"; ruff format --check .; echo "exit $?"
echo "== mypy"; mypy src; echo "exit $?"
echo "== interrogate"; interrogate src; echo "exit $?"
echo "== labdocs lint"; labdocs lint; echo "exit $?"
echo "== labdocs render --check"; labdocs render --check; echo "exit $?"
echo "== prose"; python -m ctx.report.prose; echo "exit $?"
echo "== glossary"; python -m ctx.report.glossary_lint; echo "exit $?"
