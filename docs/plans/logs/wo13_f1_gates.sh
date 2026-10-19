#!/bin/bash
# wo13-f1: every gate on the package's files, as CI runs them, against this worktree.
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
W=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a00e5034a0fa89898
export PATH=$BIN:$PATH PYTHONPATH=$W/src PYTHONDONTWRITEBYTECODE=1
cd $W || exit 1
echo "== ruff check"; ruff check .; echo "exit $?"
echo "== ruff format"; ruff format --check .; echo "exit $?"
echo "== mypy"; mypy src; echo "exit $?"
echo "== interrogate"; interrogate src; echo "exit $?"
echo "== labdocs lint"; labdocs lint; echo "exit $?"
echo "== labdocs render --check"; labdocs render --check; echo "exit $?"
echo "== prose"; python -m ctx.report.prose; echo "exit $?"
echo "== glossary"; python -m ctx.report.glossary_lint; echo "exit $?"
