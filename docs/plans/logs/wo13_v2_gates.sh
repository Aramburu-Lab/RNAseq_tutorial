#!/bin/bash
# wo13-v2: every gate on the package's files, as CI runs them, against this worktree.
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
export PATH=$BIN:$PATH PYTHONPATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-af13a714b9eb1c983/src PYTHONDONTWRITEBYTECODE=1
cd /Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-af13a714b9eb1c983 || exit 1
echo "== ruff check"; ruff check .; echo "exit $?"
echo "== ruff format"; ruff format --check .; echo "exit $?"
echo "== mypy"; mypy src; echo "exit $?"
echo "== interrogate"; interrogate src; echo "exit $?"
echo "== labdocs lint"; labdocs lint; echo "exit $?"
echo "== labdocs render --check"; labdocs render --check; echo "exit $?"
echo "== prose"; python -m ctx.report.prose; echo "exit $?"
echo "== glossary"; python -m ctx.report.glossary_lint; echo "exit $?"
echo "== lengths over 600 among changed files"
git diff --name-only main -- '*.py' '*.j2' '*.mjs' '*.js' '*.sh' | while read -r p; do [ -f "$p" ] && n=$(wc -l < "$p") && [ "$n" -gt 600 ] && echo "$p $n"; done; echo "done"
