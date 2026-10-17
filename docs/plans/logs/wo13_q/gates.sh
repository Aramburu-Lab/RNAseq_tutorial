#!/bin/bash
# Package Q: every gate on the committed files, by hand (no pixi in a worktree).
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a7158d7cd646891a9
P="$WT/docs/plans/logs/wo13_q/py.sh"
"$P" ruff check . ; echo "ruff check exit $?"
"$P" ruff format --check . ; echo "ruff format exit $?"
"$P" mypy src ; echo "mypy exit $?"
"$P" python -m interrogate -v src/ctx | tail -3 ; echo "interrogate exit ${PIPESTATUS[0]}"
"$P" labdocs lint ; echo "labdocs lint exit $?"
"$P" labdocs render --check ; echo "labdocs render --check exit $?"
"$P" python -m ctx.report.prose ; echo "prose lint exit $?"
"$P" python -m ctx.report.glossary_lint ; echo "glossary lint exit $?"
