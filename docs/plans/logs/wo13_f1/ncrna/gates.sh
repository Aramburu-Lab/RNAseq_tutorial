#!/bin/bash
# F.1 M8: every gate the brief names, each exit code into gates.log.
W=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a84b7900e11d976d1/docs/plans/logs/wo13_f1_ncrna_py.sh
L=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a84b7900e11d976d1/docs/plans/logs/wo13_f1/ncrna/gates.log
: > "$L"
for c in "ruff check ." "ruff format --check ." "mypy src" "interrogate src" "labdocs lint" \
         "labdocs render --check" "python -m ctx.report.prose" "python -m ctx.report.glossary_lint"; do
  $W $c > /dev/null 2>&1
  echo "$c: exit $?" >> "$L"
done
cat "$L"
