#!/bin/bash
# wo13-f1 M9: the commit gates on the whole tree, one exit code each.
PY=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aafb5beef9062b91f/docs/plans/logs/wo13_f1_genecard_py.sh
run() { "$PY" "$@" > /dev/null 2>&1; echo "$*: exit $?"; }
run ruff check .
run ruff format --check .
run mypy src
run interrogate src
run labdocs lint
run labdocs render --check
run python -m ctx.report.prose
run python -m ctx.report.glossary_lint
"$PY" wc -l src/ctx/figures/kinds_m9.py src/ctx/figures/draw_m9.py src/ctx/genecard/figure_tables.py \
  src/ctx/genecard/page.py tests/test_figures_m9.py docs/plans/logs/wo13_f1/genecard/mutations.py
