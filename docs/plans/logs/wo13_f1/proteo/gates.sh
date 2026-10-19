#!/bin/bash
# wo13-f1 M7: every gate the brief names, through the worktree wrapper, each exit code printed.
PY=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a0cad45c60690dabd/docs/plans/logs/wo13_f1_proteo_py.sh
"$PY" ruff check . ; echo "ruff check: $?"
"$PY" ruff format --check . ; echo "ruff format: $?"
"$PY" mypy src | tail -1 ; echo "mypy: ${PIPESTATUS[0]}"
"$PY" interrogate src | tail -1 ; echo "interrogate: ${PIPESTATUS[0]}"
"$PY" labdocs lint | tail -2 ; echo "labdocs lint: ${PIPESTATUS[0]}"
"$PY" labdocs render --check | tail -2 ; echo "labdocs render --check: ${PIPESTATUS[0]}"
"$PY" python -m ctx.report.prose | tail -2 ; echo "prose: ${PIPESTATUS[0]}"
"$PY" python -m ctx.report.glossary_lint | tail -2 ; echo "glossary: ${PIPESTATUS[0]}"
"$PY" wc -l src/ctx/figures/kinds_m7.py src/ctx/figures/draw_m7.py src/ctx/proteo/figure_tables.py tests/test_figures_m7.py
"$PY" python -m pytest -q -n 4 -p no:cacheprovider tests/test_figures_m7.py tests/test_proteo_annotate.py tests/test_proteo_budget.py tests/test_proteo_folds.py tests/test_proteo_measure.py tests/test_proteo_page.py tests/test_proteo_prepare.py tests/test_proteo_study.py tests/test_proteo_tumour_normal.py tests/test_proteo_views_browser.py tests/test_figures.py tests/test_text.py | tail -3 ; echo "pytest: ${PIPESTATUS[0]}"
