#!/bin/bash
# Package K's gates, in this worktree's code, on the files K touches; each prints its own verdict.
W="$(cd "$(dirname "$0")/../../../.." && pwd)"
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
export PATH="$BIN:$PATH" PYTHONPATH="$W/src"
cd "$W" || exit 1
FILES="scripts/python/build_kit.py scripts/python/check_kit.py tests/test_kit.py docs/plans/logs/wo13_k/K_mutate.py docs/plans/logs/wo13_k/K_kidney_check.py"
echo "== ruff check"; ruff check $FILES
echo "== ruff format --check"; ruff format --check $FILES
echo "== mypy src"; mypy src | tail -1
echo "== interrogate"; interrogate src | tail -1
echo "== labdocs lint"; labdocs lint; echo "exit $?"
echo "== labdocs render --check"; labdocs render --check; echo "exit $?"
echo "== prose lint"; python -m ctx.report.prose
echo "== glossary lint"; python -m ctx.report.glossary_lint | tail -3
