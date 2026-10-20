#!/bin/bash
# wo13-z: every gate as CI runs it, against this worktree; exits non-zero on the first failure.
W=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a6e30099a3085f55d
export PATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH
export PYTHONPATH=$W/src PYTHONDONTWRITEBYTECODE=1
cd "$W" || exit 1
set -e
echo "== ruff check"; ruff check .
echo "== ruff format"; ruff format --check .
echo "== mypy"; mypy src
echo "== interrogate"; interrogate src
echo "== labdocs lint"; labdocs lint
echo "== labdocs render --check"; labdocs render --check
echo "== prose"; python -m ctx.report.prose
echo "== glossary"; python -m ctx.report.glossary_lint
# The prek hook's scope (.pre-commit-config.yaml check-file-length): scripts/ and bin/ only.
echo "== lengths over 600 among changed files under scripts/ and bin/"
over=$(git diff --name-only main -- 'scripts/*' 'bin/*' | while read -r p; do
  if [ -f "$p" ]; then n=$(wc -l < "$p"); [ "$n" -gt 600 ] && echo "$p $n"; fi
done || true)
[ -z "$over" ] || { echo "$over"; exit 1; }
echo "== all gates passed"
