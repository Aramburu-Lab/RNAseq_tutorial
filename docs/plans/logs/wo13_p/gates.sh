#!/bin/bash
# Package P: the static gates by hand, main's environment, this worktree's src.
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a4bc57162a0d35df8
export PATH=$BIN:$PATH PYTHONPATH=$WT/src PYTHONDONTWRITEBYTECODE=1
cd "$WT" || exit 1
step() { echo "== $*"; "$@"; echo "exit $?"; }
step ruff check src tests scripts docs/plans/logs/wo13_p
step ruff format --check src tests scripts docs/plans/logs/wo13_p
step mypy src
step interrogate -q src/ctx
step labdocs lint
step labdocs render --check
step python -m ctx.report.prose
step python -m ctx.report.glossary_lint
# The prek file-length hook's own rule, over scripts/ and bin/ as it runs.
w=$(sed -n "s/^[[:space:]]*warn:[[:space:]]*//p" .labtemplate.yml)
f=$(sed -n "s/^[[:space:]]*fail:[[:space:]]*//p" .labtemplate.yml)
echo "== file length (warn $w, fail $f)"
rc=0
for p in $(git ls-files 'scripts/*' 'bin/*' | grep -E '\.(py|sh|bash|R|r|nf|rs)$'); do
  n=$(awk 'END{print NR}' "$p")
  if [ "$n" -gt "$f" ]; then echo "$p: $n > $f"; rc=1; fi
done
echo "exit $rc"
