#!/bin/bash
# wo-z5: build every session root's pages from one tree (before: f287e07's; after: this
# worktree's src and tests copied), into scratch. Usage: build_fixture_pages.sh before|after
S=/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/z5
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a6b39a2692f007e36
which=$1
tree=$S/$which
if [ "$which" = after ]; then
  rm -rf "$tree"; mkdir -p "$tree"
  git -C "$WT" ls-files -z src tests conf scripts containers knowledge settings bin pyproject.toml pixi.toml ruff.toml README.md LICENSE .labtemplate.yml CONVENTIONS.md AGENTS.md CLAUDE.md TODO.md CHANGELOG.md logs outputs .gitignore | (cd "$WT" && xargs -0 tar -cf -) | tar -xf - -C "$tree"
  cp -R "$WT"/src/. "$tree"/src/   # uncommitted edits too
fi
cp "$S/test_z5_build_pages.py" "$tree/tests/"
rm -rf "$S/pages_$which" "$S/base_$which"
export Z5_OUT=$S/pages_$which
Z5_SRC=$tree/src Z5_CWD=$tree "$S/env_z5.sh" python -m pytest -p no:cacheprovider -q -n 6 \
  --basetemp="$S/base_$which" tests/test_z5_build_pages.py 2>&1 | tail -30
