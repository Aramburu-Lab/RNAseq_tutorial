#!/usr/bin/env bash
# Package F2's runner: main's pixi default environment, this worktree's src, then the command.
#   bash docs/plans/logs/2026-10-19_f2_env.sh pytest -q tests/test_figure_panel.py
W=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ad214d746812b3316
export PATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH
export PYTHONPATH=$W/src
export CTX_CHROMIUM=/Applications/Chromium.app/Contents/MacOS/Chromium
cd "$W" || exit 1
exec "$@"
