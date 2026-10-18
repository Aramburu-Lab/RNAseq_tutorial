#!/bin/bash
# wo13-u R2.2: the guide's Part 1c lines, sourced by step. The clone is this worktree and
# `pixi run -e laptop` is main's environment with this worktree's src (never pixi here); the
# guide's `git clone` and `pixi install` are the only lines not run.
# Usage: source r22_env.sh <GSE> <slug>
export CTX_REPO=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aa1d340947c025c9d
export GSE=$1 SLUG=$2
export RAW="$HOME/ctx_$SLUG/data/raw/$SLUG"
PIXI_LAPTOP="$CTX_REPO/docs/plans/logs/wo13_u_py.sh"
export PIXI_LAPTOP
