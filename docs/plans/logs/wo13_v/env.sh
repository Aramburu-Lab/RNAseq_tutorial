#!/usr/bin/env bash
# wo13-v: main's environment binaries with this worktree's src first (no pixi in a worktree).
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ad2999f5c1b3365b1
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
export PATH=$BIN:$PATH PYTHONPATH=$WT/src PYTHONDONTWRITEBYTECODE=1
export CTX_CHROMIUM=/Applications/Chromium.app/Contents/MacOS/Chromium
cd "$WT" || exit 1
"$@"
