#!/bin/bash
# wo13-v2: one command with main's environment and this worktree's src.
export PATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH
export PYTHONPATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-af13a714b9eb1c983/src
export PYTHONDONTWRITEBYTECODE=1
export CTX_CHROMIUM=/Applications/Chromium.app/Contents/MacOS/Chromium
cd /Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-af13a714b9eb1c983 || exit 1
"$@"
