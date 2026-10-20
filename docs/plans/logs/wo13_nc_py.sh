#!/bin/bash
# wo13-nc: one command with main's environment and this worktree's src.
export PATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH
export PYTHONPATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-af2e7af7de3e7d833/src
export PYTHONDONTWRITEBYTECODE=1
export CTX_CHROMIUM=/Applications/Chromium.app/Contents/MacOS/Chromium
cd /Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-af2e7af7de3e7d833 || exit 1
"$@"
