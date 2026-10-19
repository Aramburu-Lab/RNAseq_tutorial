#!/bin/bash
# wo13-f1 M9: one command with main's environment and this worktree's src.
export PATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH
export PYTHONPATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aafb5beef9062b91f/src
export PYTHONDONTWRITEBYTECODE=1
export CTX_CHROMIUM=/Applications/Chromium.app/Contents/MacOS/Chromium
cd /Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aafb5beef9062b91f || exit 1
"$@"
