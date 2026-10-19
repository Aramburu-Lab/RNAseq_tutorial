#!/bin/bash
# wo13-s2: one command with main's environment and this worktree's src.
export PATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH
export PYTHONPATH=${S2_SRC:-/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aa11ca42a7ace51a9/src}
export PYTHONDONTWRITEBYTECODE=1
export CTX_CHROMIUM=/Applications/Chromium.app/Contents/MacOS/Chromium
cd /Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aa11ca42a7ace51a9 || exit 1
"$@"
