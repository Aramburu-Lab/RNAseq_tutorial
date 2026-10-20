#!/bin/bash
# Package PK: pytest with main's environment, this worktree's src and Chromium; arguments pass.
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aa7a952e67debab5e
export PATH=$BIN:$PATH PYTHONPATH=$WT/src PYTHONDONTWRITEBYTECODE=1
export CTX_CHROMIUM=${CTX_CHROMIUM:-/Applications/Chromium.app/Contents/MacOS/Chromium}
cd "$WT" || exit 1
exec python -m pytest "$@"
