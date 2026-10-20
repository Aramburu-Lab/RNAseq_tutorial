#!/bin/bash
# Package PK: main's environment, this worktree's src. bash py.sh <python args>
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aa7a952e67debab5e
export PATH=$BIN:$PATH PYTHONPATH=$WT/src PYTHONDONTWRITEBYTECODE=1
# Without Chromium the browser tests skip, and a skip reads as a surviving mutation.
export CTX_CHROMIUM=${CTX_CHROMIUM:-/Applications/Chromium.app/Contents/MacOS/Chromium}
cd "$WT" || exit 1
exec python "$@"
