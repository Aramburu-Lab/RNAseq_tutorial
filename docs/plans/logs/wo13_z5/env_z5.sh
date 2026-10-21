#!/bin/bash
# wo-z5 environment: the default pixi env's tools, this worktree's src; never pixi itself.
export BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
export WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a6b39a2692f007e36
export PYTHONPATH=${Z5_SRC:-$WT/src} PYTHONDONTWRITEBYTECODE=1 CTX_CHROMIUM=/opt/homebrew/bin/chromium PATH=$BIN:$PATH
cd ${Z5_CWD:-$WT}
exec "$@"
