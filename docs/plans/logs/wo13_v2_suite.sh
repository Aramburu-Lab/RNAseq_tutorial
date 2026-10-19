#!/bin/bash
# wo13-v2: the full suite once, as the brief asks, with the skip reasons.
export PATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH
export PYTHONPATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-af13a714b9eb1c983/src
cd /Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-af13a714b9eb1c983 || exit 1
CTX_CHROMIUM=/Applications/Chromium.app/Contents/MacOS/Chromium pytest -q -n auto --dist loadscope -rs -p no:cacheprovider
