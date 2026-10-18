#!/bin/bash
# wo13-r2: the full suite once, as the brief asks, with the skip reasons.
export PATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH
export PYTHONPATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a525b7995fa273222/src
cd /Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a525b7995fa273222 || exit 1
CTX_CHROMIUM=/Applications/Chromium.app/Contents/MacOS/Chromium pytest -q -n auto --dist loadscope -rs -p no:cacheprovider
