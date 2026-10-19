#!/bin/bash
# wo13-s2: the full suite once, as the brief asks, with the skip reasons.
export PATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH
export PYTHONPATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aa11ca42a7ace51a9/src
cd /Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aa11ca42a7ace51a9 || exit 1
CTX_CHROMIUM=/Applications/Chromium.app/Contents/MacOS/Chromium pytest -q -n auto --dist loadscope -rs -p no:cacheprovider
