#!/bin/bash
# wo13-t: the full suite once, as the brief asks, with the skip reasons.
export PATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH
export PYTHONPATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ae58c42b269c90dca/src
cd /Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ae58c42b269c90dca
pytest -q -n auto --dist loadscope -rs -p no:cacheprovider
