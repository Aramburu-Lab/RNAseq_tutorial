#!/bin/bash
# wo13-t: labdocs gates with their exit codes.
export PATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH
cd /Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ae58c42b269c90dca
labdocs lint; echo "lint exit $?"
labdocs render --check; echo "render --check exit $?"
