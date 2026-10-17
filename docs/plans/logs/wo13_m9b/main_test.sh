#!/bin/bash
# Package M9b: a test of main d67de7b, extracted by git archive into a scratch folder, on its
# own src, to tell a failure of this package from one main has.
#   bash main_test.sh <pytest args>
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
MAIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-af0224262630483d8/outputs/scratch_roots/m9b_main
export PATH=$BIN:$PATH PYTHONPATH=$MAIN/src PYTHONDONTWRITEBYTECODE=1
cd "$MAIN" || exit 1
exec python -m pytest -p no:cacheprovider "$@"
