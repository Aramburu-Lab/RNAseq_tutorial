#!/bin/bash
# wo13-t: the package's own tests, against this worktree's src.
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
export PATH=$BIN:$PATH PYTHONPATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ae58c42b269c90dca/src PYTHONDONTWRITEBYTECODE=1
cd /Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ae58c42b269c90dca
pytest -q -p no:cacheprovider tests/test_gdc_door.py tests/test_kit_datasets.py tests/test_kit.py tests/test_qc_gdc_door.py tests/test_gdc.py "$@"
