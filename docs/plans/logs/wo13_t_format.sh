#!/bin/bash
# wo13-t: format and fix the package's files.
export PATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH
cd /Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ae58c42b269c90dca
ruff format scripts/python/preview_genecard.py scripts/python/build_kit.py tests/test_gdc_door.py tests/test_kit_datasets.py src/ctx/gather/gdc.py src/ctx/quant/run.py src/ctx/process/run.py src/ctx/qc/run.py
ruff check --fix tests/test_kit_datasets.py tests/test_gdc_door.py src/ctx/process/run.py src/ctx/quant/run.py scripts/python/preview_genecard.py scripts/python/preview_proteo.py
