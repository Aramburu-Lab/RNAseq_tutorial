#!/bin/bash
# The kit tests alone, in this worktree's code: K_pytest.sh [pytest args].
W="$(cd "$(dirname "$0")/../../../.." && pwd)"
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
export PATH="$BIN:$PATH" PYTHONPATH="$W/src"
cd "$W" || exit 1
python -m pytest -q tests/test_kit.py "$@"
