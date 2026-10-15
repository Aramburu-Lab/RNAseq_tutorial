#!/bin/bash
# The full suite at package K's head, as CI runs it, with wall time and peak RAM.
W="$(cd "$(dirname "$0")/../../../.." && pwd)"
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
export PATH="$BIN:$PATH" PYTHONPATH="$W/src"
cd "$W" || exit 1
git log --oneline -1
/usr/bin/time -l python -m pytest -q -n auto --dist loadscope -rs "$@" 2>&1 \
    | grep -v -E "^ +[0-9]+  (page|block|messages|signals|voluntary|involuntary|instructions|cycles|swaps|average)" \
    | grep -E "FAILED|SKIPPED|ERROR|passed|failed|real|maximum resident|peak memory|^E |Error"
