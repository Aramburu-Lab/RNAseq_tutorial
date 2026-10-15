#!/bin/bash
# One planted mutation per K.9 and K.10 rule; see K_mutate.py.
W="$(cd "$(dirname "$0")/../../../.." && pwd)"
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
export PATH="$BIN:$PATH" PYTHONPATH="$W/src"
python "$W/docs/plans/logs/wo13_k/K_mutate.py"
