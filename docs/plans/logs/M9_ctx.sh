#!/bin/bash
# Package M9 (wo13-m9): one timed ctx command on a root, with this worktree's src; its log here.
#   bash M9_ctx.sh <root> <log name> <ctx arguments after --root>
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a209f728a6cb9c172
export PATH=$BIN:$PATH PYTHONPATH=$WT/src
root=$1 log=$WT/docs/plans/logs/$2.log
shift 2
echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) ctx --root $root $*" > "$log"
/usr/bin/time -l python -m ctx --root "$root" "$@" >> "$log" 2>&1
echo "exit $?" >> "$log"
tail -25 "$log"
