#!/bin/bash
# Package H: one timed ctx command on a root, with this worktree's src; its log in logs/.
#   bash H_ctx.sh <root> <log name> <ctx arguments after --root>
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a01b45b7d9228a7fd
export PATH=$BIN:$PATH PYTHONPATH=$WT/src
root=$1 log=$WT/docs/plans/logs/$2.log
shift 2
echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) ctx --root $root $*" > "$log"
/usr/bin/time -l python -m ctx --root "$root" "$@" >> "$log" 2>&1
echo "exit $?" >> "$log"
tail -25 "$log"
