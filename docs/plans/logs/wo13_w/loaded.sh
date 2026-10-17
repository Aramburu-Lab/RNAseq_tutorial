#!/bin/bash
# Package W, K.13: run pytest with every core saturated by `yes`, then kill the load.
#   bash loaded.sh <log> <pytest args>
# The load starts before pytest and is killed when pytest ends, whatever its exit.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aea85c67a862f2e82
log=$1
shift
cores=$(sysctl -n hw.ncpu)
pids=()
for _ in $(seq "$cores"); do
  yes > /dev/null &
  pids+=($!)
done
echo "== load: $cores yes processes, $(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$log"
start=$(date +%s)
bash "$WT/docs/plans/logs/wo13_w/pytest.sh" "$@" >> "$log" 2>&1
rc=$?
kill "${pids[@]}" 2> /dev/null
echo "== exit $rc after $(($(date +%s) - start)) s under load" >> "$log"
exit $rc
