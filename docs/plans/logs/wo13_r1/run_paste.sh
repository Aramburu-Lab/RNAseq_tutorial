#!/usr/bin/env bash
# wo13-r1 proof: run pasted guide blocks in one shell, as a student's terminal would, with the
# pixi shim first on PATH and no project environment active. Usage: run_paste.sh <paste>...
L=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a00f69bc2d903ff45/docs/plans/logs/wo13_r1
export PATH="$L/shim:/usr/bin:/bin:/usr/sbin:/sbin" PIXI_SHIM_LOG="$L/pixi_shim.log"
cd "$HOME" || exit 1
for paste in "$@"; do
    echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) $paste"
    # shellcheck disable=SC1090
    source "$L/$paste"
done
