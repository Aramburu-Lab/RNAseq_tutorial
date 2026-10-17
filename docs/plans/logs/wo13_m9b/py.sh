#!/bin/bash
# Package M9b: python from main's environment on this worktree's src, run from the worktree.
. "$(dirname "$0")/env.sh"
cd "$WT" || exit 1
exec python "$@"
