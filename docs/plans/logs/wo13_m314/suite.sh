#!/bin/sh
# M3.14: one probed suite run. Usage: suite.sh <checkout> <run folder> <pytest args...>
# Runs in a clean clone, so the session guard's git status is the clone's; writes time.txt,
# out.txt, the load average every 30 s, bt/ (the basetemp) and probe/ into the run folder.
# Then: python m314_table.py <run folder> <label> "<command>" <out.md>
W=$1
O=$2
shift 2
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=${BIN:-/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin}
rm -rf "$O"
mkdir -p "$O/probe"
export M314_OUT="$O/probe"
export PATH="$BIN:$PATH"
export PYTHONPATH="$W/src:$HERE"
export CTX_CHROMIUM=${CTX_CHROMIUM:-/opt/homebrew/bin/chromium}
uptime > "$O/load_before.txt"
( while true; do uptime >> "$O/load_during.txt"; sleep 30; done ) &
LOADS=$!
cd "$W" || exit 1
/usr/bin/time -l "$BIN/python" -m pytest -p m314_probe --basetemp="$O/bt" "$@" > "$O/out.txt" 2> "$O/time.txt"
echo "exit $?" >> "$O/out.txt"
kill $LOADS
uptime > "$O/load_after.txt"
