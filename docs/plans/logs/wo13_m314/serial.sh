#!/bin/sh
# M3.14: each of the ten slowest files run serially (no xdist), the ten at once, each under
# /usr/bin/time -l. Usage: serial.sh <checkout> <out folder>. The ten are the files of most
# summed setup and call time in the run before the change (before.md).
W=$1
O=$2
BIN=${BIN:-/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin}
rm -rf "$O"
mkdir -p "$O"
export PATH="$BIN:$PATH"
export PYTHONPATH="$W/src"
export CTX_CHROMIUM=${CTX_CHROMIUM:-/opt/homebrew/bin/chromium}
cd "$W" || exit 1
uptime > "$O/load_before.txt"
( while true; do uptime >> "$O/load_during.txt"; sleep 30; done ) &
LOADS=$!
for f in test_figures_m2 test_cli test_qc_run test_figures_m7 test_outputs_contract \
         test_figures test_process_run test_pipeline test_figures_m8 test_figures_m9; do
  /usr/bin/time -l "$BIN/python" -m pytest -q -p no:xdist -p no:cacheprovider \
    --basetemp="$O/bt_$f" "tests/$f.py" > "$O/$f.out" 2> "$O/$f.time" &
done
sleep 5
while pgrep -f "basetemp=$O/bt_" > /dev/null; do sleep 5; done
kill $LOADS
uptime > "$O/load_after.txt"
