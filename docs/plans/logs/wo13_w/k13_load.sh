#!/bin/bash
# Package W, K.13: run a command with n `yes` processes busy, then kill them.
#   bash k13_load.sh <n> <command ...>
n=$1
shift
pids=()
for _ in $(seq "$n"); do
  yes > /dev/null &
  pids+=($!)
done
"$@"
rc=$?
kill "${pids[@]}" 2> /dev/null
wait 2> /dev/null
exit $rc
