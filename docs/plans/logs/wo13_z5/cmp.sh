#!/bin/bash
# Usage: cmp.sh <before dir> <after dir> <log> [compare.py flags]
S=/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/z5
b=$1 a=$2 log=$3
shift 3
"$S/env_z5.sh" python docs/plans/logs/wo13_z5/compare.py "$b" "$a" "$@" > "$log" 2>&1
echo "exit $?"
tail -4 "$log"
grep -cE "^   ok" "$log"
grep -E "FAIL|MISSING" "$log" | head -20
