#!/bin/bash
# Package M9b: the full suite once, as main runs it, into a log beside this script.
. "$(dirname "$0")/env.sh"
cd "$WT" || exit 1
log=/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/m9b_suite.log
echo "== $(date -u +%Y-%m-%dT%H:%M:%SZ) $(git rev-parse --short HEAD) $(uptime | sed 's/.*load/load/')" > "$log"
python -m pytest -q -n auto --dist loadscope -rs -p no:cacheprovider >> "$log" 2>&1
echo "exit $?" >> "$log"
