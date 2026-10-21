#!/bin/bash
# wo13-z1: measure.py on each cohort's page before and after, and compare.py with Chromium.
HERE=$(dirname "$0")
ROOTS=/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/roots
for c in brca kirc luad lusc hnsc; do
  b=$(ls "$ROOTS/before_$c"/outputs/12_clinical/*_overview_v1.html)
  a=$(ls "$ROOTS/after_$c"/outputs/12_clinical/*_overview_v1.html)
  python "$HERE/measure.py" "$b" > "$HERE/measure_${c}_before.log" 2>&1
  python "$HERE/measure.py" "$a" > "$HERE/measure_${c}_after.log" 2>&1
  [ "$c" = brca ] && continue
  python "$HERE/compare.py" "$b" "$a" --browser > "$HERE/compare_$c.log" 2>&1
  echo "$c: $(tail -1 "$HERE/compare_$c.log")"
done
