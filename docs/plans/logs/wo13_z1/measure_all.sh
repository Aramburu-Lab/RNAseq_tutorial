#!/bin/bash
# wo13-z1: measure.py on each cohort's scratch page, before and after.
HERE=$(dirname "$0")
ROOTS=/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/roots
for c in brca kirc luad lusc hnsc; do
  for s in before after; do
    python "$HERE/measure.py" "$(ls "$ROOTS/${s}_$c"/outputs/12_clinical/*_overview_v1.html)" \
      > "$HERE/measure_${c}_$s.log" 2>&1 || echo "failed: $c $s"
  done
done
