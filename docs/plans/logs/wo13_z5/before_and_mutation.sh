#!/bin/bash
# The new tests on f287e07's src (must fail), and on this src with one planted mutation.
S=/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/z5
Z5_SRC=$S/fb/src Z5_CWD=$S/fb "$S/env_z5.sh" python -m pytest -p no:cacheprovider -q \
  tests/test_packed_figures_z5.py 2>&1 | tail -6 > "$S/tests_before.log"
cat "$S/tests_before.log"
echo "== mutation: figure_html returns the figure's script unpacked"
Z5_SRC=$S/mut/src Z5_CWD=$S/mut "$S/env_z5.sh" python -m pytest -p no:cacheprovider -q \
  tests/test_packed_figures.py -k "not chromium and not unpack_inflates" 2>&1 \
  | grep -E "passed|failed|FAILED|^E " | head -14 > "$S/mutation.log"
cat "$S/mutation.log"
