#!/bin/bash
# Mutation 2: the stored block's rewind of the bytes the table decode read ahead taken out.
S=/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/z5
echo "== mutation 2: stored() without the read-ahead rewind (mut tree)"
Z5_SRC=$S/mut/src Z5_CWD=$S/mut "$S/env_z5.sh" python -m pytest -p no:cacheprovider -q \
  tests/test_packed_figures.py -k unpack_inflates 2>&1 | grep -E "passed|failed|^E  " | head -6
echo "== the same test on this worktree"
"$S/env_z5.sh" python -m pytest -p no:cacheprovider -q tests/test_packed_figures.py \
  -k unpack_inflates 2>&1 | tail -2
