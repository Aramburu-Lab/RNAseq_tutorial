#!/bin/bash
# Every before/after pair compared, text, Chromium and timing, one after another.
S=/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/z5
R=$S/real
"$S/cmp.sh" "$S/pages_before" "$S/pages_after" "$S/compare_fixture.log" --browser --timing 5 \
  --same 'scratchpad/z5/base_(before|after)/[^/]+' --same 'scratchpad/z5/(before|after)/'
for root in tier1 tcga_brca; do
  mkdir -p "$R/pages_${root}_before" "$R/pages_${root}_after"
  for w in before after; do
    rm -rf "$R/pages_${root}_$w"; mkdir -p "$R/pages_${root}_$w"
    for d in "$R/${root}_$w"/outputs/*/; do
      [ -L "${d%/}" ] && continue   # linked from the real root, not rebuilt here
      n=$(basename "$d"); mkdir -p "$R/pages_${root}_$w/$n"
      cp "$d"/*.html "$R/pages_${root}_$w/$n/" 2>/dev/null
    done
  done
  "$S/cmp.sh" "$R/pages_${root}_before" "$R/pages_${root}_after" "$R/compare_$root.log" --browser \
    --timing 5 --same "scratchpad/z5/real/${root}_(before|after)" --same 'scratchpad/z5/(before/src|after)'
done
echo ALLDONE
