#!/bin/bash
# V2: the one M6 mutation the fixture cannot catch (its contrasts test the same sets in the same
# order), planted, built into Step-1's M6 page and walked by m6_real.py, then restored.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-af13a714b9eb1c983
cd "$WT" || exit 1
file=src/ctx/report/templates/_page_running.js.j2
cp "$file" /tmp/v2_page_running.bak
old='      if (c !== shownContrast || (e \&\& e.target === picker)) relist(c);'
grep -c "if (c !== shownContrast || (e && e.target === picker)) relist(c);" "$file"
sed -i '' "/${old}/d" "$file"
grep -c "relist(c);" "$file"
DE_FROM=a docs/plans/logs/wo13_v2/build.sh relist tier1 function
cp /tmp/v2_page_running.bak "$file"
git -C "$WT" diff --stat -- "$file"
docs/plans/logs/wo13_v2_py.sh python docs/plans/logs/wo13_v2/m6_real.py \
  outputs/scratch_roots/v2_tier1_relist tier1 outputs/scratch_roots/v2_downloads/relist 2>&1 | tail -9
