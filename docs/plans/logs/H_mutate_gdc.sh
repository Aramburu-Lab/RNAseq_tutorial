#!/bin/bash
# Package H: plant three mutations in the second-entity check, run its test, restore each.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a01b45b7d9228a7fd
F=$WT/src/ctx/gather/gdc.py
export PYTHONDONTWRITEBYTECODE=1
cp "$F" "$F.orig"
run() {
  python3 - "$F" "$1" "$2" <<'PY'
import sys
path, old, new = sys.argv[1:]
text = open(path).read()
assert text.count(old) == 1, old
open(path, "w").write(text.replace(old, new))
PY
  if bash "$WT/docs/plans/logs/H_pytest.sh" -q -p no:cacheprovider \
      "tests/test_gdc.py::test_an_aliquot_listed_under_two_samples_of_its_case_is_one_file" \
      "tests/test_gdc.py::test_a_file_of_two_cases_or_another_project_is_refused" >/dev/null 2>&1; then
    echo "SURVIVED: $3"
  else
    echo "caught: $3"
  fi
  cp "$F.orig" "$F"
}
run 'bad = [r for r in rows if r.get(column)' 'bad = [r for r in rows if True' "refused by header alone"
run '(not pooled and r[column] != r[first])' '(False and r[column] != r[first])' "second listing never compared"
run 'pooled = POOLED.search(column) is not None' 'pooled = False' "pooled sample IDs refused"
rm "$F.orig"
find "$WT/src" -name __pycache__ -type d -exec rm -rf {} +
