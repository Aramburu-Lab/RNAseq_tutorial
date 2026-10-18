#!/bin/bash
# wo13-u R2.3: every command that loads a dataset file, on a file that is not valid TOML.
S=/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/u
PY=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aa1d340947c025c9d/docs/plans/logs/wo13_u_py.sh
mkdir -p "$S"
printf '[dataset]\nslug = "x"\nlabel = \n' > "$S/bad.toml"
for c in "run --dataset x --dataset-file $S/bad.toml --check" \
         "gather --dataset x --counts $S/bad.toml --dataset-file $S/bad.toml" \
         "qc --dataset x --dataset-file $S/bad.toml" "de --dataset x --dataset-file $S/bad.toml" \
         "function --dataset x --dataset-file $S/bad.toml" \
         "ncrna --dataset x --dataset-file $S/bad.toml" "clinical --dataset x --dataset-file $S/bad.toml" \
         "proteo --dataset x --dataset-file $S/bad.toml" "genecard --dataset x --dataset-file $S/bad.toml" \
         "process --dataset x --dataset-file $S/bad.toml --genome $S" \
         "quant --dataset x --dataset-file $S/bad.toml --rnaseq-dir $S" \
         "report --dataset x --dataset-file $S/bad.toml" \
         "setup --dataset-file $S/bad.toml" "guide --dataset-file $S/bad.toml"; do
  # shellcheck disable=SC2086
  "$PY" python -m ctx --root "$S/root" $c > "$S/out.txt" 2>&1
  echo "== ctx $c: exit $?, traceback lines $(grep -c 'Traceback\|╭' "$S/out.txt")"
  tail -3 "$S/out.txt"
done
