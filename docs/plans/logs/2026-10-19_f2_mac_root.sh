#!/usr/bin/env bash
# F2's cross-machine proof: the fixture root drawn on this Mac (osx-arm64, main's pixi default
# environment), ctx imported from a copy of src outside any checkout, as an installed ctx is;
# then the volcano's Docker tab redraws it in the figures image.
set -euo pipefail
W=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-ad214d746812b3316
S=/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad
export PATH=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin:$PATH
export PYTHONPATH=$S/f2_src
cd "$W"
python -c 'import ctx, ctx.runtime as r; print(ctx.__file__, r.source())'
/usr/bin/time -l python containers/figures_root.py "$S/f2_mac/root"
/usr/bin/time -l /usr/bin/python3 containers/prove_figures.py --image ctx-figures:f2 \
  --work "$S/f2_mac/proof" --root "$S/f2_mac/root"
