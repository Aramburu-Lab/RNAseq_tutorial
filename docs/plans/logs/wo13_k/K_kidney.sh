#!/bin/bash
# K.9's real check on ~/ctx_kidney (read only), under the worktree's .scratch_k folder, never
# ~/ctx_kit. First the whole kit, which D-072 refuses on this root today; then build_kit's own
# select() and carry_matrices() on it, which write LEFT_OUT.json (K_kidney_check.py).
W="$(cd "$(dirname "$0")/../../../.." && pwd)"
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
export PATH="$BIN:$PATH" PYTHONPATH="$W/src"
OUT="$W/.scratch_k/2026-10-15_course_kit_kidney_k9"
mkdir -p "$W/.scratch_k"
cd "$W" || exit 1
python scripts/python/build_kit.py --out "$OUT" --real kidney="$HOME/ctx_kidney" \
    --note "kidney=M9 on TCGA-KIRC shows the [clinical] genes only."
echo "build_kit exit $?"
python "$W/docs/plans/logs/wo13_k/K_kidney_check.py" "$HOME/ctx_kidney" \
    "$W/.scratch_k/kidney_left_out"
