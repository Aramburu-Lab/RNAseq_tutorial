#!/bin/bash
# K.9's real check: a kit from ~/ctx_kidney alone (read only), built under the worktree's
# .scratch_k folder, never ~/ctx_kit; then check_kit on it and the LEFT_OUT.json commands.
W="$(cd "$(dirname "$0")/../../../.." && pwd)"
BIN=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin
export PATH="$BIN:$PATH" PYTHONPATH="$W/src"
OUT="$W/.scratch_k/2026-10-15_course_kit_kidney_k9"
mkdir -p "$W/.scratch_k"
cd "$W" || exit 1
/usr/bin/time -l python scripts/python/build_kit.py --out "$OUT" --real kidney="$HOME/ctx_kidney" \
    --note "kidney=M9 on TCGA-KIRC shows the [clinical] genes only." 2>&1 \
    | grep -E "course_kit|real|maximum resident|refused|fail|^/|left out"
python "$OUT/check_kit.py" "$OUT"
python "$W/docs/plans/logs/wo13_k/K_kidney_check.py" "$OUT"
python -c "import json,sys; print(json.load(open(sys.argv[1]))['roots'])" "$OUT/MANIFEST.json"
grep -n "kidney/\|builder's" "$OUT/README.md"
