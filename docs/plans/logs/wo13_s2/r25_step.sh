#!/bin/bash
# wo13-s2 R2.5: one step of the guide's Part 1c for one study in a fresh root ~/ctx_$SLUG, its
# lines as the guide prints them with `pixi run -e laptop` replaced by this worktree's wrapper.
# The matrix and series matrix are copied from R2.2's root (~/ctx_u_<gse>), never moved; the
# guide's find, size and download steps are not run again, and nothing is deleted.
# Usage: r25_step.sh <gse lowercase> <stage|table|edit|setup|check|run>
set -uo pipefail
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-aa11ca42a7ace51a9
PIXI_LAPTOP="$WT/docs/plans/logs/wo13_s2_py.sh"
gse=$1 step=$2
SLUG=s2_$gse GSE=$(echo "$gse" | tr a-z A-Z)
OLD="$HOME/ctx_u_$gse/data/raw/u_$gse"
RAW="$HOME/ctx_$SLUG/data/raw/$SLUG"
FILE=$(cd "$OLD" && ls | grep -v series_matrix)
echo "## $step $GSE $SLUG $(date '+%Y-%m-%d %H:%M:%S')"
case $step in
  stage)
    if [ -e "$HOME/ctx_$SLUG" ]; then echo "exists: ~/ctx_$SLUG"; exit 1; fi
    mkdir -p "$RAW" && cp -p "$OLD/$FILE" "$OLD/${GSE}_series_matrix.txt.gz" "$RAW/"
    # R2.2's resources, copied: ctx setup checks each against its SHA-256 and fetches the rest.
    # R2.2's manifest stays behind, so this root's setup writes its own.
    mkdir -p "$HOME/ctx_$SLUG/resources"
    for d in gencode hgnc msigdb; do cp -Rp "$HOME/ctx_u_$gse/resources/$d" "$HOME/ctx_$SLUG/resources/"; done
    shasum -a 256 "$OLD"/* "$RAW"/*
    ;;
  table)
    cd "$WT" || exit 1
    "$PIXI_LAPTOP" python -m ctx.gather.geo_table "$RAW/${GSE}_series_matrix.txt.gz" ~/ctx_"$SLUG"/"$SLUG".csv
    cp conf/dataset_template_counts.toml ~/ctx_"$SLUG"/"$SLUG".toml
    ;;
  edit)
    sed 's/u_gse/s2_gse/g' "$WT/docs/plans/logs/wo13_u/r22_${gse}_edits.json" > "$WT/docs/plans/logs/wo13_s2/r25_${gse}_edits.json"
    "$PIXI_LAPTOP" python "$WT/docs/plans/logs/wo13_u/r22_edit.py" ~/ctx_"$SLUG"/"$SLUG".toml "$WT/docs/plans/logs/wo13_s2/r25_${gse}_edits.json"
    ;;
  setup)
    cd "$WT" || exit 1
    "$PIXI_LAPTOP" ctx --root ~/ctx_"$SLUG" setup --dataset-file ~/ctx_"$SLUG"/"$SLUG".toml | tail -8
    ;;
  check)
    cd "$WT" || exit 1
    "$PIXI_LAPTOP" ctx --root ~/ctx_"$SLUG" run --dataset "$SLUG" --dataset-file ~/ctx_"$SLUG"/"$SLUG".toml --counts "$RAW/$FILE" --check
    ;;
  run)
    cd "$WT" || exit 1
    /usr/bin/time -l "$PIXI_LAPTOP" ctx --root ~/ctx_"$SLUG" run --dataset "$SLUG" --dataset-file ~/ctx_"$SLUG"/"$SLUG".toml --counts "$RAW/$FILE"
    ;;
esac
echo "## exit $?"
