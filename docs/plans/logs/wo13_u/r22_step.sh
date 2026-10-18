#!/bin/bash
# wo13-u R2.2: one step of the guide's Part 1c for one study, its lines as the guide prints
# them with `pixi run -e laptop` replaced by "$PIXI_LAPTOP" (r22_env.sh).
# Usage: r22_step.sh <GSE> <slug> <find|size|download|table|check|run> [FILE] [URL]
set -uo pipefail
here=$(dirname "$0")
# shellcheck source=/dev/null
source "$here/r22_env.sh" "$1" "$2"
step=$3
FILE=${4:-} URL=${5:-}
echo "## $step $GSE $SLUG $(date '+%Y-%m-%d %H:%M:%S')"
case $step in
  find)
    mkdir -p "$RAW" && cd "$RAW" || exit 1
    curl -fsS "https://www.ncbi.nlm.nih.gov/geo/download/?acc=$GSE" | grep -o 'file=[^"&]*' | sort -u
    curl -fsS "https://ftp.ncbi.nlm.nih.gov/geo/series/${GSE%???}nnn/$GSE/suppl/" | grep -o 'href="[^"?/]*"' | sort -u
    ;;
  size)
    curl -fsSIL "$URL" | grep -i content-length
    ;;
  download)
    cd "$RAW" && curl -fL -o "$FILE" "$URL"
    curl -fL -O "https://ftp.ncbi.nlm.nih.gov/geo/series/${GSE%???}nnn/$GSE/matrix/${GSE}_series_matrix.txt.gz"
    ls -l "$RAW"
    ;;
  table)
    cd "$CTX_REPO" || exit 1
    "$PIXI_LAPTOP" python -m ctx.gather.geo_table "$RAW/${GSE}_series_matrix.txt.gz" ~/ctx_"$SLUG"/"$SLUG".csv "${@:6}"
    [ -e ~/ctx_"$SLUG"/"$SLUG".toml ] || cp conf/dataset_template_counts.toml ~/ctx_"$SLUG"/"$SLUG".toml
    ;;
  check)
    cd "$CTX_REPO" || exit 1
    "$PIXI_LAPTOP" ctx --root ~/ctx_"$SLUG" setup --dataset-file ~/ctx_"$SLUG"/"$SLUG".toml | tail -4
    "$PIXI_LAPTOP" ctx --root ~/ctx_"$SLUG" run --dataset "$SLUG" --dataset-file ~/ctx_"$SLUG"/"$SLUG".toml --counts "$RAW/$FILE" --check
    ;;
  run)
    cd "$CTX_REPO" || exit 1
    /usr/bin/time -l "$PIXI_LAPTOP" ctx --root ~/ctx_"$SLUG" run --dataset "$SLUG" --dataset-file ~/ctx_"$SLUG"/"$SLUG".toml --counts "$RAW/$FILE"
    ;;
esac
echo "## exit $?"
