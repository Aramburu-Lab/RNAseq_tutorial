#!/usr/bin/env bash
# wo13-r1: the new fetch script against the kidney project's existing folders, fetching nothing.
# 1. A scratch copy of each folder (records copied, data files hard-linked, never written) is run
#    through the script with GDC's API replaced by a local status file: every record is kept,
#    every data file's md5 matches, so nothing is fetched, and the end check must pass.
# 2. An empty folder run against a local API whose endpoints return one byte: the script saves
#    each query before it posts, and each must equal the existing folder's byte for byte.
# CTX_REPO is main's checkout for its fetch environment; its site_laptop.sh equals this branch's.
set -uo pipefail
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a00f69bc2d903ff45
SCRIPT=$WT/scripts/bash/fetch_gdc.sbatch
SRC=$HOME/ctx_kidney/data/raw
S=/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/r1_layout
export CTX_REPO=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial CTX_SITE=laptop
[[ ! -e "$S" ]] || { echo "move $S aside first"; exit 1; }
mkdir -p "$S/api" "$S/kept" "$S/empty" "$S/api1"
cp "$SRC/tcga_kirc/gdc_status.json" "$S/api/status"
cp "$SRC/tcga_kirc/gdc_status.json" "$S/api1/status"
for e in files cases annotations; do printf 'x' >"$S/api1/$e"; done
for f in tcga_kirc tcga_kirc_mirna tcga_kirc_isoform; do
    mkdir -p "$S/kept/$f"
    find "$SRC/$f" -maxdepth 1 -type f -exec cp -p {} "$S/kept/$f/" \;
    for d in "$SRC/$f"/*/; do
        id=$(basename "$d")
        mkdir -p "$S/kept/$f/$id"
        ln "$d"* "$S/kept/$f/$id/"
    done
done
listing() { (cd "$1" && find . -type f -not -path './*/*/*' | sort | while read -r p; do printf '%s %s\n' "$(md5 -q "$p")" "$p"; done); }
for f in tcga_kirc tcga_kirc_mirna tcga_kirc_isoform; do listing "$S/kept/$f" >"$S/before_$f.txt"; done
echo "== 1. existing folders, nothing fetched"
GDC_API="file://$S/api" GDC_WEB="file://$S/api" bash "$SCRIPT" TCGA-KIRC "$S/kept/tcga_kirc"
echo "exit $?"
for f in tcga_kirc tcga_kirc_mirna tcga_kirc_isoform; do
    listing "$S/kept/$f" >"$S/after_$f.txt"
    echo "-- $f: files added, removed or changed by the run"
    diff "$S/before_$f.txt" "$S/after_$f.txt" | grep '^[<>]' | awk '{print $1, $3}'
    echo "-- $f: retrieval rows before, after: $(wc -l <"$SRC/$f/_retrieval_log.tsv") $(wc -l <"$S/kept/$f/_retrieval_log.tsv")"
done
echo "== 2. the queries of an empty folder against the existing folders' queries"
GDC_API="file://$S/api1" GDC_WEB="file://$S/api1" RETRIES=0 bash "$SCRIPT" TCGA-KIRC "$S/empty/tcga_kirc" >"$S/empty.log" 2>&1
echo "exit $? (1 expected: the stand-in endpoints return no manifest)"
for f in tcga_kirc tcga_kirc_mirna tcga_kirc_isoform; do
    for q in "$SRC/$f"/*.query.json; do
        n=$(basename "$q")
        if cmp -s "$q" "$S/empty/$f/$n"; then echo "equal $f/$n"; else echo "DIFFERS $f/$n"; fi
    done
    for q in "$S/empty/$f"/*.query.json; do
        [[ -f "$SRC/$f/$(basename "$q")" ]] || echo "NEW $f/$(basename "$q")"
    done
done
