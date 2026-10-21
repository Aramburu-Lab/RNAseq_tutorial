#!/bin/bash
# wo13-z1: local M10 cohorts into scratch roots with one src. Usage: build_all.sh <side>
# <src dir> <slug> [slug ...]; side names the roots, e.g. before (HEAD's src) or after.
# macOS bash 3.2: no associative arrays, so the real root comes from a case.
side=$1 src=$2
shift 2
HERE=$(dirname "$0")
ROOTS=/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/roots
for s in "$@"; do
  case $s in
    tcga_brca) real=ctx_breast ;;
    tcga_kirc) real=ctx_kidney ;;
    tcga_luad | tcga_lusc) real=ctx_lung ;;
    tcga_hnsc) real=ctx_tcga_hnsc ;;
    *) echo "no root for $s"; continue ;;
  esac
  bash "$HERE/build.sh" "/Users/rytis/$real" "$s" "$src" "$ROOTS/${side}_${s#tcga_}"
done
