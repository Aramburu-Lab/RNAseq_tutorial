#!/bin/bash
# Package M9b: a figure PDF of a scratch root as PNG for an eye check, into the scratchpad.
#   bash render.sh <root label> <stem part, e.g. contrasts_isg15> <out name>
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-af0224262630483d8
SCRATCH=/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad
pdf=$(ls "$WT"/outputs/scratch_roots/m9b_$1/outputs/11_genecard/*_"$2"_*v1_white.pdf | tail -1)
sips -s format png -s dpiWidth 144 -s dpiHeight 144 "$pdf" --out "$SCRATCH/m9b_$3.png" > /dev/null
echo "$pdf"
