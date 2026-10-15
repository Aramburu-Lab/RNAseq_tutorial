#!/bin/bash
# Package P, GS.3: the unthinned volcano's pages. Step-1 (de, function, report) and TCGA-BRCA
# from the gs3 src; LUSC and CPTAC ccRCC, the two largest M5 pages, gathered and built from
# main's src and from the gs3 src, so each page is measured both ways.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a4bc57162a0d35df8
B=$WT/docs/plans/logs/wo13_p/build.sh
MAIN=$WT/outputs/scratch_roots/src_main/src GS3=$WT/outputs/scratch_roots/src_gs3/src
bash "$B" gs3 "$GS3" tier1 de function report &
bash "$B" gs3 "$GS3" breast de &
wait
bash "$B" base "$MAIN" lusc gather &
bash "$B" base "$MAIN" ccrcc gather &
wait
bash "$B" before "$MAIN" lusc de &
bash "$B" before "$MAIN" ccrcc de &
wait
bash "$B" gs3 "$GS3" lusc de &
bash "$B" gs3 "$GS3" ccrcc de &
wait
echo "gs3 builds done"
