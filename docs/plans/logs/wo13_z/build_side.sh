#!/bin/bash
# wo13-z: build.py's groups on one side, three lanes in parallel, each lane's log beside this file.
D=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a6e30099a3085f55d/docs/plans/logs/wo13_z
PY=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin/python
side=$1
"$PY" "$D/build.py" "$side" tier1 hnsc s2_gse327166 s2_gse327906 r2_gse130605 > "$D/build_${side}_lane1.log" 2>&1 &
"$PY" "$D/build.py" "$side" breast lung > "$D/build_${side}_lane2.log" 2>&1 &
"$PY" "$D/build.py" "$side" kidney > "$D/build_${side}_lane3.log" 2>&1 &
wait
echo "side $side done"
