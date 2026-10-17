#!/bin/bash
# Package Q: a scratch ctx root whose gather outputs are symlinks to ~/ctx_lung's newest
# tcga_lung gather (2026-10-16); ctx qc writes only under the scratch root's outputs/02_qc.
WT=/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a7158d7cd646891a9
ROOT=$WT/.scratch_q/root
mkdir -p "$ROOT/outputs/01_gather"
for f in "$HOME"/ctx_lung/outputs/01_gather/2026-10-16_gather_tcga_lung_*; do
  ln -sf "$f" "$ROOT/outputs/01_gather/$(basename "$f")"
done
ls -l "$ROOT/outputs/01_gather" | head -3
