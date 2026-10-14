#!/bin/bash
# Package H: the full suite once; its log in the scratch folder so the checkout stays untouched.
S=/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad
start=$(date +%s)
bash /Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a01b45b7d9228a7fd/docs/plans/logs/H_pytest.sh \
  -q -n auto --dist loadscope -rs -p no:cacheprovider > "$S/H_suite.log" 2>&1
echo "exit $? wall $(( $(date +%s) - start )) s" >> "$S/H_suite.log"
