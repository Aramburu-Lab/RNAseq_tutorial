# /// codebase-meta
# exempt: sourced by every launcher; names Slurm jobs, runs no step
# ///
# Slurm job names (Rytis 2026-10-23): <ANALYSIS>_<PROJECT>_<PROCESS>_<VERSION>, four fields joined
# by _, none holding _ itself, so a queue listing splits into exactly four. Inside a field a colon
# and parentheses are dropped, any other character outside A-Za-z0-9.- becomes - and a run of -
# is one: the rule that made each dataset's name from its label (D-114). PROJECT is the label of
# conf/datasets/<run folder>.toml, which for a course dataset reads as the folder's own name
# (Breast-tissue-BRCA); a folder with no dataset file shows its own name. PROCESS is nf-head
# for a Nextflow head job, nf-<task> for its tasks, the tag kept without its parentheses, or the
# launcher's own step; VERSION is NA, since no run folder carries a release number. The
# #SBATCH -J line can only say NA for PROJECT: Slurm reads it before the arguments, so a pending
# job shows NA until it starts and renames itself.
JOB_ANALYSIS=RNAseq-Tutorial

job_field() {
    local s=${1:-NA}
    s=${s//[:()]/}
    s=${s//[^A-Za-z0-9.-]/-}
    while [[ "$s" == *--* ]]; do s=${s//--/-}; done
    printf '%s' "$s"
}

# A run folder's PROJECT: its dataset's label where the course has a dataset file of that name.
job_project() {
    local file="${CTX_REPO:-}/conf/datasets/$1.toml" label=""
    [[ -f "$file" ]] && label=$(sed -n 's/^label = "\(.*\)"$/\1/p' "$file" | head -n 1)
    printf '%s' "${label:-$1}"
}

job_name() {
    printf '%s_%s_%s_%s' "$JOB_ANALYSIS" "$(job_field "$1")" "$(job_field "$2")" "$(job_field "${3:-NA}")"
}

# The running job renamed for its project; a failed rename only warns.
job_rename() {
    [[ -n "${SLURM_JOB_ID:-}" ]] && command -v scontrol >/dev/null || return 0
    scontrol update JobId="$SLURM_JOB_ID" JobName="$(job_name "$@")" ||
        echo "warning: job $SLURM_JOB_ID keeps its name" >&2
}

# The executor line that names each Nextflow task the same way: the last part of its task name,
# by job_field's rule.
job_nextflow_name() {
    printf '    jobName = { "%s_%s_nf-" + task.name.tokenize('"':'"')[-1].replaceAll(/[()]/, '"''"').replaceAll(/[^A-Za-z0-9.-]/, '"'-'"').replaceAll(/-+/, '"'-'"') + "_%s" }\n' \
        "$JOB_ANALYSIS" "$(job_field "$1")" "$(job_field "${2:-NA}")"
}
