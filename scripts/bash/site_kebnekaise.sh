# /// codebase-meta
# exempt: sourced by fetchngs.sbatch and rnaseq.sbatch; says what one computer changes, runs no step
# ///
# What fetchngs.sbatch and rnaseq.sbatch do differently on Kebnekaise, chosen by CTX_SITE (this
# file is the default). Another Slurm site copies this file to site_<name>.sh and changes four
# things: the container profile (apptainer, or the site's nf-core/configs profile, which both
# pipelines load by name), the executor block, whether a run folder may sit in home, and the activation,
# whose environment file names the site's scratch variable. The #SBATCH lines stay in the
# launchers.

# Images run under Apptainer, which Kebnekaise has and Docker it does not.
SITE_FETCHNGS_PROFILE=singularity
SITE_RNASEQ_PROFILE=singularity
SITE_FETCHNGS_ARGS=()
# Home is 25 GB, and a run folder takes every FASTQ and image.
SITE_HOME_ALLOWED=0
# The index and the references published once per run folder, which later runs on the storage
# project can reuse.
SITE_SAVE_REFERENCE=true
# STAR on the genome, Salmon on its alignment and Salmon alone beside it (Rytis 2026-10-23,
# D-113): M2 and M3 read nf-core's own trimming, alignment and QC.
SITE_RNASEQ_ALIGN=true

# Nextflow from the fetch environment, every task a Slurm job (kebnekaise_env.sh).
site_activate() {
    source "$CTX_REPO/scripts/bash/kebnekaise_env.sh"
}

# Each task a Slurm job; HPC2N refused 100 jobs at once from one head job. Each named for its
# project and task (job_name.sh); $1 is the run folder's name.
site_executor() {
    printf 'executor {\n    queueSize = 50\n'
    job_nextflow_name "$1"
    printf '}\n'
}

# Nothing more: Slurm gives each task what nf-core/rnaseq asks for.
site_rnaseq_config() {
    :
}
