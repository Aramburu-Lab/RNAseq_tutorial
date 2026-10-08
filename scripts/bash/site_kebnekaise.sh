# /// codebase-meta
# exempt: sourced by fetchngs.sbatch and rnaseq.sbatch; says what one computer changes, runs no step
# ///
# What fetchngs.sbatch and rnaseq.sbatch do differently on Kebnekaise, chosen by CTX_SITE (this
# file is the default). Another Slurm site copies this file to site_<name>.sh and changes four
# things: the container profile (apptainer, or the site's nf-core/configs profile, which both
# pipelines load by name), queueSize, whether a run folder may sit in home, and the activation,
# whose environment file names the site's scratch variable. The #SBATCH lines stay in the
# launchers.

# Images run under Apptainer, which Kebnekaise has and Docker it does not.
SITE_FETCHNGS_PROFILE=singularity
SITE_RNASEQ_PROFILE=singularity
# HPC2N refused 100 jobs at once from one head job.
SITE_QUEUE_SIZE=50
SITE_FETCHNGS_ARGS=()
# Home is 25 GB, and a run folder takes every FASTQ and image.
SITE_HOME_ALLOWED=0

# Nextflow from the fetch environment, every task a Slurm job (kebnekaise_env.sh).
site_activate() {
    source "$CTX_REPO/scripts/bash/kebnekaise_env.sh"
}

# Nothing beyond queueSize: Slurm gives each task what nf-core/rnaseq asks for.
site_rnaseq_config() {
    :
}
