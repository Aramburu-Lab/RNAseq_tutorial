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

# Each task's time limit at half of nf-core/rnaseq 3.27.0's conf/base.config (default and
# process_single and process_low 4 h, process_medium 8 h, process_high 16 h, process_long 20 h,
# each times task.attempt), so tasks wait less in the queue (Rytis 2026-10-23). A task past its
# half limit is retried once at twice it, the original, by the pipeline's own retry rule on exit
# 130 to 145. The head job's #SBATCH -t is not this.
site_rnaseq_config() {
    cat <<'EOF'
process {
    time = { 2.h * task.attempt }
    withLabel: process_single { time = { 2.h * task.attempt } }
    withLabel: process_low { time = { 2.h * task.attempt } }
    withLabel: process_medium { time = { 4.h * task.attempt } }
    withLabel: process_high { time = { 8.h * task.attempt } }
    withLabel: process_long { time = { 10.h * task.attempt } }
    // Measured (sacct, 2026-10-23, every completed task of the six aligned runs: 1,057 Trim
    // Galore, 61 STAR alignments, 6 index builds), with headroom (Rytis 2026-10-23, after the
    // exact maxima ran out of retries): time 2x the longest, memory 1.25x the largest peak,
    // cores as used; 25 % more memory and time per retry, cores never; three retries.
    withName: '.*:TRIMGALORE' { cpus = 3; memory = { 8000.MB * Math.pow(1.25, task.attempt - 1) }; time = { 33.min * Math.pow(1.25, task.attempt - 1) }; maxRetries = 3 }
    withName: '.*:FQ_LINT_AFTER_TRIMMING' { cpus = 1; memory = { 3520.MB * Math.pow(1.25, task.attempt - 1) }; time = { 12.min * Math.pow(1.25, task.attempt - 1) }; maxRetries = 3 }
    withName: '.*:FQ_LINT' { cpus = 1; memory = { 1984.MB * Math.pow(1.25, task.attempt - 1) }; time = { 11.min * Math.pow(1.25, task.attempt - 1) }; maxRetries = 3 }
    withName: '.*:FASTQC' { cpus = 2; memory = { 6208.MB * Math.pow(1.25, task.attempt - 1) }; time = { 35.min * Math.pow(1.25, task.attempt - 1) }; maxRetries = 3 }
    withName: '.*:SALMON_QUANT' { cpus = 6; memory = { 27328.MB * Math.pow(1.25, task.attempt - 1) }; time = { 31.min * Math.pow(1.25, task.attempt - 1) }; maxRetries = 3 }
    withName: '.*:FQ_SUBSAMPLE' { cpus = 1; memory = { 1473.MB * Math.pow(1.25, task.attempt - 1) }; time = { 11.min * Math.pow(1.25, task.attempt - 1) }; maxRetries = 3 }
    withName: '.*:STAR_ALIGN' { cpus = 11; memory = { 71616.MB * Math.pow(1.25, task.attempt - 1) }; time = { 159.min * Math.pow(1.25, task.attempt - 1) }; maxRetries = 3 }
    withName: '.*:SAMTOOLS_SORT' { cpus = 5; memory = { 16704.MB * Math.pow(1.25, task.attempt - 1) }; time = { 9.min * Math.pow(1.25, task.attempt - 1) }; maxRetries = 3 }
    withName: '.*:SAMTOOLS_INDEX' { cpus = 2; memory = { 1088.MB * Math.pow(1.25, task.attempt - 1) }; time = { 61.min * Math.pow(1.25, task.attempt - 1) }; maxRetries = 3 }
    withName: '.*:CAT_FASTQ' { cpus = 1; memory = { 6336.MB * Math.pow(1.25, task.attempt - 1) }; time = { 2.min * Math.pow(1.25, task.attempt - 1) }; maxRetries = 3 }
    withName: '.*:STAR_GENOMEGENERATE' { cpus = 6; memory = { 85824.MB * Math.pow(1.25, task.attempt - 1) }; time = { 111.min * Math.pow(1.25, task.attempt - 1) }; maxRetries = 3 }
    withName: '.*:GUNZIP_FASTA' { cpus = 1; memory = { 2880.MB * Math.pow(1.25, task.attempt - 1) }; time = { 2.min * Math.pow(1.25, task.attempt - 1) }; maxRetries = 3 }
    withName: '.*:GUNZIP_GTF' { cpus = 1; memory = { 4160.MB * Math.pow(1.25, task.attempt - 1) }; time = { 2.min * Math.pow(1.25, task.attempt - 1) }; maxRetries = 3 }
    withName: '.*:SAMTOOLS_FAIDX' { cpus = 1; memory = { 832.MB * Math.pow(1.25, task.attempt - 1) }; time = { 1.min * Math.pow(1.25, task.attempt - 1) }; maxRetries = 3 }
    withName: '.*:CUSTOM_GTFFILTER' { cpus = 1; memory = { 4672.MB * Math.pow(1.25, task.attempt - 1) }; time = { 2.min * Math.pow(1.25, task.attempt - 1) }; maxRetries = 3 }
    withName: '.*:MAKE_TRANSCRIPTS_FASTA' { cpus = 1; memory = { 8128.MB * Math.pow(1.25, task.attempt - 1) }; time = { 6.min * Math.pow(1.25, task.attempt - 1) }; maxRetries = 3 }
    withName: '.*:EAUTILS_GTF2BED' { cpus = 1; memory = { 14784.MB * Math.pow(1.25, task.attempt - 1) }; time = { 4.min * Math.pow(1.25, task.attempt - 1) }; maxRetries = 3 }
    withName: '.*:SALMON_INDEX' { cpus = 3; memory = { 33728.MB * Math.pow(1.25, task.attempt - 1) }; time = { 55.min * Math.pow(1.25, task.attempt - 1) }; maxRetries = 3 }
    withName: '.*:SE_TRANSCRIPT_UNIFIED' { cpus = 1; memory = { 9408.MB * Math.pow(1.25, task.attempt - 1) }; time = { 23.min * Math.pow(1.25, task.attempt - 1) }; maxRetries = 3 }
    withName: '.*:SE_GENE_UNIFIED' { cpus = 1; memory = { 3904.MB * Math.pow(1.25, task.attempt - 1) }; time = { 8.min * Math.pow(1.25, task.attempt - 1) }; maxRetries = 3 }
    // Preseq in defect mode from the first try: 3.27.0's retry, meant to add -defects, runs
    // `args.join(' -defects')`, which is Java's static String.join with no elements and gives ''
    // (printed '[]' under Nextflow 26.04.6), so the retry drops every argument and fails again;
    // Pancreas-tissue-PAAD stopped on it 2026-10-24. The rest is nf-core's own ext.args, verbatim
    // from conf/modules/preseq.config; the same selector, so this file's value is the one used.
    // No retry at all: every retry runs with no arguments, -bam among them, and fails on "reads
    // unsorted", hiding the first try's own cause (on 2026-10-25 the full disk killed it). A
    // failed Preseq ends the run on its real error, and a resume runs it again whole. Not
    // maxRetries = 0: Nextflow 26.04.6 still retried once with it (a toy pipeline, 2026-10-25).
    // 32 GB on its one try: every first try on Pancreas-tissue-PAAD was killed OUT_OF_MEMORY at
    // nf-core's 6 GB (sacct, 2026-10-25, 10 tasks, 1 to 18 min in), so its real peak is not
    // measured; a choice to replace with 1.25x the peak sacct gives once one completes.
    withName: 'PRESEQ_LCEXTRAP' { ext.args = '-verbose -bam -seed 1 -seg_len 100000000 -defects'; errorStrategy = 'terminate'; memory = 32.GB }
}
EOF
}
