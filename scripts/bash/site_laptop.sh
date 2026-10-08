# /// codebase-meta
# exempt: sourced by fetchngs.sbatch and rnaseq.sbatch; says what one computer changes, runs no step
# ///
# What fetchngs.sbatch and rnaseq.sbatch do differently on a laptop (CTX_SITE=laptop): every task
# runs here under Docker, Nextflow comes from the repository's fetch environment, and the run
# folder may sit in home. Run the launchers with bash, not sbatch; Docker must be running.

# nf-core's arm profiles: fetchngs 1.12.0 names it arm, rnaseq 3.27.0 arm64; both run amd64
# images under emulation where no arm image exists.
case "$(uname -m)" in
arm64 | aarch64) SITE_FETCHNGS_PROFILE=docker,arm SITE_RNASEQ_PROFILE=docker,arm64 ;;
*) SITE_FETCHNGS_PROFILE=docker SITE_RNASEQ_PROFILE=docker ;;
esac
# Nextflow's own default for local tasks: one per core.
SITE_QUEUE_SIZE=$(getconf _NPROCESSORS_ONLN)
# sra-tools' fasterq-dump asks for 36 GB by default and failed on a 24 GB Mac.
SITE_FETCHNGS_ARGS=(--max_cpus 4 --max_memory 8.GB)
SITE_HOME_ALLOWED=1

site_activate() {
    if [[ ! -x "$CTX_REPO/.pixi/envs/fetch/bin/nextflow" ]]; then
        echo "no nextflow in $CTX_REPO/.pixi/envs/fetch/bin; run pixi install -e fetch --frozen in $CTX_REPO" >&2
        return 1
    fi
    export PATH="$CTX_REPO/.pixi/envs/fetch/bin:$PATH" NXF_ANSI_LOG=false
}

# A task may ask for no more than Docker has: nf-core/rnaseq's largest ask more than a laptop
# holds, and Nextflow refuses a task that asks for more than the machine.
site_rnaseq_config() {
    local found cpus bytes
    if ! found=$(docker info --format '{{.NCPU}} {{.MemTotal}}' 2>/dev/null) || [[ -z "$found" ]]; then
        echo "docker info failed: start Docker first" >&2
        return 1
    fi
    read -r cpus bytes <<<"$found"
    [[ "$cpus" =~ ^[1-9][0-9]*$ && "$bytes" =~ ^[1-9][0-9]*$ ]] || {
        echo "docker info gave $found, not CPUs and bytes" >&2
        return 1
    }
    printf 'process {\n    resourceLimits = [cpus: %s, memory: %s.MB]\n}\n' "$cpus" "$((bytes / 1048576))"
}
