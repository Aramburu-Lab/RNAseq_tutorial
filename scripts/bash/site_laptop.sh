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
# sra-tools' fasterq-dump asks for 36 GB by default and failed on a 24 GB Mac.
SITE_FETCHNGS_ARGS=(--max_cpus 4 --max_memory 8.GB)
SITE_HOME_ALLOWED=1
# ctx reads none of the published references, and -resume reuses them from work/: saved, they
# were copied into results/ beside work/'s copy (U.10, 2026-10-08).
SITE_SAVE_REFERENCE=false

site_activate() {
    if [[ ! -x "$CTX_REPO/.pixi/envs/fetch/bin/nextflow" ]]; then
        echo "no nextflow in $CTX_REPO/.pixi/envs/fetch/bin; run pixi install -e fetch --frozen in $CTX_REPO" >&2
        return 1
    fi
    export PATH="$CTX_REPO/.pixi/envs/fetch/bin:$PATH" NXF_ANSI_LOG=false
}

# Docker's CPUs and memory, from one docker info, as "<cpus> <MB>": every task runs inside them.
_site_docker() {
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
    echo "$cpus $((bytes / 1048576))"
}

# The local executor sized to Docker, not to the host: sized to the host (14 CPUs, 24 GB), it ran
# two tasks each allowed all of Docker's 11.67 GiB side by side, and Docker killed one (exit 137,
# GTF2BED beside SALMON_INDEX, U.10).
site_executor() {
    local cpus mb
    read -r cpus mb <<<"$(_site_docker)" || return 1
    [[ -n "$mb" ]] || return 1
    printf 'executor {\n    queueSize = %s\n    cpus = %s\n    memory = %s.MB\n}\n' "$cpus" "$cpus" "$mb"
}

# A task may ask for no more than Docker has: nf-core/rnaseq's largest ask more than a laptop
# holds, and Nextflow refuses a task that asks for more than the executor.
site_rnaseq_config() {
    local cpus mb
    read -r cpus mb <<<"$(_site_docker)" || return 1
    [[ -n "$mb" ]] || return 1
    printf 'process {\n    resourceLimits = [cpus: %s, memory: %s.MB]\n}\n' "$cpus" "$mb"
}
