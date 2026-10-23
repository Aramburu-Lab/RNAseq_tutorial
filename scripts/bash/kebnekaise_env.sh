# /// codebase-meta
# exempt: sourced by the cluster scripts beside it and by a login shell; runs no pipeline step
# ///
# Nextflow for the course's cluster jobs on Kebnekaise, from the repository's pinned fetch
# environment, with every cache on the student's storage project. A login shell sources it as
# well as a job, so a failure returns 1 instead of exiting, and a refusal leaves the caller's
# shell options, PATH and prompt as they were, with no helper variable behind. The Usage lines
# run once, on a login node; every later login needs the two export lines again before any
# sbatch.
#
# Usage:
#   export CTX_STORE=<your folder on the storage project>
#   export CTX_REPO="$CTX_STORE/RNAseq_tutorial"
#   mkdir -p "$CTX_STORE"
#   git clone <course repository> "$CTX_REPO"
#   curl -fsSL https://pixi.sh/install.sh | PIXI_VERSION=v0.76.2 PIXI_NO_PATH_UPDATE=1 bash
#   PIXI_CACHE_DIR="$CTX_STORE/pixi_cache" "${PIXI_HOME:-$HOME/.pixi}/bin/pixi" install --manifest-path "$CTX_REPO/pixi.toml" -e fetch --frozen
#   source "$CTX_REPO/scripts/bash/kebnekaise_env.sh"
#   nextflow -version

# Run as a program, every return below is a no-op and the checks fall through; this exit ends
# only that child process.
[[ "${BASH_SOURCE[0]}" != "$0" ]] || { echo "kebnekaise_env.sh: source this file, do not run it" >&2; exit 1; }

CTX_REPO="${CTX_REPO:-}"
CTX_STORE="${CTX_STORE:-}"
PIXI="${PIXI:-${PIXI_HOME:-$HOME/.pixi}/bin/pixi}"

# Every refusal passes through here, so none leaves a helper in the student's shell.
_ctx_fail() {
    echo "kebnekaise_env.sh: $*" >&2
    unset _ctx_store _ctx_home _ctx_hook _ctx_want _ctx_nf _ctx_jvm
    unset -f _ctx_fail
}

if [[ -z "$CTX_REPO" || -z "$CTX_STORE" ]]; then
    _ctx_fail "export CTX_REPO and CTX_STORE first (this file's Usage)"
    return 1
fi
if ! _ctx_store=$(CDPATH='' cd -- "$CTX_STORE" 2>/dev/null && pwd -P); then
    _ctx_fail "CTX_STORE=$CTX_STORE is not a folder; mkdir -p it first"
    return 1
fi
# Canonical paths on both sides: either may be reached through a symlink.
_ctx_home=$(CDPATH='' cd -- "${HOME:-/nonexistent}" 2>/dev/null && pwd -P) ||
    _ctx_home="${HOME:-/nonexistent}"
case "$_ctx_store/" in
"$_ctx_home"/*)
    _ctx_fail "CTX_STORE=$CTX_STORE is in your home folder (25 GB), which" \
        "Nextflow's images fill; use your folder on the storage project"
    return 1
    ;;
esac
# pixi 0.76.2's shell-hook --as-is exits 0 with a full hook for an environment never installed;
# evaluated, it would change PATH and the prompt before the check after it could refuse.
if [[ ! -x "$CTX_REPO/.pixi/envs/fetch/bin/nextflow" ]]; then
    _ctx_fail "no nextflow in $CTX_REPO/.pixi/envs/fetch/bin; install the fetch environment" \
        "(this file's Usage)"
    return 1
fi

# The lab's login environment exports NXF_ASSETS; a set NXF_VER makes the launcher fetch another
# Nextflow; NXF_JAVA_HOME overrides the environment's Java. module purge clears none of them.
unset NXF_JAVA_HOME NXF_VER NXF_ASSETS
# The hook is captured before it is evaluated: an empty hook would leave a login Nextflow on
# PATH. --as-is because a job never installs; the login node did that once. The prompt stays the
# student's own.
if ! _ctx_hook=$("$PIXI" shell-hook --shell bash --manifest-path "$CTX_REPO/pixi.toml" -e fetch \
    --as-is --change-ps1 false); then
    _ctx_fail "'$PIXI shell-hook' failed; install pixi and the fetch environment" \
        "(this file's Usage)"
    return 1
fi
eval "$_ctx_hook"
# Folders compared canonically: pixi prints the physical path of a storage project entered
# through its symlink, so the two strings differ for the same file.
_ctx_want=$(CDPATH='' cd -- "$CTX_REPO/.pixi/envs/fetch/bin" 2>/dev/null && pwd -P || true)
_ctx_nf=$(command -v nextflow || true)
if [[ -z "$_ctx_want" || -z "$_ctx_nf" ||
    "$(CDPATH='' cd -- "${_ctx_nf%/*}" 2>/dev/null && pwd -P)" != "$_ctx_want" ]]; then
    _ctx_fail "nextflow is ${_ctx_nf:-not on PATH}, not the one in" \
        "$CTX_REPO/.pixi/envs/fetch/bin; install the fetch environment (this file's Usage)"
    return 1
fi

# Assets and images off the home folder; the head job has 8 GB.
export NXF_HOME="$_ctx_store/nextflow" NXF_OPTS='-Xms500M -Xmx4G'
# Plain logs in slurm-*.out, and every task its own Slurm job rather than a process on the head job.
export NXF_ANSI_LOG=false NXF_EXECUTOR=slurm
# One image folder under both names: a run that reads only the other pulls every image again.
export NXF_SINGULARITY_CACHEDIR="$_ctx_store/images" NXF_APPTAINER_CACHEDIR="$_ctx_store/images"
# Apptainer's blob cache outside the image folder (nested, one reached 222 GB), and its scratch
# node-local: unpacking an image on Lustre fails with unpriv.lsetxattr.
export APPTAINER_CACHEDIR="$_ctx_store/apptainer_cache" APPTAINER_TMPDIR="${SNIC_TMP:-/tmp}"
# Child jobs are billed to the compute project the head job was submitted with (sbatch -A).
if [[ -n "${SLURM_JOB_ACCOUNT:-}" ]]; then
    export SBATCH_ACCOUNT="$SLURM_JOB_ACCOUNT"
fi
if ! mkdir -p "$NXF_HOME" "$NXF_SINGULARITY_CACHEDIR" "$APPTAINER_CACHEDIR"; then
    _ctx_fail "cannot create the cache folders under $CTX_STORE"
    return 1
fi
# In a job, Java runs from a copy on the node's own disk. Run from Lustre, its code could not be
# read back on 2026-10-23: five head jobs on one node died of SIGBUS (BUS_ADRERR at the program
# counter, in libjvm.so) while the file stood unchanged since 2026-10-22 (NF.13). A copy that
# fails or does not start leaves Nextflow on the environment's own Java, said in the job's log.
_ctx_jvm="$CTX_REPO/.pixi/envs/fetch/lib/jvm"
if [[ -n "${SLURM_JOB_ID:-}" && -n "${SNIC_TMP:-}" && -d "$_ctx_jvm" ]]; then
    if mkdir -p "$SNIC_TMP/ctx_jvm" && cp -a "$_ctx_jvm/." "$SNIC_TMP/ctx_jvm/" &&
        "$SNIC_TMP/ctx_jvm/bin/java" -version >/dev/null 2>&1; then
        export NXF_JAVA_HOME="$SNIC_TMP/ctx_jvm"
        echo "kebnekaise_env.sh: Java from the node's disk, $NXF_JAVA_HOME"
    else
        echo "kebnekaise_env.sh: Java not copied to $SNIC_TMP; it runs from $_ctx_jvm" >&2
    fi
fi
unset _ctx_store _ctx_home _ctx_hook _ctx_want _ctx_nf _ctx_jvm
unset -f _ctx_fail
