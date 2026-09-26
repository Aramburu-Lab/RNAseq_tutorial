#!/usr/bin/env bash
# /// codebase-meta
# exempt: >-
#   Build tooling, not a pipeline step. Its artifact is a registry digest recorded in
#   containers/digests.lock, not a file under outputs/, and it runs on a laptop rather than as
#   part of any analysis.
# ///
# Build the analysis image. Two targets, two registries, two different jobs (D-004).
#
# Usage:
#   containers/build.sh dev     [--tag TAG] [--push] [--no-cache] [--ssh]
#   containers/build.sh release --version X.Y.Z [--no-cache] [--ssh] [--skip-sif]
#
#   dev      fast loop. Builds on this machine, pushes to ${DEV_REGISTRY} (default
#            docker.io/stakaitis). Consumed with `apptainer pull docker://…` on the login node.
#            NOT reproducible, and that is fine — never cite a dev tag in a manuscript.
#   release  built by CI on native linux-64. Pushes the OCI image to GHCR, converts it once with
#            `apptainer build --reproducible`, pushes that .sif as an ORAS artifact, and records
#            both digests in containers/digests.lock (ADR-7).
#
#   --tag       dev tag, default: short git SHA, or `dirty-<SHA>` with uncommitted changes
#   --version   release version, e.g. 0.1.0. Written to the image tag and to digests.lock
#   --push      dev only. release always pushes; there is nothing to convert otherwise
#   --no-cache  pass --no-cache to the docker build
#   --ssh       forward the ssh-agent (`--ssh default`) so a PRIVATE labcore can be cloned during
#               the build. See containers/README.md; omit it when labcore is public
#   --skip-sif  release only. Push the OCI image, skip the .sif leg. Escape hatch for a runner
#               without apptainer — the digests.lock it writes is then INCOMPLETE
#
# The Mac is arm64 and Kebnekaise is x86_64. A default `docker build` here produces linux/arm64,
# which `apptainer pull` accepts without complaint and which then fails to exec on a compute node
# with a format error that names neither docker nor the architecture. So the platform is hardcoded
# below and every push is gated on the built image really being amd64.

set -euo pipefail

readonly PLATFORM="linux/amd64"
readonly IMAGE_NAME="cancer_transcriptomics_teaching-analysis"
readonly GHCR_REPO="ghcr.io/aramburu-lab/${IMAGE_NAME}"
readonly SIF_REPO="ghcr.io/aramburu-lab/${IMAGE_NAME}-sif"
DEV_REGISTRY="${DEV_REGISTRY:-docker.io/stakaitis}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
readonly DOCKERFILE="containers/analysis.Dockerfile"
readonly LOCKFILE="${REPO_ROOT}/containers/digests.lock"

die() { printf '%s: %s\n' "${0##*/}" "$*" >&2; exit 1; }
warn() { printf '%s: warning: %s\n' "${0##*/}" "$*" >&2; }
say() { printf '==> %s\n' "$*"; }

TARGET=""
TAG=""
VERSION=""
PUSH=0
SKIP_SIF=0
BUILD_ARGS=()
SSH_ARGS=()

[[ $# -gt 0 ]] || { awk 'NR==1{next} /^#/{sub(/^# ?/,"");print;next}{exit}' "$0"; exit 1; }

case "$1" in
    dev|release) TARGET="$1"; shift ;;
    -h|--help)   awk 'NR==1{next} /^#/{sub(/^# ?/,"");print;next}{exit}' "$0"; exit 0 ;;
    *) die "first argument must be 'dev' or 'release' (got '$1')" ;;
esac

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tag)      [[ $# -ge 2 ]] || die "--tag needs a value";     TAG="$2";     shift 2 ;;
        --version)  [[ $# -ge 2 ]] || die "--version needs a value"; VERSION="$2"; shift 2 ;;
        --push)     PUSH=1; shift ;;
        --skip-sif) SKIP_SIF=1; shift ;;
        --no-cache) BUILD_ARGS+=(--no-cache); shift ;;
        --ssh)      SSH_ARGS=(--ssh default); shift ;;
        -h|--help)  awk 'NR==1{next} /^#/{sub(/^# ?/,"");print;next}{exit}' "$0"; exit 0 ;;
        *) die "unknown arg: $1 (see --help)" ;;
    esac
done

command -v docker >/dev/null 2>&1 || die "docker is not on PATH."
[[ -f "${REPO_ROOT}/pixi.toml" ]] || die "no pixi.toml at ${REPO_ROOT} — the build context is the repo root."
[[ -f "${REPO_ROOT}/pixi.lock" ]] \
    || die "no pixi.lock at ${REPO_ROOT}. The Dockerfile installs with --locked and will fail without one. Run 'pixi install' first and commit the lock."

VCS_REF="$(git -C "${REPO_ROOT}" rev-parse --short HEAD 2>/dev/null || echo unknown)"
if [[ -n "$(git -C "${REPO_ROOT}" status --porcelain 2>/dev/null)" ]]; then
    VCS_REF="dirty-${VCS_REF}"
fi

# Emulation is the default failure mode on this laptop, not an edge case: QEMU runs any x86 binary
# the BUILD itself executes — conda post-link scripts, aligner self-tests, anything that compiles
# and then runs a check. Those images are CI-only. This warning is what stops someone spending an
# hour on a QEMU segfault that a native runner would never produce.
HOST_ARCH="$(docker version --format '{{.Server.Arch}}' 2>/dev/null || echo unknown)"
if [[ "${HOST_ARCH}" != "amd64" ]]; then
    warn "host is ${HOST_ARCH}; building ${PLATFORM} under QEMU emulation. Expect 3-10x slower."
    warn "if this image's build RUNS compiled x86 binaries (conda post-link scripts, tool test"
    warn "suites) it may segfault here and succeed in CI. Such images are CI-only — do not chase it."
fi

# The one thing that must never be wrong. An arm64 image pushed to either registry looks fine
# everywhere until it reaches a compute node.
assert_amd64() {
    local ref="$1" got
    got="$(docker image inspect --format '{{.Os}}/{{.Architecture}}' "${ref}" 2>/dev/null)" \
        || die "cannot inspect ${ref} — refusing to push an image whose architecture is unknown."
    [[ "${got}" == "${PLATFORM}" ]] \
        || die "${ref} is ${got}, not ${PLATFORM}. Refusing to push: it would fail to exec on Kebnekaise."
}

# bash and ps are Nextflow's hard container requirements, so a missing one is a build failure here
# rather than a cryptic process-metrics error inside a pipeline run three weeks later.
smoke_test() {
    local ref="$1"
    say "smoke test: bash, ps"
    docker run --rm --platform "${PLATFORM}" --entrypoint /bin/bash "${ref}" -lc 'bash --version | head -1' \
        || die "bash is missing or broken in ${ref}."
    docker run --rm --platform "${PLATFORM}" --entrypoint /bin/bash "${ref}" -lc 'ps --version' \
        || die "ps is missing in ${ref} — install procps (Nextflow requires it)."
    # Two numbers because they differ by ~3x and people quote the wrong one. `image inspect .Size`
    # is roughly what a registry pull transfers; `image ls` is what the unpacked image occupies,
    # which is also roughly the .sif size on the cluster and what fills $APPTAINER_CACHEDIR.
    # `docker image ls` stores Docker Hub images under the short name, so the docker.io/ prefix has
    # to come off or the reference filter silently matches nothing and prints an empty size.
    say "size: $(docker image inspect --format '{{.Size}}' "${ref}" | awk '{printf "%.0f MB compressed", $1/1048576}')" \
        "/ $(docker image ls --format '{{.Size}}' "${ref#docker.io/}" | head -1) unpacked"
}

build_image() {
    local ref="$1" mode="$2"
    say "building ${ref} for ${PLATFORM}"
    docker buildx build \
        --platform "${PLATFORM}" \
        -f "${DOCKERFILE}" \
        -t "${ref}" \
        --build-arg "VCS_REF=${VCS_REF}" \
        --build-arg "IMAGE_VERSION=${VERSION:-${TAG:-dev}}" \
        "${SSH_ARGS[@]+"${SSH_ARGS[@]}"}" \
        "${BUILD_ARGS[@]+"${BUILD_ARGS[@]}"}" \
        "${mode}" \
        "${REPO_ROOT}"
}

if [[ "${TARGET}" == "dev" ]]; then
    TAG="${TAG:-${VCS_REF}}"
    REF="${DEV_REGISTRY}/${IMAGE_NAME}:${TAG}"

    ( cd "${REPO_ROOT}" && build_image "${REF}" --load )
    assert_amd64 "${REF}"
    smoke_test "${REF}"

    if [[ "${PUSH}" -eq 1 ]]; then
        say "pushing ${REF}"
        docker push "${REF}"
        cat <<EOF

Pushed a DEV image. It is not reproducible and its digest does not belong in digests.lock.
On the HPC2N login node:

    apptainer pull ${IMAGE_NAME}_${TAG}.sif docker://${REF}

EOF
    else
        say "built locally. Add --push to publish to ${DEV_REGISTRY}."
    fi
    exit 0
fi

[[ -n "${VERSION}" ]] || die "release needs --version (e.g. --version 0.1.0)."
[[ "${VERSION}" =~ ^[0-9]+\.[0-9]+\.[0-9]+([-+][0-9A-Za-z.]+)?$ ]] \
    || die "--version '${VERSION}' is not semver; digests.lock entries are keyed on it."
[[ "${VCS_REF}" != dirty-* ]] || die "working tree is dirty — a release digest must be reproducible from a commit."

REF="${GHCR_REPO}:${VERSION}"

# --load first, so the architecture gate and the smoke test run BEFORE anything reaches a registry.
# The second build is a cache hit.
( cd "${REPO_ROOT}" && build_image "${REF}" --load )
assert_amd64 "${REF}"
smoke_test "${REF}"

say "pushing ${REF}"
( cd "${REPO_ROOT}" && build_image "${REF}" --push )
OCI_DIGEST="$(docker buildx imagetools inspect "${REF}" --format '{{.Manifest.Digest}}')"
[[ "${OCI_DIGEST}" == sha256:* ]] || die "could not read the pushed manifest digest for ${REF}."

SIF_DIGEST="skipped"
if [[ "${SKIP_SIF}" -eq 1 ]]; then
    warn "--skip-sif: no ORAS artifact. digests.lock will be incomplete and the ADR-7 byte-identity guarantee does not apply to this release."
else
    command -v apptainer >/dev/null 2>&1 \
        || die "apptainer is not on PATH. There is none on macOS (D-003) — run the release target in CI or on the HPC2N login node, or pass --skip-sif to publish the OCI image alone."
    SIF="${REPO_ROOT}/containers/${IMAGE_NAME}_${VERSION}.sif"
    # Converting from the DIGEST, not the tag, is what makes the .sif correspond to exactly the
    # bytes just pushed even if the tag moves a minute later.
    say "converting ${GHCR_REPO}@${OCI_DIGEST} -> ${SIF##*/}"
    apptainer build --reproducible --force "${SIF}" "docker://${GHCR_REPO}@${OCI_DIGEST}"
    say "pushing oras://${SIF_REPO}:${VERSION}"
    apptainer push "${SIF}" "oras://${SIF_REPO}:${VERSION}"
    # The checksum of the .sif FILE, not of a registry manifest. That is the number the ADR-7
    # acceptance test compares: two `apptainer pull oras://` must reproduce exactly this.
    if command -v sha256sum >/dev/null 2>&1; then
        SIF_DIGEST="sha256:$(sha256sum "${SIF}" | awk '{print $1}')"
    else
        SIF_DIGEST="sha256:$(shasum -a 256 "${SIF}" | awk '{print $1}')"
    fi
fi

{
    printf '# Published container digests for cancer_transcriptomics_teaching. Generated by containers/build.sh.\n'
    printf '# GHCR and ORAS only (D-004): a docker.io/ dev tag must never appear here.\n'
    printf '# Pin by digest, never by tag, in anything archived with a manuscript.\n'
    printf '# digest: for oci, the registry manifest digest; for sif, the sha256 of the .sif file,\n'
    printf '#         which is what two `apptainer pull oras://` must reproduce byte for byte.\n'
    printf '#\n'
    printf '# kind\tversion\treference\tdigest\tbuilt_utc\tvcs_ref\n'
    printf 'oci\t%s\t%s\t%s\t%s\t%s\n' \
        "${VERSION}" "${GHCR_REPO}" "${OCI_DIGEST}" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${VCS_REF}"
    printf 'sif\t%s\t%s\t%s\t%s\t%s\n' \
        "${VERSION}" "oras://${SIF_REPO}" "${SIF_DIGEST}" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${VCS_REF}"
} > "${LOCKFILE}"

cat <<EOF

Wrote ${LOCKFILE}. Commit it — it is the record of what this version actually is.

    apptainer pull oras://${SIF_REPO}:${VERSION}
    apptainer pull docker://${GHCR_REPO}@${OCI_DIGEST}

EOF
