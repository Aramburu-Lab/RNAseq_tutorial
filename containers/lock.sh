#!/usr/bin/env bash
# /// codebase-meta
# exempt: >-
#   Build tooling, not a pipeline step. Sourced by containers/build.sh and by
#   .github/workflows/containers.yml, the two writers of containers/digests.lock.
# ///
# One definition for both writers of the lock (TODO.md F.2d). A row is the artifact (kind,
# version, reference) and its digest. A writer that rewrote the file wholesale dropped what it
# did not build: build.sh release, which builds the analysis image alone, dropped the figures
# image's row the Reproduce panels pin by.

# lock_write LOCK ROW...: each ROW replaces its artifact's row; every other line stays as it was.
lock_write() {
    local lock="$1" tmp
    shift
    tmp="${lock}.tmp.$$"
    {
        if [[ -f "${lock}" ]]; then
            # A comment's first fields are never a row's key, so comments stay too.
            awk -F'\t' '
                NR == FNR { new[$1 FS $2 FS $3] = 1; next }
                !(($1 FS $2 FS $3) in new)
            ' <(printf '%s\n' "$@") "${lock}"
        else
            printf '# Published container digests. Written by containers/lock.sh.\n'
            printf '# kind\tversion\treference\tdigest\tbuilt_utc\tvcs_ref\n'
        fi
        printf '%s\n' "$@"
    } > "${tmp}"
    mv "${tmp}" "${lock}"
}
