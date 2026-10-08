#!/usr/bin/env bash
#SBATCH -J ctx_subsample
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 0-04:00:00
#SBATCH -o slurm-%x-%j.out
# /// codebase-meta
# exempt: runs on Kebnekaise into the student's run folder, never into outputs/
# ///
# A teaching subset: a copy of one fetchngs download with every FASTQ cut to READS reads by
# seqtk, fetchngs's file names kept and the md5 files rewritten, so ctx gather reads it as a
# download. The download itself is left as it was. The run folder is one of its own, whose
# fetchngs.sbatch download holds only the chosen runs.
#
# Usage:
#   sbatch -A <compute project> "$CTX_REPO/scripts/bash/subsample.sbatch" "$RUN"
#   # after the job: its time and peak memory
#   sacct -j <job id> --format=JobID,Elapsed,MaxRSS,State
set -euo pipefail
# An unmatched pattern is then an empty list, never a file named *.fastq.gz.
shopt -s nullglob

RUN="${1:?usage: subsample.sbatch <run folder>}"
READS="${READS:-1000000}"
SEED="${SEED:-42}"
CTX_STORE="${CTX_STORE:?export CTX_STORE=<your folder on the storage project>}"
# Kebnekaise only: Apptainer, and a run folder never in home.
SITE_HOME_ALLOWED=0

die() {
    printf 'subsample.sbatch: %s\n' "$*" >&2
    exit 1
}
[[ "$READS" =~ ^[1-9][0-9]*$ && "$SEED" =~ ^[0-9]+$ ]] || die "READS and SEED are whole numbers"
[[ "$RUN" == /* ]] || RUN="$PWD/$RUN"
# Home is 25 GB, and the run folder takes every FASTQ. Compared canonically, from the nearest
# folder that exists: the run folder may not exist yet, and either side may be a symlink.
here=$RUN
# dirname ends at . or /, both folders, so the walk ends for any path.
while [[ ! -d "$here" ]]; do here=$(dirname -- "$here"); done
here=$(CDPATH='' cd -- "$here" && pwd -P) || die "cannot enter the folder that holds $RUN"
home=$(CDPATH='' cd -- "${HOME:-/nonexistent}" 2>/dev/null && pwd -P) || home=${HOME:-/nonexistent}
case "$here/" in
"$home"/*) [[ "${SITE_HOME_ALLOWED:-0}" == 1 ]] ||
    die "the run folder $RUN is in your home folder (25 GB); use a folder under CTX_STORE" ;;
esac
[[ "$CTX_STORE" == /* ]] || CTX_STORE="$PWD/$CTX_STORE"
R="$RUN/nfcore_fetchngs/results"
O="$RUN/subset"
# A download still running leaves no marker; subsampling it wrote a hollow outdir once.
MARKER="$RUN/nfcore_fetchngs/results.complete"
[[ -f "$MARKER" ]] || die "no $MARKER: fetchngs.sbatch has not finished this download"
fastq=("$R"/fastq/*.fastq.gz)
[[ ${#fastq[@]} -gt 0 ]] || die "no FASTQ in $R/fastq; a metadata-only run has none"
[[ ! -e "$O" ]] || die "$O exists; move it aside first"

IMAGE="$CTX_STORE/images/seqtk_1.4.sif"
if [[ ! -s "$IMAGE" ]]; then
    mkdir -p "${IMAGE%/*}"
    # To .part first: an interrupted download never looks like an image.
    wget -nv -O "$IMAGE.part" https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_1 ||
        die "the seqtk image did not download; submit again"
    mv "$IMAGE.part" "$IMAGE"
fi
# seqtk with no command prints its usage and exits 1.
{ apptainer exec "$IMAGE" seqtk 2>&1 || true; } | grep -q '^Version: 1\.4-' ||
    die "$IMAGE is not seqtk 1.4; move it aside and submit again"

rsync -a --exclude=/fastq/ "$R/" "$O/"
mkdir -p "$O/fastq/md5"
for f in "${fastq[@]}"; do
    n=${f##*/}
    # One seed for every file, so the two files of a pair keep the same reads. The FASTQ goes in
    # on stdin, so the container needs no bind of the run folder.
    apptainer exec "$IMAGE" seqtk sample -s"$SEED" - "$READS" < "$f" | gzip -n > "$O/fastq/$n"
    (cd "$O/fastq" && md5sum "$n" > "md5/$n.md5")
done
{
    printf '# Teaching subset\n\n'
    # Relative to the run folder: this README travels, the storage path stays behind.
    printf 'The fetchngs outdir nfcore_fetchngs/results, each FASTQ replaced by `seqtk sample -s%s - %s | gzip -n`' "$SEED" "$READS"
    printf ' (seqtk 1.4) and fastq/md5/ rewritten. read_count, base_count, fastq_bytes and'
    printf ' fastq_md5 in samplesheet/ and metadata/ describe the ORIGINAL runs.\n\n'
    printf '| file | reads | bytes | md5 |\n|---|---|---|---|\n'
    for f in "$O"/fastq/*.fastq.gz; do
        n=${f##*/}
        printf '| %s | %s | %s | %s |\n' "$n" "$(gzip -dc "$f" | awk 'END{print NR / 4}')" \
            "$(($(wc -c < "$f")))" "$(cut -d' ' -f1 "$O/fastq/md5/$n.md5")"
    done
} > "$O/README.md"
echo "subset: ${#fastq[@]} FASTQ files of at most $READS reads each in $O"
