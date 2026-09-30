#!/usr/bin/env bash
#SBATCH -J ctx_fetchngs
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 3-00:00:00
#SBATCH -o slurm-%x-%j.out
# /// codebase-meta
# exempt: runs on Kebnekaise into the student's run folder, never into outputs/
# ///
# One accession list through nf-core/fetchngs into a fresh outdir, laid out as ctx gather reads
# it. The account is never written here: sbatch -A names the compute project. Run the Usage
# lines on a login node with CTX_REPO and CTX_STORE exported (kebnekaise_env.sh).
#
# Usage:
#   RUN="$CTX_STORE/<run name, no spaces>"
#   mkdir -p "$RUN" && cd "$RUN"
#   printf '%s\n' <accession> <accession> > ids.txt
#   DRY_RUN=1 bash "$CTX_REPO/scripts/bash/fetchngs.sbatch" "$RUN" ids.txt
#   SKIP_DOWNLOAD=1 sbatch -A <compute project> "$CTX_REPO/scripts/bash/fetchngs.sbatch" "$RUN" ids.txt
#   # once that job has ended: the size of the download, read before it is submitted
#   grep 'bytes of FASTQ at ENA' slurm-ctx_fetchngs-<job id>.out
#   sbatch -A <compute project> "$CTX_REPO/scripts/bash/fetchngs.sbatch" "$RUN" ids.txt
#   # a failed download: the dry run prints the mv that moves its results aside to a free name;
#   # after that mv, the same sbatch line resumes the rest
#   DRY_RUN=1 bash "$CTX_REPO/scripts/bash/fetchngs.sbatch" "$RUN" ids.txt
#   sbatch -A <compute project> "$CTX_REPO/scripts/bash/fetchngs.sbatch" "$RUN" ids.txt
set -euo pipefail

RUN="${1:?usage: fetchngs.sbatch <run folder> <accession list>}"
IDS="${2:?usage: fetchngs.sbatch <run folder> <accession list>}"
SKIP_DOWNLOAD="${SKIP_DOWNLOAD:-0}"
MAX_FORKS="${MAX_FORKS:-5}"
DRY_RUN="${DRY_RUN:-0}"

die() {
    printf 'fetchngs.sbatch: %s\n' "$*" >&2
    exit 1
}
# From CTX_REPO, never from $0: sbatch runs a copy of this file from its spool folder.
ENV_FILE="${CTX_REPO:?export CTX_REPO=<your clone of the course repository>}/scripts/bash/kebnekaise_env.sh"
[[ -f "$ENV_FILE" ]] || die "no $ENV_FILE; is CTX_REPO your clone of the course repository?"
[[ "$SKIP_DOWNLOAD" == [01] && "$DRY_RUN" == [01] ]] || die "SKIP_DOWNLOAD and DRY_RUN take 0 or 1"
[[ "$MAX_FORKS" =~ ^[1-9][0-9]*$ ]] || die "MAX_FORKS is a whole number, not $MAX_FORKS"
[[ "$RUN" == /* ]] || RUN="$PWD/$RUN"
# fetchngs's schema takes an --input path only without whitespace (^\S+\.(csv|tsv|txt)$), and
# the job passes an absolute one, so the folder a relative run folder was given from counts.
[[ "$RUN" != *[[:space:]]* ]] || die "the run folder has a space, which fetchngs rejects: $RUN"
# Home is 25 GB, and the run folder takes every FASTQ. Compared canonically, from the nearest
# folder that exists: the run folder may not exist yet, and either side may be a symlink.
here=$RUN
# dirname ends at . or /, both folders, so the walk ends for any path.
while [[ ! -d "$here" ]]; do here=$(dirname -- "$here"); done
here=$(CDPATH='' cd -- "$here" && pwd -P) || die "cannot enter the folder that holds $RUN"
home=$(CDPATH='' cd -- "${HOME:-/nonexistent}" 2>/dev/null && pwd -P) || home=${HOME:-/nonexistent}
case "$here/" in
"$home"/*) die "the run folder $RUN is in your home folder (25 GB); use a folder under CTX_STORE" ;;
esac
[[ -s "$IDS" ]] || die "no accession list, or an empty one: $IDS"
# A metadata-only run keeps its own outdir: one params file per outdir is what ctx gather reads.
LAUNCH="$RUN/nfcore_fetchngs"
if [[ "$SKIP_DOWNLOAD" == 1 ]]; then LAUNCH="$RUN/nfcore_fetchngs_meta"; fi
OUT="$LAUNCH/results"
if [[ -e "$OUT" ]]; then
    # A name already taken would put results inside it, or stop mv.
    try=1
    while [[ -e "${OUT}_try$try" ]]; do try=$((try + 1)); done
    die "$OUT exists; move it aside (mv $OUT ${OUT}_try$try) and submit again"
fi
mkdir -p "$LAUNCH"

# One accession per line, as fetchngs's isSraId reads them: one line outside its pattern stops
# the whole run, so a byte-order mark, a Windows line end, a blank line or a comment is removed
# here and anything else is refused with its line number.
ID_RE='^(((SR|ER|DR)[APRSX])|(SAM(N|EA|D))|(PRJ(NA|EB|DB))|(GS[EM]))[0-9]+$'
SEPARATOR='[[:space:],;]'
n=0
bad=0
: > "$LAUNCH/ids.csv.part"
while IFS= read -r line || [[ -n "$line" ]]; do
    n=$((n + 1))
    if [[ $n -eq 1 ]]; then line=${line#$'\xef\xbb\xbf'}; fi
    # [:space:] takes the carriage return of a Windows line end.
    line=${line#"${line%%[![:space:]]*}"}
    line=${line%"${line##*[![:space:]]}"}
    if [[ -z "$line" || "$line" == "#"* ]]; then continue; fi
    if [[ "$line" =~ $SEPARATOR ]]; then
        echo "line $n: one accession per line, one column: $line" >&2
        bad=1
    elif [[ "$line" == E-MTAB-* ]]; then
        echo "line $n: $line is an ArrayExpress study, which fetchngs rejects;" \
            "give its ENA study (PRJEB<n> or ERP<n>) from the study's BioStudies page" >&2
        bad=1
    elif [[ ! "$line" =~ $ID_RE ]]; then
        echo "line $n: $line is not an SRA, ENA, DDBJ or GEO accession fetchngs takes" >&2
        bad=1
    else
        printf '%s\n' "$line" >> "$LAUNCH/ids.csv.part"
    fi
done < "$IDS"
[[ $bad -eq 0 ]] || die "fix the lines above in $IDS"
[[ -s "$LAUNCH/ids.csv.part" ]] || die "no accession left in $IDS"
mv "$LAUNCH/ids.csv.part" "$LAUNCH/ids.csv"

# D-035: wget retries a refused connection, with a growing pause, MAX_FORKS downloads at once.
cat > "$LAUNCH/nextflow.config" <<CONF
executor {
    queueSize = 50
}
process {
    withName: 'SRA_FASTQ_FTP' {
        ext.args      = '-t 5 -nv -c -T 60 --retry-connrefused --waitretry 60'
        time          = { 30.min * task.attempt }
        maxForks      = $MAX_FORKS
        maxRetries    = 6
        errorStrategy = { sleep(Math.pow(2, task.attempt) * 30000 as long); return 'retry' }
    }
}
CONF
# The input list sits inside the outdir, where ctx gather looks for it. -resume finds nothing on
# a first launch; after results/ is moved aside it reuses every finished download.
cmd=(nextflow run nf-core/fetchngs -r 1.12.0 -profile singularity --input "$OUT/ids.csv" --outdir "$OUT"
    --nf_core_pipeline rnaseq --monochrome_logs true --publish_dir_mode link -resume)
if [[ "$SKIP_DOWNLOAD" == 1 ]]; then cmd+=(--skip_fastq_download); fi
printf '%s accessions in %s; the job runs, in %s:\n%s\n' \
    "$(awk 'END{print NR}' "$LAUNCH/ids.csv")" "$LAUNCH/ids.csv" "$LAUNCH" "${cmd[*]}"
if [[ "$DRY_RUN" == 1 ]]; then exit 0; fi

source "$ENV_FILE" || die "the fetch environment did not activate; see the line above"
rm -f "$LAUNCH/results.complete"
mkdir -p "$OUT"
cp "$LAUNCH/ids.csv" "$OUT/ids.csv"
cd "$LAUNCH"
# fetchngs 1.12.0's config fails Nextflow 26's strict parser; v1 reads it.
NXF_SYNTAX_PARSER=v1 "${cmd[@]}"

# fetchngs's guarantee is one metadata file per input line; a download must also hold every
# FASTQ its samplesheet names. The sheet is read by column name, never by position.
python3 - "$OUT" "$SKIP_DOWNLOAD" <<'CHECK' || die "the outdir is incomplete; no marker written"
import csv
import sys
from pathlib import Path

out, download = Path(sys.argv[1]), sys.argv[2] == "0"
ids = (out / "ids.csv").read_text().split()
problems = [f"no metadata/{i}.runinfo_ftp.tsv" for i in ids
            if not (out / "metadata" / f"{i}.runinfo_ftp.tsv").is_file()]
sheet = out / "samplesheet" / "samplesheet.csv"
rows = list(csv.DictReader(sheet.open(newline=""))) if sheet.is_file() else []
problems += [] if rows else ["no rows in samplesheet/samplesheet.csv"]
if download:
    names = [Path(r[k]).name for r in rows for k in ("fastq_1", "fastq_2") if r[k]]
    problems += [f"no fastq/{n}" for n in names if not (out / "fastq" / n).is_file()]
size = sum(int(b) for r in rows for b in r["fastq_bytes"].split(";") if b)
print(f"{len(ids)} accessions; {len(rows)} runs; {size:,} bytes of FASTQ at ENA")
if problems:
    sys.exit("\n".join(problems))
CHECK
if [[ "$SKIP_DOWNLOAD" == 0 ]]; then
    echo "FASTQ files not hard-linked to work/: $(($(find "$OUT/fastq" -name '*.fastq.gz' -links 1 | wc -l)))"
fi
# rnaseq.sbatch and subsample.sbatch start only from a download this check passed.
date > "$LAUNCH/results.complete"
echo "complete: $LAUNCH/results.complete"
