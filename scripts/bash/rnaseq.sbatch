#!/usr/bin/env bash
#SBATCH -J ctx_rnaseq
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 7-00:00:00
#SBATCH -o slurm-%x-%j.out
# /// codebase-meta
# exempt: runs on Kebnekaise into the student's run folder, never into outputs/
# ///
# Gene counts with nf-core/rnaseq, Salmon only, from the FASTQs fetchngs.sbatch downloaded. It
# starts only once fetchngs.sbatch has written nfcore_fetchngs/results.complete in the run folder.
#
# Usage:
#   mkdir -p "$RUN/refs" && cd "$RUN/refs"
#   wget -nv <genome FASTA URL> <GTF URL>
#   md5sum <genome FASTA> <GTF>
#   # each sum equal to the publisher's MD5 list, or download that file again
#   cd "$RUN"
#   DRY_RUN=1 bash "$CTX_REPO/scripts/bash/rnaseq.sbatch" "$RUN" refs/<genome FASTA> refs/<GTF> <true for a GENCODE GTF, else false>
#   sbatch -A <compute project> "$CTX_REPO/scripts/bash/rnaseq.sbatch" "$RUN" refs/<genome FASTA> refs/<GTF> <true or false>
#   # keep nfcore_rnaseq/work until ctx gather on the laptop has read the copied results
set -euo pipefail

RUN="${1:?usage: rnaseq.sbatch <run folder> <genome FASTA> <GTF> <true|false, a GENCODE GTF>}"
FASTA="${2:?usage: rnaseq.sbatch <run folder> <genome FASTA> <GTF> <true|false, a GENCODE GTF>}"
GTF="${3:?usage: rnaseq.sbatch <run folder> <genome FASTA> <GTF> <true|false, a GENCODE GTF>}"
GENCODE="${4:?usage: rnaseq.sbatch <run folder> <genome FASTA> <GTF> <true|false, a GENCODE GTF>}"
DRY_RUN="${DRY_RUN:-0}"

die() {
    printf 'rnaseq.sbatch: %s\n' "$*" >&2
    exit 1
}
# From CTX_REPO, never from $0: sbatch runs a copy of this file from its spool folder.
ENV_FILE="${CTX_REPO:?export CTX_REPO=<your clone of the course repository>}/scripts/bash/kebnekaise_env.sh"
[[ -f "$ENV_FILE" ]] || die "no $ENV_FILE; is CTX_REPO your clone of the course repository?"
[[ "$DRY_RUN" == [01] ]] || die "DRY_RUN takes 0 or 1"
# No default: nf-core/rnaseq builds the Salmon index differently for a GENCODE GTF, so the
# student states which kind of GTF this is.
[[ "$GENCODE" == true || "$GENCODE" == false ]] || die "GENCODE is true or false, not $GENCODE"
[[ "$RUN" == /* ]] || RUN="$PWD/$RUN"
[[ "$FASTA" == /* ]] || FASTA="$PWD/$FASTA"
[[ "$GTF" == /* ]] || GTF="$PWD/$GTF"
# nf-core/rnaseq 3.27.0's schema takes input ^\S+\.csv$, fasta ^\S+\.fn?a(sta)?(\.gz)?$ and gtf
# ^\S+\.gtf(\.gz)?$, and says so only once the job has left the queue.
for p in "$RUN" "$FASTA" "$GTF"; do
    [[ "$p" != *[[:space:]]* ]] || die "$p has a space, which nf-core/rnaseq rejects"
done
FASTA_RE='\.fn?a(sta)?(\.gz)?$'
GTF_RE='\.gtf(\.gz)?$'
[[ "$FASTA" =~ $FASTA_RE ]] || die "nf-core/rnaseq takes a genome FASTA only if its name matches $FASTA_RE: $FASTA"
[[ "$GTF" =~ $GTF_RE ]] || die "nf-core/rnaseq takes a GTF only if its name matches $GTF_RE: $GTF"
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
[[ -s "$FASTA" ]] || die "no genome FASTA, or an empty one: $FASTA"
[[ -s "$GTF" ]] || die "no GTF, or an empty one: $GTF"
# A download still running, or one that failed its end check, would be quantified as it stands.
MARKER="$RUN/nfcore_fetchngs/results.complete"
[[ -f "$MARKER" ]] || die "no $MARKER: fetchngs.sbatch has not finished this download"
SHEET="$RUN/nfcore_fetchngs/results/samplesheet/samplesheet.csv"
[[ -s "$SHEET" ]] || die "no samplesheet at $SHEET"
LAUNCH="$RUN/nfcore_rnaseq"
mkdir -p "$LAUNCH"

# publish_dir_mode is stated: under NXF_SYNTAX_PARSER=v2 the config below cannot see a param only
# the pipeline sets, and its reports were published as symlinks into work/ (D-025).
cat > "$LAUNCH/params.yml" <<EOF
input: $LAUNCH/samplesheet.csv
outdir: $LAUNCH/results
fasta: $FASTA
gtf: $GTF
gencode: $GENCODE
pseudo_aligner: salmon
skip_alignment: true
save_reference: true
publish_dir_mode: copy
EOF
# queueSize 50: HPC2N refused 100 jobs at once. Then the Trim Galore JSON reports, the only ones
# MultiQC 1.35 reads for Trim Galore 2.x, which ctx qc needs (D-025).
cat > "$LAUNCH/kebnekaise.config" <<'EOF'
executor {
    queueSize = 50
}
EOF
cat >> "$LAUNCH/kebnekaise.config" <<'EOF'
process {
    withName: '.*:FASTQ_FASTQC_UMITOOLS_TRIMGALORE:TRIMGALORE' {
        publishDir = [
            [
                path: { "${params.outdir}/${task.ext.publish_prefix}fastqc/trim" },
                mode: params.publish_dir_mode,
                pattern: "*.{html,zip}"
            ],
            [
                path: { params.save_trimmed ? "${params.outdir}/${task.ext.publish_prefix}${params.trimmer}" : params.outdir },
                mode: params.publish_dir_mode,
                pattern: "*.fq.gz",
                saveAs: { filename -> params.save_trimmed ? filename : null }
            ],
            [
                path: { "${params.outdir}/${task.ext.publish_prefix}${params.trimmer}" },
                mode: params.publish_dir_mode,
                pattern: "*_trimming_report.{txt,json}"
            ]
        ]
    }
}
EOF
cmd=(nextflow run nf-core/rnaseq -r 3.27.0 -profile singularity -params-file "$LAUNCH/params.yml" -c "$LAUNCH/kebnekaise.config" -resume)
printf 'the job runs, in %s:\n%s\n' "$LAUNCH" "${cmd[*]}"
if [[ "$DRY_RUN" == 1 ]]; then exit 0; fi

# Sourced first, so the cut below runs on the environment's Python, not the system's.
source "$ENV_FILE" || die "the fetch environment did not activate; see the line above"
python3 - "$SHEET" "$LAUNCH/samplesheet.csv" <<'CUT'
import csv
import sys

keep = ("sample", "fastq_1", "fastq_2", "strandedness")
with open(sys.argv[1], newline="") as fh:
    reader = csv.DictReader(fh)
    missing = [k for k in keep if k not in (reader.fieldnames or [])]
    if missing:
        sys.exit(f"{sys.argv[1]}: no column {missing}; was fetchngs run with --nf_core_pipeline rnaseq?")
    rows = [{k: r[k] for k in keep} for r in reader]
with open(sys.argv[2], "w", newline="") as fh:
    out = csv.DictWriter(fh, fieldnames=keep, lineterminator="\n")
    out.writeheader()
    out.writerows(rows)
print(f"rnaseq samplesheet: {len(rows)} rows, {len({r['sample'] for r in rows})} samples")
CUT
cd "$LAUNCH"
NXF_SYNTAX_PARSER=v2 "${cmd[@]}"

# A link into work/ dangles once work/ is cleaned, and on every other machine (D-025).
links=$(($(find "$LAUNCH/results" -type l | wc -l)))
[[ $links -eq 0 ]] || die "$links published files are links into work/; keep work/ and copy results/ with rsync -L"
M="$LAUNCH/results/salmon/salmon.merged.gene_counts_length_scaled.tsv"
[[ -s "$M" ]] || die "no count matrix at $M"
echo "count matrix: $(awk 'END{print NR - 1}' "$M") genes, $(head -n 1 "$M" | awk -F'\t' '{print NF - 2}') samples"
