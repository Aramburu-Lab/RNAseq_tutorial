#!/usr/bin/env bash
#SBATCH -J ctx_rnaseq
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=12G
#SBATCH -t 1-00:00:00
#SBATCH -o slurm-%x-%j.out
# /// codebase-meta
# exempt: runs on Kebnekaise into the student's run folder, never into outputs/
# ///
# One day: Step-1's download of 210 runs took 1 h 39 min and its count, the Salmon index included,
# 1 h 52 min (2026-09-30). A 7-day head job waited on ReqNodeNotAvail until its limit was lowered to
# 1 day, which cleared it at once (2026-10-01).
# 12 GB: the head job peaked at 6.16, 6.61 and 7.49 GiB of its 8 on three breast studies (P4.5).
# Gene counts with nf-core/rnaseq from the FASTQs fetchngs.sbatch downloaded: STAR and Salmon
# where the site aligns (SITE_RNASEQ_ALIGN, D-113), Salmon alone always, so M1 to M5 read one
# matrix on every computer. It
# starts only once fetchngs.sbatch has written nfcore_fetchngs/results.complete in the run folder,
# and only with references equal to the publisher's MD5SUMS where one sits beside them. Another
# genome: its FASTA and GTF in place of GENCODE's, and false for a GTF not from GENCODE.
#
# Usage:
#   mkdir -p "$RUN/refs" && cd "$RUN/refs"
#   B=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50
#   curl -fL --remote-name-all "$B/MD5SUMS" "$B/GRCh38.primary_assembly.genome.fa.gz" "$B/gencode.v50.primary_assembly.annotation.gtf.gz"
#   cd "$RUN"
#   DRY_RUN=1 bash "$CTX_REPO/scripts/bash/rnaseq.sbatch" "$RUN" refs/GRCh38.primary_assembly.genome.fa.gz refs/gencode.v50.primary_assembly.annotation.gtf.gz true
#   sbatch -A <compute project> "$CTX_REPO/scripts/bash/rnaseq.sbatch" "$RUN" refs/GRCh38.primary_assembly.genome.fa.gz refs/gencode.v50.primary_assembly.annotation.gtf.gz true
#   # keep nfcore_rnaseq/work until ctx gather on the laptop has read the copied results
#   # on a laptop, Docker running: the references where ctx setup finds the GTF, then the count
#   REFS="$HOME/ctx_<run name, no spaces>/resources/gencode"
#   mkdir -p "$REFS" && cd "$REFS"
#   cd "$RUN" && CTX_SITE=laptop bash "$CTX_REPO/scripts/bash/rnaseq.sbatch" "$RUN" "$REFS/GRCh38.primary_assembly.genome.fa.gz" "$REFS/gencode.v50.primary_assembly.annotation.gtf.gz" true
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
# From CTX_REPO, never from $0: sbatch runs a copy of this file from its spool folder. The site
# file says what this computer changes; read before the dry run, which prints its profile.
CTX_SITE="${CTX_SITE:-kebnekaise}"
[[ "$CTX_SITE" =~ ^[a-z0-9_]+$ ]] || die "CTX_SITE is a site name such as kebnekaise or laptop, not $CTX_SITE"
SITE_FILE="${CTX_REPO:?export CTX_REPO=<your clone of the course repository>}/scripts/bash/site_$CTX_SITE.sh"
[[ -f "$SITE_FILE" ]] || die "no $SITE_FILE; is CTX_REPO your clone of the course repository, and CTX_SITE a site it has?"
source "$SITE_FILE"
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
"$home"/*) [[ "${SITE_HOME_ALLOWED:-0}" == 1 ]] ||
    die "the run folder $RUN is in your home folder (25 GB); use a folder under CTX_STORE" ;;
esac
[[ -s "$FASTA" ]] || die "no genome FASTA, or an empty one: $FASTA"
[[ -s "$GTF" ]] || die "no GTF, or an empty one: $GTF"
# A download still running, or one that failed its end check, would be quantified as it stands.
MARKER="$RUN/nfcore_fetchngs/results.complete"
[[ -f "$MARKER" ]] || die "no $MARKER: fetchngs.sbatch has not finished this download"
SHEET="$RUN/nfcore_fetchngs/results/samplesheet/samplesheet.csv"
[[ -s "$SHEET" ]] || die "no samplesheet at $SHEET"
# Asked of the site before anything is written: on a laptop, what Docker has.
executor=$(site_executor) || die "the site's executor is unknown; see the line above"
limits=$(site_rnaseq_config) || die "the site's limits are unknown; see the line above"
LAUNCH="$RUN/nfcore_rnaseq"
mkdir -p "$LAUNCH"
[[ "${SITE_RNASEQ_ALIGN:-}" == true || "${SITE_RNASEQ_ALIGN:-}" == false ]] || die "the site file sets no SITE_RNASEQ_ALIGN"
if [[ "$SITE_RNASEQ_ALIGN" == true ]]; then route="aligner: star_salmon"; else route="skip_alignment: true"; fi

# publish_dir_mode is stated: under NXF_SYNTAX_PARSER=v2 the config below cannot see a param only
# the pipeline sets, and its reports were published as symlinks into work/ (D-025).
cat > "$LAUNCH/params.yml" <<EOF
input: $LAUNCH/samplesheet.csv
outdir: $LAUNCH/results
fasta: $FASTA
gtf: $GTF
gencode: $GENCODE
pseudo_aligner: salmon
$route
save_reference: $SITE_SAVE_REFERENCE
publish_dir_mode: copy
EOF
# The site's executor and limits. Then the Trim Galore JSON reports, the only ones MultiQC 1.35
# reads for Trim Galore 2.x, which ctx qc needs (D-025).
printf '%s\n%s\n' "$executor" "$limits" > "$LAUNCH/site.config"
cat >> "$LAUNCH/site.config" <<'EOF'
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
cmd=(nextflow run nf-core/rnaseq -r 3.27.0 -profile "$SITE_RNASEQ_PROFILE" -params-file "$LAUNCH/params.yml" -c "$LAUNCH/site.config" -resume)
printf 'the job runs, in %s:\n%s\n' "$LAUNCH" "${cmd[*]}"
if [[ "$DRY_RUN" == 1 ]]; then exit 0; fi

# Activated first, so the checks below run on the environment's Python and md5sum.
site_activate || die "the fetch environment did not activate; see the line above"
# The publisher's MD5SUMS beside a reference, where there is one, compared hash to hash: macOS's
# md5sum -c passes an empty list, so a list that names neither file would pass it.
for ref in "$FASTA" "$GTF"; do
    sums="${ref%/*}/MD5SUMS"
    if [[ ! -f "$sums" ]]; then
        echo "${ref##*/}: no MD5SUMS beside it, so not checked against its publisher"
        continue
    fi
    want=$(awk -v n="${ref##*/}" '$2 == n || $2 == "*" n {print $1}' "$sums")
    [[ "$want" =~ ^[0-9a-f]{32}$ ]] || die "$sums has no single line for ${ref##*/}; download both again"
    got=$(md5sum <"$ref" | awk '{print $1}')
    [[ "$got" == "$want" ]] || die "${ref##*/}: MD5 $got, not the $want in $sums; download it again"
    echo "${ref##*/}: MD5 equal to $sums"
done
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
if [[ "$SITE_RNASEQ_ALIGN" == true ]]; then
    A="$LAUNCH/results/star_salmon/salmon.merged.gene_counts_length_scaled.tsv"
    [[ -s "$A" ]] || die "no count matrix from the alignment at $A"
    echo "alignment count matrix: $(awk 'END{print NR - 1}' "$A") genes, $(head -n 1 "$A" | awk -F'\t' '{print NF - 2}') samples"
fi
