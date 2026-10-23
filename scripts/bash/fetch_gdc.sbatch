#!/usr/bin/env bash
#SBATCH -J RNAseq-Tutorial_NA_gdc-download_NA
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 0-06:00:00
#SBATCH -o slurm-%x-%j.out
# /// codebase-meta
# exempt: runs on Kebnekaise or a laptop into a ctx root's data/raw folder, never into outputs/
# ///
# Six hours: one head and neck project's three kinds, 2.65 GB, took 5,965 s on a laptop
# (2026-10-18), where the kidney project's STAR files, 2.6 GB, had taken 665 s (2026-10-01).
# 2 GB: the end check holds one project's clinical record, 10.7 MB for breast, the largest.
# The open-access GDC files of one project, laid out as ctx gather --gdc-dir, --mirna-dir and
# --isoform-dir read them: STAR - Counts in <folder>, with the file metadata, manifest, clinical
# record, annotations, read groups, release and retrieval log, and for a TCGA project the TCGA-CDR;
# miRNA and isoform files in <folder>_mirna and <folder>_isoform, each with its own metadata.
# Every data file lands as .part and is kept only when its MD5 equals the manifest's. Run it again
# to resume: metadata saved by an earlier run is kept, so a rerun stays on its first manifest, and
# a file is skipped once its MD5 matches. It ends with exit 1, naming each problem, while any file
# is missing, partial or mismatched. A folder of another project or workflow is refused.
# DATA=0 saves the metadata alone and prints the size; KINDS picks among star, mirna and isoform;
# JOBS sets the parallel downloads (4); CDR=1 or 0 fetches the CDR or not, whatever the project.
# For tests only: GDC_API, GDC_DATA, GDC_WEB, RETRIES (6).
#
# Usage:
#   # on Kebnekaise, after kebnekaise_env.sh's once-per-login block: the metadata, which prints the size
#   RAW="$CTX_STORE/ctx_<slug>/data/raw"
#   DATA=0 bash "$CTX_REPO/scripts/bash/fetch_gdc.sbatch" <project> "$RAW/<slug>"
#   # its last line is the size at GDC: read it, then download
#   sbatch -A <compute project> "$CTX_REPO/scripts/bash/fetch_gdc.sbatch" <project> "$RAW/<slug>"
#   # once that job has ended: every check passed, or the problems it names
#   tail -n 3 slurm-RNAseq-Tutorial_NA_gdc-download_NA-<job id>.out
#   # on a laptop: the clone, then the metadata, which prints the size
#   export CTX_REPO=<your clone of the course repository> RAW="$HOME/ctx_<slug>/data/raw"
#   git clone <course repository> "$CTX_REPO"
#   pixi install --manifest-path "$CTX_REPO/pixi.toml" -e fetch --frozen
#   CTX_SITE=laptop DATA=0 bash "$CTX_REPO/scripts/bash/fetch_gdc.sbatch" <project> "$RAW/<slug>"
#   # its last line is the size at GDC: read it, then download
#   CTX_SITE=laptop bash "$CTX_REPO/scripts/bash/fetch_gdc.sbatch" <project> "$RAW/<slug>"
set -uo pipefail
# No -e: one failed file must not stop the others; the end check counts every failure.

PROJECT="${1:?usage: fetch_gdc.sbatch <GDC project> <folder>}"
RUN="${2:?usage: fetch_gdc.sbatch <GDC project> <folder>}"
DATA="${DATA:-1}"
KINDS="${KINDS:-star mirna isoform}"
JOBS="${JOBS:-4}"
RETRIES="${RETRIES:-6}"
CDR="${CDR:-auto}"
GDC_API="${GDC_API:-https://api.gdc.cancer.gov}"
GDC_DATA="${GDC_DATA:-$GDC_API/data}"
GDC_WEB="${GDC_WEB:-https://gdc.cancer.gov}"
# The TCGA-CDR (Liu et al. 2018), one pan-cancer workbook, by its UUID in GDC's PanCan manifest.
CDR_ID=1b5f413e-a8d1-4d10-92eb-7c4ae739ed81
PANCAN=PanCan-General_Open_GDC-Manifest_2.txt
# follow_ups.molecular_tests for every project, so no project's clinical record lacks a field
# another has, and a count of zero from it is a count.
EXPAND=demographic,diagnoses,diagnoses.treatments,exposures,follow_ups,follow_ups.molecular_tests
# preservation_method tells FFPE files from frozen ones, which the one-file-per-case rule needs.
FIELDS=file_id,file_name,md5sum,file_size,cases.case_id,cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type,cases.samples.preservation_method,cases.samples.portions.analytes.aliquots.submitter_id
READ_FIELDS=file_id,cases.samples.submitter_id,cases.samples.portions.analytes.aliquots.submitter_id,analysis.metadata.read_groups.library_selection,analysis.metadata.read_groups.read_group_name

die() {
    printf 'fetch_gdc.sbatch: %s\n' "$*" >&2
    exit 1
}
# From CTX_REPO, never from $0: sbatch runs a copy of this file from its spool folder. The site
# file says whether a folder may sit in home; its activation puts the fetch environment's
# Python, which the end check runs on, first on PATH.
CTX_SITE="${CTX_SITE:-kebnekaise}"
[[ "$CTX_SITE" =~ ^[a-z0-9_]+$ ]] || die "CTX_SITE is a site name such as kebnekaise or laptop, not $CTX_SITE"
SITE_FILE="${CTX_REPO:?export CTX_REPO=<your clone of the course repository>}/scripts/bash/site_$CTX_SITE.sh"
[[ -f "$SITE_FILE" ]] || die "no $SITE_FILE; is CTX_REPO your clone of the course repository, and CTX_SITE a site it has?"
source "$SITE_FILE"
source "$CTX_REPO/scripts/bash/job_name.sh"
[[ "$PROJECT" =~ ^[A-Z0-9]+(-[A-Z0-9]+)+$ ]] || die "a GDC project ID is capitals and digits joined by hyphens, not $PROJECT"
[[ "$DATA" == [01] ]] || die "DATA takes 0 or 1"
[[ "$CDR" == auto || "$CDR" == [01] ]] || die "CDR takes auto, 0 or 1"
[[ "$CDR" != auto ]] || { [[ "$PROJECT" == TCGA-* ]] && CDR=1 || CDR=0; }
[[ "$JOBS" =~ ^[1-9][0-9]?$ ]] || die "JOBS takes 1 to 99"
[[ "$RETRIES" =~ ^[0-9]$ ]] || die "RETRIES takes 0 to 9"
for kind in $KINDS; do
    [[ "$kind" == star || "$kind" == mirna || "$kind" == isoform ]] || die "KINDS takes star, mirna and isoform, not $kind"
done
[[ "$RUN" == /* ]] || RUN="$PWD/$RUN"
job_rename "${RUN##*/}" gdc-download
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
site_activate || die "the site's activation failed"
PY=$(command -v python3) || die "no python3 on PATH after the site's activation"

fetch() { curl -sS --fail --retry "$RETRIES" --retry-all-errors --connect-timeout 30 --max-time 900 "$@"; }
md5of() { if command -v md5sum >/dev/null; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }
# One row per retrieval, so every file carries its own retrieval date and the release it came from.
log() { printf '%s\t%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$NOW" "$1" "$2" >>"$OUT/_retrieval_log.tsv"; }
stamp() { date -u +%Y%m%dT%H%M%SZ; }

# download <file_id> <file_name> <md5>: into gdc-client's layout <file_id>/<file_name>, through
# .part, moved into place only when the md5 matches. A finished file whose md5 no longer
# matches is moved aside, never deleted.
download() {
    local dest="$OUT/$1/$2"
    if [[ -f "$dest" && "$(md5of "$dest")" == "$3" ]]; then return 0; fi
    mkdir -p "$OUT/$1"
    if [[ -f "$dest" ]]; then mv "$dest" "$dest.md5_mismatch_$(stamp)"; fi
    fetch -o "$dest.part" "$GDC_DATA/$1" || { echo "FAIL download $1"; return 1; }
    if [[ "$(md5of "$dest.part")" == "$3" ]]; then
        mv "$dest.part" "$dest" && log "$1/$2" "GET $GDC_DATA/$1"
        return 0
    fi
    echo "FAIL md5 $dest.part: $(md5of "$dest.part"), the manifest says $3"
    return 1
}
export -f download fetch md5of log stamp
export GDC_DATA RETRIES

# post <file> <endpoint> <json>: fetched once; its query is saved beside it as <file>.query.json.
post() {
    if [[ -s "$OUT/$1" ]]; then echo "kept $1 from an earlier run"; return 0; fi
    printf '%s\n' "$3" >"$OUT/$1.query.json"
    fetch -H 'Content-Type: application/json' --data-binary @"$OUT/$1.query.json" -o "$OUT/$1.part" "$GDC_API/$2" || { echo "FAIL $1"; return 1; }
    mv "$OUT/$1.part" "$OUT/$1" && log "$1" "POST $GDC_API/$2 with $1.query.json"
}
get() {
    if [[ -s "$OUT/$1" ]]; then echo "kept $1 from an earlier run"; return 0; fi
    fetch -o "$OUT/$1.part" "$2" || { echo "FAIL $1"; return 1; }
    mv "$OUT/$1.part" "$OUT/$1" && log "$1" "GET $2"
}
release() { sed -n 's/.*"data_release": *"\([^"]*\)".*/\1/p' "$1"; }

# one_kind <kind>: metadata, then data unless DATA=0, then the end check; its exit status.
one_kind() {
    local kind=$1 data_type workflow name
    case "$kind" in
    star) data_type="Gene Expression Quantification" workflow="STAR - Counts" name=star_counts OUT=$RUN ;;
    mirna) data_type="miRNA Expression Quantification" workflow="BCGSC miRNA Profiling" name=mirna OUT=${RUN}_mirna ;;
    isoform) data_type="Isoform Expression Quantification" workflow="BCGSC miRNA Profiling" name=isoform OUT=${RUN}_isoform ;;
    esac
    export OUT
    local files='{"op":"and","content":[{"op":"in","content":{"field":"cases.project.project_id","value":["'"$PROJECT"'"]}},{"op":"in","content":{"field":"data_type","value":["'"$data_type"'"]}},{"op":"in","content":{"field":"analysis.workflow_type","value":["'"$workflow"'"]}},{"op":"in","content":{"field":"access","value":["open"]}}]}'
    local cases='{"op":"in","content":{"field":"project.project_id","value":["'"$PROJECT"'"]}}'
    local reads='{"op":"and","content":[{"op":"in","content":{"field":"cases.project.project_id","value":["'"$PROJECT"'"]}},{"op":"in","content":{"field":"experimental_strategy","value":["RNA-Seq"]}},{"op":"in","content":{"field":"data_type","value":["Aligned Reads"]}}]}'
    local query="{\"filters\":$files,\"size\":10000,\"format\":\"TSV\",\"fields\":\"$FIELDS\"}"
    # A folder fetched for another project or workflow would mix two file lists.
    local saved="$OUT/gdc_files_$name.tsv.query.json"
    if [[ -s "$saved" ]] && ! cmp -s "$saved" <(printf '%s\n' "$query"); then
        echo "FAIL $OUT holds another query's files ($saved); give this project its own folder"
        return 1
    fi
    echo "PASTE start $(date -u +%Y-%m-%dT%H:%M:%SZ): $PROJECT $kind into $OUT"
    # This script is kept beside what it fetched; an older, different copy is moved aside.
    if [[ -f "$OUT/_fetch_gdc.sh" ]] && ! cmp -s "$0" "$OUT/_fetch_gdc.sh"; then mv "$OUT/_fetch_gdc.sh" "$OUT/_fetch_gdc.sh.superseded_$(stamp)"; fi
    [[ -f "$OUT/_fetch_gdc.sh" ]] || cp -p "$0" "$OUT/_fetch_gdc.sh"
    # The release GDC serves now, and the one the saved metadata came from (the first run's).
    fetch -o "$OUT/_status_last_run.json" "$GDC_API/status" || { echo "FAIL: no answer from $GDC_API/status"; return 1; }
    NOW=$(release "$OUT/_status_last_run.json")
    export NOW
    [[ -s "$OUT/_retrieval_log.tsv" ]] || printf 'retrieved_utc\tgdc_release\tfile\tsource\n' >"$OUT/_retrieval_log.tsv"
    [[ -s "$OUT/gdc_status.json" ]] || { cp "$OUT/_status_last_run.json" "$OUT/gdc_status.json" && log gdc_status.json "GET $GDC_API/status"; }
    [[ "$NOW" == "$(release "$OUT/gdc_status.json")" ]] || echo "WARNING: GDC changed release since the first run; the saved manifest still decides what is fetched. For the new release, move $OUT aside and run again."

    # Manifest and file metadata first, then the clinical side, so all metadata is from one moment.
    post "gdc_manifest_$name.txt" files "{\"filters\":$files,\"size\":10000,\"return_type\":\"manifest\"}"
    post "gdc_files_$name.tsv" files "$query"
    if [[ "$kind" == star ]]; then
        post gdc_cases_clinical.json cases "{\"filters\":$cases,\"size\":10000,\"format\":\"JSON\",\"expand\":\"$EXPAND\"}"
        post gdc_annotations.json annotations "{\"filters\":$cases,\"size\":10000,\"format\":\"JSON\"}"
        post gdc_read_groups.json files "{\"filters\":$reads,\"size\":10000,\"format\":\"JSON\",\"fields\":\"$READ_FIELDS\"}"
        # The CDR covers TCGA projects alone (CDR=auto); its md5 is the PanCan manifest's.
        if [[ "$CDR" == 1 ]]; then
            get "$PANCAN" "$GDC_WEB/system/files/public/file/$PANCAN"
            local cdr_name cdr_md5
            cdr_name=$(awk -F'\t' -v id="$CDR_ID" '$1 == id {print $2}' "$OUT/$PANCAN" 2>/dev/null)
            cdr_md5=$(awk -F'\t' -v id="$CDR_ID" '$1 == id {print $3}' "$OUT/$PANCAN" 2>/dev/null)
            if [[ -n "$cdr_md5" && "$cdr_name" =~ ^[A-Za-z0-9._-]+$ ]]; then
                # Beside the metadata, where [gdc] cdr_file names it, not in a file_id folder.
                if [[ ! -f "$OUT/$cdr_name" || "$(md5of "$OUT/$cdr_name")" != "$cdr_md5" ]]; then
                    if [[ -f "$OUT/$cdr_name" ]]; then mv "$OUT/$cdr_name" "$OUT/$cdr_name.md5_mismatch_$(stamp)"; fi
                    fetch -o "$OUT/$cdr_name.part" "$GDC_DATA/$CDR_ID" &&
                        [[ "$(md5of "$OUT/$cdr_name.part")" == "$cdr_md5" ]] &&
                        mv "$OUT/$cdr_name.part" "$OUT/$cdr_name" && log "$cdr_name" "GET $GDC_DATA/$CDR_ID"
                fi
                [[ -f "$OUT/$cdr_name" ]] || echo "FAIL: the TCGA-CDR $cdr_name was not saved with the md5 $PANCAN lists"
            else
                echo "FAIL: $CDR_ID is not in $PANCAN"
            fi
        fi
    fi

    # The data files. xargs treats blanks, quotes and backslashes specially: allow only safe names.
    local m="$OUT/gdc_manifest_$name.txt"
    [[ -s "$m" ]] || { echo "FAIL: no manifest at $m"; return 1; }
    [[ "$(head -n 1 "$m")" == "$(printf 'id\tfilename\tmd5\tsize\tstate')" ]] || { echo "FAIL: the manifest header is not id, filename, md5, size, state"; return 1; }
    if awk -F'\t' 'NR > 1 && (NF != 5 || $1 !~ /^[0-9a-f-]+$/ || $2 ~ /[^A-Za-z0-9._-]/)' "$m" | grep -q .; then
        echo "FAIL: a manifest row is not 5 fields, or its id or name holds a character other than letters, digits, dot, underscore or hyphen"
        return 1
    fi
    local n bytes
    n=$(awk 'NR > 1' "$m" | wc -l | tr -d ' ')
    bytes=$(awk -F'\t' 'NR > 1 {s += $4} END {printf "%d", s}' "$m")
    [[ "$n" -gt 0 ]] || { echo "FAIL: the manifest lists no file of $PROJECT $kind"; return 1; }
    if [[ "$DATA" == 0 ]]; then
        echo "PASTE bytes of GDC files, $kind: $n files, $bytes bytes"
        TOTAL=$((TOTAL + bytes))
        return 0
    fi
    echo "fetching $n data files, $bytes bytes, $JOBS at a time"
    awk -F'\t' 'NR > 1 {print $1, $2, $3}' "$m" | xargs -n 3 -P "$JOBS" bash -c 'download "$@"' _
    find "$OUT" -name '*.md5_mismatch_*' | sed 's/^/moved aside for inspection: /'
    "$PY" - "$OUT" "$kind" "$name" "$RUN/gdc_files_star_counts.tsv" "$CDR_ID" "$CDR" <<'PY'
# Every check a folder must pass before gather reads it; each failure is named, then exit 1.
import collections, csv, hashlib, json, os, sys

out, kind, name, star_tsv, cdr_id, want_cdr = sys.argv[1:7]
got, bad = {}, []


def path(f):
    return os.path.join(out, f)


def counts(c):
    return ", ".join("%s %d" % kv for kv in sorted(c.items())) or "none"


gates = []


def must(check, ok, value):
    got[check] = value
    gates.append(check)
    if not ok:
        bad.append(check)


def md5(f):
    h = hashlib.md5()
    with open(f, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def release(f):
    try:
        return json.load(open(path(f)))["data_release"]
    except (OSError, ValueError, KeyError):
        return "MISSING " + f


now, saved = release("_status_last_run.json"), release("gdc_status.json")
must("release now, of the saved metadata", now == saved, "%s; %s" % (now, saved))
manifest = list(csv.DictReader(open(path("gdc_manifest_%s.txt" % name), newline=""), delimiter="\t"))
ok = mismatched = gone = on_disk = listed = 0
for row in manifest:
    listed += int(row["size"])
    f = path(os.path.join(row["id"], row["filename"]))
    if not os.path.isfile(f):
        gone += 1
        print("missing: " + f)
    elif md5(f) == row["md5"]:
        ok += 1
        on_disk += os.path.getsize(f)
    else:
        mismatched += 1
        print("bad md5: " + f)
must("data files md5 OK, bad, missing", mismatched == gone == 0, "%d, %d, %d" % (ok, mismatched, gone))
must("bytes on disk, in the manifest", on_disk == listed, "%d, %d" % (on_disk, listed))
parts = sum(f.endswith(".part") for _, _, fs in os.walk(out) for f in fs)
must(".part leftovers", parts == 0, str(parts))


def pairs_of(tsv):
    """Per file: case, sample type, preservation; by column name (GDC ends lines with CRLF)."""
    rows = list(csv.DictReader(open(tsv, newline=""), delimiter="\t"))
    tum = {r["cases.0.submitter_id"] for r in rows if r["cases.0.samples.0.sample_type"] == "Primary Tumor"}
    nor = {r["cases.0.submitter_id"] for r in rows if r["cases.0.samples.0.sample_type"] == "Solid Tissue Normal"}
    return rows, tum & nor


rows, pairs = pairs_of(path("gdc_files_%s.tsv" % name))
ids = {r["id"] for r in manifest}
outside = sum(r["file_id"] not in ids for r in rows)
must("metadata files not in the manifest", outside == 0 and len(rows) == len(manifest), "%d of %d" % (outside, len(rows)))
got["metadata files, cases"] = "%d, %d" % (len(rows), len({r["cases.0.submitter_id"] for r in rows}))
got["sample types"] = counts(collections.Counter(r["cases.0.samples.0.sample_type"] for r in rows))
got["pairs, cases with a primary tumour and a solid tissue normal"] = str(len(pairs))
got["preservation"] = counts(collections.Counter(r["cases.0.samples.0.preservation_method"] for r in rows))
if kind != "star":
    try:
        got["pairs also paired in STAR - Counts"] = str(len(pairs & pairs_of(star_tsv)[1]))
    except (OSError, KeyError) as e:
        got["pairs also paired in STAR - Counts"] = "not measured: %s (%s)" % (star_tsv, e)
else:
    for record in ("gdc_cases_clinical.json", "gdc_annotations.json", "gdc_read_groups.json"):
        try:
            found = json.load(open(path(record)))
            hits, total = len(found["data"]["hits"]), found["data"]["pagination"]["total"]
            warnings = found.get("warnings") or {}
            must(record + ", warnings", hits == total and not warnings, "%d of %d; %s" % (hits, total, warnings))
        except (OSError, ValueError, KeyError) as e:
            must(record + ", warnings", False, "MISSING (%s)" % e)
    if want_cdr == "1":
        try:
            cdr = [l.split("\t") for l in open(path("PanCan-General_Open_GDC-Manifest_2.txt")).read().splitlines()]
            n, want = [r for r in cdr if r[0] == cdr_id][0][1:3]
            have = md5(path(n))
            must("CDR", have == want, "%s md5 %s bytes %d" % (n, have, os.path.getsize(path(n))))
        except (OSError, IndexError) as e:
            must("CDR", False, "MISSING (%s)" % e)
for check, value in got.items():
    print("PASTE %s: %s%s" % (check, value, "   FAILS" if check in bad else ""))
print("PASTE checks: %d of %d pass" % (len(gates) - len(bad), len(gates)))
sys.exit(1 if bad else 0)
PY
}

T0=$(date +%s)
TOTAL=0
STATUS=0
for kind in $KINDS; do
    folder=$RUN
    [[ "$kind" == star ]] || folder="${RUN}_$kind"
    mkdir -p "$folder" || die "cannot make $folder"
    one_kind "$kind" 2>&1 | tee -a "$folder/_fetch_gdc.log"
    [[ "${PIPESTATUS[0]}" == 0 ]] || STATUS=1
done
if [[ "$DATA" == 0 ]]; then
    # one_kind ran in a pipe's subshell, so the total is read back from each folder's manifest.
    for kind in $KINDS; do
        case "$kind" in star) f="$RUN/gdc_manifest_star_counts.txt" ;; *) f="${RUN}_$kind/gdc_manifest_$kind.txt" ;; esac
        [[ -s "$f" ]] && TOTAL=$((TOTAL + $(awk -F'\t' 'NR > 1 {s += $4} END {printf "%d", s}' "$f")))
    done
fi
echo "PASTE wall time $(($(date +%s) - T0)) s"
if [[ "$DATA" == 0 ]]; then
    echo "PASTE bytes of GDC files, $KINDS: $TOTAL bytes"
else
    echo "PASTE $PROJECT: $([[ $STATUS == 0 ]] && echo 'every check passed' || echo 'a check failed; run the same command again, then read the lines marked FAIL')"
fi
exit $STATUS
