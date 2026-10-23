# ED.6 M2 and M3 on the four core projects: Kebnekaise blocks (written 2026-10-22; block A run)

Superseded 2026-10-23 by D-113 and D-116: the read-level pages read nf-core/rnaseq's outputs; `ctx.process.cluster`, `cluster_route`, `fc_compact` and `ctx.qc.cluster`, which these blocks fed, are removed.

Rytis decided on 2026-10-22: in the core edition, all four projects get M2 (QC, its read stage
included) and M3 (alignment) on their real runs. TODO.md ED.6 holds the items; this file holds the
blocks. Block A was run by Rytis on 2026-10-22 at 12:27; its numbers are below. No other block
has run. Every number below that was not measured is marked "estimate" with its basis.

pancreas_gse192647 (PRJNA792551, TODO.md ED.8) was added on 2026-10-22, after block A had run, as
a sixth slug: Rytis decided the core edition covers breast, lung and pancreas. Its run folder is
`/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/pancreas_gse192647`, its download and
count queued by Rytis that day; block A never saw it.

Each login-node block writes one `set -euo pipefail` script into the run folder and runs it with
`bash`, so a failure stops that script and never the login shell (ED.4b's pattern). The run folder
is `/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/ed6_align`, on the project storage,
never `$SNIC_NOBACKUP`. The compute project is `hpc2n2026-161`.

## What is fixed before any block runs

**The projects.** Run folders as the earlier blocks made them; FASTQs in
`<run folder>/nfcore_fetchngs/results/fastq/`, with fetchngs's `md5/<file>.md5` beside them and its
`samplesheet/samplesheet.csv`, whose `sample` column is the experiment accession (equal on all 210
Step-1 rows and all 102 breast_mono rows, read on the Mac 2026-10-22).

| slug | BioProject | runs | samples | layout | read length | FASTQ bytes | basis of the bytes |
|---|---|---|---|---|---|---|---|
| tier1 | PRJEB49000 | 210 | 210 | single | 75, 76 or 101 nt | 287.1 GB | block A, 2026-10-22 |
| breast_tissue | PRJNA393443 | 103 | 103 | single | 50 nt (80), 100 nt (23) | 135,964,972,806 | ENA filereport, 2026-09-30 (P4.0) |
| breast_mac | PRJNA483877 | 22 | 22 | paired | 2 x 100 nt | 137,937,197,084 | ENA filereport, 2026-09-30 (P4.0) |
| breast_mono | PRJNA483877 | 102 | 94 | paired | 2 x 100 nt | 580,591,689,686 | ENA fastq_bytes, 2026-10-01 (P4.4) |
| lung_gse229705 | PRJNA955664 | 246 | 246 | paired | 2 x 51 nt | 385,618,590,287 | ENA fastq_bytes, 2026-10-22 (ED.4b) |
| pancreas_gse192647 | PRJNA792551 | 28 | 28 | paired | 2 x 101 nt | 127,118,223,285 | ENA fastq_bytes, 2026-10-22 (ED.8a) |

Read lengths are the dataset files' `[expect.counts.read_length]`; the lung length is ED.4a's, the
pancreas length ED.8a's.

**Block A's measurements (Rytis, Kebnekaise, 2026-10-22 12:27).** The storage
`/pfs/proj/nobackup/fs/projnb10` had 21T free of 36T. Every project's FASTQs are in
`/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/<slug>/nfcore_fetchngs/results/fastq`:

| slug | FASTQ files | bytes | layout |
|---|---|---|---|
| tier1 | 210 | 287.1 GB | single-end |
| breast_tissue | 103 | 136.0 GB | single-end |
| breast_mac | 44 (22 read 2) | 137.9 GB | paired |
| breast_mono | 204 (102 read 2) | 580.6 GB | paired |
| lung_gse229705 | 34 of 492, downloading | | paired |
| pancreas_gse192647 | not checked: added after block A | | paired |

The four complete projects match the table above to its rounding. The lung and pancreas pilots
use a sample already downloaded: while `results.complete` is absent, block B1 lists only the
samples whose every file and md5 file is present, so each is the largest downloaded sample, not
the project's. Lung's download goes through ctx's `fetchngs.sbatch`, which writes
`results.complete` (ED.4b; the lab launcher writes none, so the four complete projects are not
tested for it). How pancreas's was launched is not recorded here: if its folder never gets
`results.complete`, block C1 refuses to write its full task list, and the marker is checked by hand.

**The tools.** The course's `reads` pixi environment, installed from the committed `pixi.lock`
on Kebnekaise (`pixi install -e reads --frozen`), so every tool is the build ctx itself runs and
`ctx.tools` accepts: STAR 2.7.11b (`h5ca1c30_8`, the build M3.11's depot image carries), fastp 1.1.0
(`heae3180_0`), FastQC 0.12.1, seqkit 2.13.0, samtools 1.22.1, subread 2.1.1. This replaces M3.11's
depot image with the same STAR build and adds the other five tools from one source.

**fastp, as M2 and M3 run it (M3.9: one setting for both sides).** `ctx.qc.Fastp()`'s defaults,
the defaults of both `ctx qc` and `ctx process`: `-q 15 -u 40 -l 15`, no `--cut_tail`. The call is
`ctx.qc.reads.run_fastp`'s, per run, titled by the run:
`fastp -i <r1> [-I <r2>] -R <run> -j <run>.fastp.json -h <html> -w <threads> -q 15 -u 40 -l 15 [-o ... [-O ...]]`.
M3 writes the trimmed reads (`-o`, `-O`), M2 does not; one call with `-o` serves both, and the pilot
checks that the JSON is the same either way. One departure: `-h <run>.fastp.html` where ctx writes
`-h /dev/null`, because ED.6 keeps the HTML; it names a report file and changes no read and no JSON
number (the JSON's `command` string differs).

**FastQC and seqkit, as M2's read stage runs them (`reads.run_tools`).** `fastqc --quiet --threads
<n> --outdir <dir> <raw files>`; `seqkit stats -a -T -j <n> -o <tsv> <raw files>`. ctx writes one
seqkit table for all files; here one per sample (ED.6e-3). M2 reads FastQC through MultiQC's
FastQC module, which reads the `_fastqc.zip` files (`ctx.qc.fastqc`, D-028), so the zips are kept.

**The grid (M2's `reads.fastp_grid`).** fastp, no `-o`, on one run at every `-q` in
`GRID_Q = (10, 15, 20, 25, 30, 35)` and `-l` in `GRID_L = (15, 25, 36, 50, 75, 102)`, `-u 40`
fixed, titled `grid_q<q>_l<l>`: 36 runs per project.

**STAR, as M3 runs it.** The index is M3.11's block 2: `--runMode genomeGenerate --runThreadN 8
--sjdbOverhang 100 --genomeSAindexNbases 14` from the MD5-checked GENCODE 50 genome and GTF,
unzipped, plus M3.11's one departure `--limitGenomeGenerateRAM 42000000000`; 100 is
`reference.OVERHANG`, 14 is `sa_index_bases` for the 3.1 Gb genome. Mapping is `align.star_align`
with `STAR_OPTIONS`: `--runThreadN <n> --genomeDir <index> --readFilesIn <r1,...> [<r2,...>]
--readFilesCommand gzip -dcf --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate
--outSAMunmapped Within --outFileNamePrefix <dir>/<sample>.`, no `--outFilterMultimapNmax`
(D-036). One STAR run per sample, its runs comma-joined in run-accession order (breast_mono's 8
two-run samples), as `align.by_sample` joins them. The BAM is renamed `<sample>.bam`
(`align.bam_name`), then `samtools index`, and `samtools flagstat` and `idxstats` into
`<sample>.flagstat` and `<sample>.idxstats` (`align.bam_stats`).

**featureCounts, as M3 runs it (`count.FcOptions.args`, `count.feature_counts`).** `featureCounts
-T <n> -a <gtf> -s <k> [-p --countReadPairs] [-O] [-M] [--fraction] -o <out> <sample>.bam`, run in
the BAM's folder so the column is headed `<sample>.bam`. M3 counts `-s 0`, `-s 1`, `-s 2`
(`plain_settings`), then `-O`, `-O --fraction` and `-M` at each sample's called strand
(`variant_settings`). The strand is called by ctx from the plain counts through the registry's
`gene_strand_bias` bands (`count.called_strands`); a second copy of that rule on the cluster would
be a second definition, so the cluster counts each variant at all three strands, 12 runs per
sample, and the Mac takes the one at the called strand (ED.6e-5). The pilot measures their cost.

**What M3 reads from a BAM, and so what the BAM is kept for.** featureCounts (all 12 settings) and
samtools (index, flagstat, idxstats). Nothing else: M3 reads `Log.final.out`, `ReadsPerGene.out.tab`,
the flagstat and idxstats and the featureCounts matrices and summaries, and no code under
`src/ctx` reads `SJ.out.tab` (grep, 2026-10-22). So each task writes those tables and then deletes
its BAM, its index and its trimmed reads; the pilot keeps them. `SJ.out.tab`, `Log.out` and
`Log.progress.out` stay on the cluster and are not copied.

**The reference.** GENCODE 50 in `lung_gse229705/refs/` (ED.4b block 4 downloaded `MD5SUMS`, the
genome and the GTF there). If M3.11's index already exists in `tier1/refs/`, it is used and none is
built; block A says which. featureCounts reads the unzipped GTF the index was built from.

## Costs: what is measured and what is not

**The pilot's numbers (Rytis, Kebnekaise, 2026-10-22), as the coordinator relayed them.** STAR index:
45 min, 44.7 GB MaxRSS, 30.7 GB on disk. breast_tissue pilot (3.3 GB single-end): 18 min, 41.0 GiB
MaxRSS, STAR 4.5 min, featureCounts 12 x 30 s, peak disk 16.1 GB, strand reverse from Assigned
(`-s 2` 45.8 M against `-s 1` 4.0 M). breast_mono pilot (10.9 GB paired, two runs): 38.7 min,
40.1 GiB, STAR 12.4 min, featureCounts 12 x 73 s, peak disk 36.4 GB, unstranded (`-s 0` 13.4 M,
`-s 1` 8.1 M, `-s 2` 8.0 M). The other three pilots' numbers are not in this file yet.

**The featureCounts matrices.** One sample's matrix at one setting is 133 MB, 1.6 GB per sample at
12 settings, about 1.1 TB for 683 samples to copy back. Nearly all of it is the annotation columns
(Chr, Start, End, Strand, Length), the same in every file of one annotation, so blocks B1, C1 and E
keep each matrix compact (`ctx.process.fc_compact`): the command line, Geneid and the count, gzipped,
and the annotation columns once per project (and once for block E's half), and delete the original
only after the compact file reads back equal to it, row by row and in its sum. `ctx process
--align-dir` rebuilds featureCounts' matrix from the two. Measured on the test fixture (10 genes):
63 compact files of 369 to 473 bytes; the annotation copy 856 bytes, 6,463 uncompressed (423 for the half). Estimated on
GENCODE 50, not measured on a cluster file: 78,941 gene IDs from the GTF with one tier1 sample's
gene counts (ERX6994339, the gather counts, Salmon's, standing in for featureCounts') make a compact
file of 1,620,719 bytes, 361,462 gzipped; so about 4.3 MB per sample at 12 settings and about 3.0 GB
for 683 samples, plus one annotation copy per project. Split counts (`-O --fraction`) print decimals
and will be larger; the first compact pilot file replaces this estimate.

Nothing in this table is measured; the pilot replaces every row.

| Step | Pilot request | Basis |
|---|---|---|
| STAR index | 8 cores, 48G, 3 h | M3.11 block 2's request, itself unmeasured: RAM from STAR's README (32 GB), time a guess |
| One sample, every step | 8 cores, 48G, 12 h | 8 cores and 48G copy M3.11's alignment request; 12 h is a ceiling guess for the largest sample of each project, since no full-size STAR run has been measured in this project |
| Trimmed reads per task | about the FASTQ's bytes | estimate: fastp recompresses at its default level; unmeasured |
| BAM per task | unknown | M3's plan gives about 31 MB per 1M reads on chr22, unmeasured on the full genome |
| Kept per sample | unknown | 12 featureCounts matrices carry GENCODE's per-exon coordinates per gene; unmeasured |

The pilot runs each project's **largest** sample by FASTQ bytes, so its wall time, RSS and disk
bound the project's other samples, and its bytes give the rate for the totals.

---

## Block A (ED.6a). Kebnekaise login node: what is on the cluster

Produces `ed6_align/ed6a_check.log`. Changes nothing. Paste back the log.

```
mkdir -p /proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/ed6_align && cd /proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/ed6_align
cat > ed6a_check.sh <<'SCRIPT'
set -euo pipefail
S=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial
df -h /proj/nobackup/hpc2n2024-150
df -B1 --output=avail /proj/nobackup/hpc2n2024-150 | tail -n 1 | sed 's/^ */free bytes: /'
for p in tier1 breast_tissue breast_mac breast_mono lung_gse229705 pancreas_gse192647; do
    R="$S/$p/nfcore_fetchngs/results"
    echo "== $p"
    if [[ ! -d "$R/fastq" ]]; then
        echo "$p: NO $R/fastq"
        ls -la "$S/$p" "$S/$p/nfcore_fetchngs" 2>&1 | head -n 30 || true
        continue
    fi
    if [[ -f "$S/$p/nfcore_fetchngs/results.complete" ]]; then echo "$p: results.complete $(cat "$S/$p/nfcore_fetchngs/results.complete")"; else echo "$p: no results.complete (the lab launcher writes none; ctx's fetchngs.sbatch does)"; fi
    find "$R/fastq" -maxdepth 1 -name '*.fastq.gz' -printf '%s %n\n' | awk -v p="$p" '{n++; s+=$1; if ($2 == 1) one++} END {printf "%s: %d FASTQ files, %.0f bytes, %d not hard-linked\n", p, n, s, one}'
    printf '%s: %s md5 files\n' "$p" "$(find "$R/fastq/md5" -maxdepth 1 -name '*.md5' 2>/dev/null | wc -l)"
    python3 - "$R" "$p" <<'PY'
import csv, os, sys
r, p = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(os.path.join(r, "samplesheet", "samplesheet.csv"))))
paired = sum(bool(x["fastq_2"]) for x in rows)
files = [os.path.basename(x[k]) for x in rows for k in ("fastq_1", "fastq_2") if x[k]]
present = sum(os.path.isfile(os.path.join(r, "fastq", f)) for f in files)
md5 = sum(os.path.isfile(os.path.join(r, "fastq", "md5", f + ".md5")) for f in files)
print(f"{p}: sheet {len(rows)} runs, {len({x['sample'] for x in rows})} samples, {paired} paired, "
      f"{len(rows) - paired} single; {present} of {len(files)} sheet files on disk, {md5} with an md5 file")
PY
    f=$(find "$R/fastq" -maxdepth 1 -name '*.fastq.gz' | sort | head -n 1)
    if [[ -n "$f" ]]; then
        printf '%s: first read of %s is %s nt (one read, not a survey)\n' "$p" "$(basename "$f")" "$( { gzip -dc "$f" 2>/dev/null || true; } | awk 'NR == 2 {print length($0); exit}')"
    fi
    if [[ -d "$S/$p/nfcore_rnaseq/results/fastqc/raw" ]]; then
        printf '%s: %s raw-read FastQC zips in the nf-core outdir\n' "$p" "$(find "$S/$p/nfcore_rnaseq/results/fastqc/raw" -name '*_fastqc.zip' | wc -l)"
    fi
done
echo "== references and indices"
for d in "$S/tier1/refs" "$S/lung_gse229705/refs"; do ls -la "$d" 2>&1 | head -n 20 || true; done
for i in "$S/tier1/refs/star_2.7.11b_gencode50_sjdb100" "$S/lung_gse229705/refs/star_2.7.11b_gencode50_sjdb100"; do
    if [[ -s "$i/SA" && -s "$i/genomeParameters.txt" ]]; then echo "STAR index present: $i, $(du -sb "$i" | cut -f1) bytes"; else echo "no STAR index at $i"; fi
done
echo "== the course clone and its environments"
git -C "$S/RNAseq_tutorial" log -1 --format='%h %s' || echo "no clone at $S/RNAseq_tutorial"
ls "$S/RNAseq_tutorial/.pixi/envs" 2>/dev/null || echo "no pixi environments in the clone"
squeue -u "$USER"
SCRIPT
bash ed6a_check.sh 2>&1 | tee ed6a_check.log
```

Expected, per project present: 210, 103, 44, 204 and 492 FASTQ files (lung only once ED.4b's
download has ended), as many md5 files, every sheet file on disk, the layout of the table above,
and the free space. **Stop and paste** if a project's folder is missing or a count differs: a
missing project needs its download again (P2A.3, P4.0 or P4.4 for the breast and Step-1 runs,
through `fetchngs.sbatch` as ED.4b does), which is Rytis's decision, not this plan's.

## Block B1 (ED.6b). Kebnekaise login node: tools, task lists, the index and the pilots

Takes the projects to pilot as arguments, so lung runs later by the same script. Produces
`ed6_align/{bin.txt,versions.txt,refs.txt,jobs.txt,ed6_index.sbatch,ed6_align.sbatch}`,
`ed6_align/<slug>/tasks.tsv`, and at most one index job and one pilot job per project.

```
mkdir -p /proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/ed6_align && cd /proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/ed6_align
cat > ed6b_1_pilot.sh <<'SCRIPT'
set -euo pipefail
S=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial
E="$S/ed6_align"
A=hpc2n2026-161
CTX_REPO="$S/RNAseq_tutorial"
PIXI="${PIXI_HOME:-$HOME/.pixi}/bin/pixi"
[[ $# -ge 1 ]] || { echo "usage: bash ed6b_1_pilot.sh <slug> [<slug> ...]"; exit 1; }
cd "$E" && mkdir -p logs && touch jobs.txt
df -h /proj/nobackup/hpc2n2024-150
git -C "$CTX_REPO" pull --ff-only
git -C "$CTX_REPO" log -1 --format='%h %s'
PIXI_CACHE_DIR="$S/pixi_cache" "$PIXI" install --manifest-path "$CTX_REPO/pixi.toml" -e reads --frozen
BIN="$("$PIXI" run --manifest-path "$CTX_REPO/pixi.toml" -e reads --frozen printenv CONDA_PREFIX)/bin"
echo "$BIN" > bin.txt
{
    echo "STAR $("$BIN/STAR" --version)"
    "$BIN/fastp" --version 2>&1
    "$BIN/fastqc" --version
    "$BIN/seqkit" version
    "$BIN/samtools" --version | sed -n 1p
    { "$BIN/featureCounts" -v 2>&1 || true; } | sed -n '/featureCounts v/p'
} | tee versions.txt
# The pins of pixi.toml's reads feature; any other version is refused, as ctx.tools refuses it.
for want in 'STAR 2.7.11b' 'fastp 1.1.0' 'FastQC v0.12.1' 'seqkit v2.13.0' 'samtools 1.22.1' 'featureCounts v2.1.1'; do
    grep -q "$want" versions.txt || { echo "NOT THE PIN: no '$want' in versions.txt"; exit 1; }
done

# One index for every project: M3.11's if it exists, else built in the lung run folder's refs.
if [[ -s "$S/tier1/refs/star_2.7.11b_gencode50_sjdb100/SA" ]]; then REFS="$S/tier1/refs"; else REFS="$S/lung_gse229705/refs"; fi
IDX="$REFS/star_2.7.11b_gencode50_sjdb100"
GTF="$REFS/unzipped/gencode.v50.primary_assembly.annotation.gtf"
printf 'REFS=%s\nIDX=%s\nGTF=%s\n' "$REFS" "$IDX" "$GTF" | tee refs.txt
# GENCODE's own MD5SUMS, and the two sums M3.11 recorded, must agree.
(cd "$REFS" && grep -E '  (GRCh38\.primary_assembly\.genome\.fa\.gz|gencode\.v50\.primary_assembly\.annotation\.gtf\.gz)$' MD5SUMS | md5sum -c -)
grep -q '^da1a11258be075cfa7af718162c894e7  GRCh38.primary_assembly.genome.fa.gz$' "$REFS/MD5SUMS"
grep -q '^289b91e5e95e8b0450d223246f10a12e  gencode.v50.primary_assembly.annotation.gtf.gz$' "$REFS/MD5SUMS"

cat > ed6_index.sbatch <<'EOF'
#!/usr/bin/env bash
#SBATCH -J ed6_star_index
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 0-03:00:00
#SBATCH -o slurm-%x-%j.out
# M3.11 block 2's index, the binary from the reads environment (same STAR build).
set -euo pipefail
: "${E:?}" "${REFS:?}" "${IDX:?}" "${BIN:?}"
U="$REFS/unzipped"
FA=GRCh38.primary_assembly.genome.fa
GTF=gencode.v50.primary_assembly.annotation.gtf
[ ! -e "$IDX" ] || { echo "$IDX exists: move it aside first"; exit 1; }
[ ! -e "$IDX.part" ] || { echo "$IDX.part exists, from a failed build: move it aside first"; exit 1; }
(cd "$REFS" && printf '%s  %s\n' da1a11258be075cfa7af718162c894e7 "$FA.gz" 289b91e5e95e8b0450d223246f10a12e "$GTF.gz" | md5sum -c -)
"$BIN/STAR" --version
mkdir -p "$U"
for f in "$FA" "$GTF"; do
    if [ ! -s "$U/$f" ]; then gzip -dc "$REFS/$f.gz" > "$U/$f.part"; mv "$U/$f.part" "$U/$f"; fi
done
ls -l "$U"
mkdir -p "$IDX.part"
/usr/bin/time -v -o "$E/logs/star_index.time" "$BIN/STAR" \
    --runMode genomeGenerate --runThreadN 8 --genomeDir "$IDX.part" \
    --genomeFastaFiles "$U/$FA" --sjdbGTFfile "$U/$GTF" \
    --sjdbOverhang 100 --genomeSAindexNbases 14 \
    --limitGenomeGenerateRAM 42000000000 --outFileNamePrefix "$IDX.part/"
mv "$IDX.part" "$IDX"
du -sb "$IDX"
EOF

cat > ed6_align.sbatch <<'EOF'
#!/usr/bin/env bash
#SBATCH -J ed6_align
#SBATCH -n 1
#SBATCH -o slurm-%x-%A_%a.out
# One sample per array task: fastp per run, FastQC and seqkit on the raw files, STAR, samtools,
# featureCounts at 12 settings; then, unless KEEP=1, the BAM and trimmed reads are deleted.
set -euo pipefail
: "${E:?}" "${PROJ:?}" "${IDX:?}" "${GTF:?}" "${BIN:?}" "${KEEP:?}" "${PILOT:?}" "${CTX_REPO:?}"
T="${SLURM_CPUS_PER_TASK:?submit with sbatch -c}"
[ -x /usr/bin/time ] || { echo "no /usr/bin/time on $(hostname)"; exit 1; }
FQ="/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/$PROJ/nfcore_fetchngs/results/fastq"
OUT="$E/$PROJ"
TASKS="$OUT/${TASKS:-tasks.tsv}"
line=$(awk -F'\t' -v n="$SLURM_ARRAY_TASK_ID" 'NR == n + 1' "$TASKS")
[ -n "$line" ] || { echo "no task $SLURM_ARRAY_TASK_ID in $TASKS"; exit 1; }
IFS=$'\t' read -r sample mates runs r1 r2 bytes <<< "$line"
D="$OUT/runs/$sample"
P="$D.part"
if [ -e "$D" ]; then echo "$sample: $D exists; nothing run"; exit 0; fi
# A failed try is moved aside, never overwritten, so its logs can be read; Rytis removes it.
if [ -e "$P" ]; then mv "$P" "$P.failed_${SLURM_JOB_ID}"; fi
mkdir -p "$P"/{fastp,fastqc,seqkit,star,featurecounts,trimmed,time,check}
( while :; do du -sb "$P" 2>/dev/null | cut -f1; sleep 60; done ) > "$P/time/du_samples.txt" &
SAMPLER=$!
trap 'kill "$SAMPLER" 2>/dev/null || true' EXIT
tm() { local name=$1; shift; /usr/bin/time -v -o "$P/time/$name.time" "$@"; }
IFS=, read -ra RUN <<< "$runs"
IFS=, read -ra F1 <<< "$r1"
F2=()
if [ "$mates" = 2 ]; then IFS=, read -ra F2 <<< "$r2"; fi
raw=()
for f in "${F1[@]}" "${F2[@]}"; do raw+=("$FQ/$f"); done
(cd "$FQ" && for f in "${F1[@]}" "${F2[@]}"; do md5sum -c "md5/$f.md5"; done) | tee "$P/check/md5.txt"
t1=(); t2=()
for i in "${!RUN[@]}"; do
    run=${RUN[$i]}
    given=(-i "$FQ/${F1[$i]}")
    if [ "$mates" = 2 ]; then given+=(-I "$FQ/${F2[$i]}"); fi
    # ctx.qc.reads.run_fastp's call; Fastp() = -q 15 -u 40 -l 15 (M2 and M3 alike, M3.9).
    if [ "$mates" = 2 ]; then
        trimmed=(-o "$P/trimmed/${run}_1.fastq.gz" -O "$P/trimmed/${run}_2.fastq.gz")
        t1+=("$P/trimmed/${run}_1.fastq.gz"); t2+=("$P/trimmed/${run}_2.fastq.gz")
    else
        trimmed=(-o "$P/trimmed/$run.fastq.gz")
        t1+=("$P/trimmed/$run.fastq.gz")
    fi
    tm "fastp_$run" "$BIN/fastp" "${given[@]}" -R "$run" -j "$P/fastp/$run.fastp.json" -h "$P/fastp/$run.fastp.html" \
        -w "$T" -q 15 -u 40 -l 15 "${trimmed[@]}" > "$P/fastp/$run.fastp.log" 2>&1
    if [ "$PILOT" = 1 ]; then
        # M2's own call, no reads written: its JSON must equal the one above but for the command.
        "$BIN/fastp" "${given[@]}" -R "$run" -j "$P/check/$run.fastp.json" -h /dev/null -w "$T" -q 15 -u 40 -l 15 > "$P/check/$run.fastp.log" 2>&1
        python3 - "$P/fastp/$run.fastp.json" "$P/check/$run.fastp.json" <<'PY' | tee -a "$P/check/fastp_equal.txt"
import json, sys
a, b = (json.load(open(p)) for p in sys.argv[1:])
a.pop("command", None); b.pop("command", None)
diff = sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))
print(f"{sys.argv[1]}: fastp JSON equal with and without -o, command aside: {not diff}; differing keys {diff}")
PY
    fi
done
tm fastqc "$BIN/fastqc" --quiet --threads "$T" --outdir "$P/fastqc" "${raw[@]}" > "$P/fastqc/fastqc.log" 2>&1
tm seqkit "$BIN/seqkit" stats -a -T -j "$T" -o "$P/seqkit/$sample.seqkit_stats.tsv" "${raw[@]}" > "$P/seqkit/seqkit.log" 2>&1
reads=("$(IFS=,; echo "${t1[*]}")")
if [ "$mates" = 2 ]; then reads+=("$(IFS=,; echo "${t2[*]}")"); fi
cd "$P/star"
# ctx.process.align.star_align with STAR_OPTIONS; no --outFilterMultimapNmax (D-036).
tm star "$BIN/STAR" --runThreadN "$T" --genomeDir "$IDX" --readFilesIn "${reads[@]}" \
    --readFilesCommand gzip -dcf --quantMode GeneCounts \
    --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
    --outFileNamePrefix "$P/star/$sample." > "$P/star/$sample.star.log" 2>&1
for n in Log.final.out ReadsPerGene.out.tab Aligned.sortedByCoord.out.bam; do
    test -s "$sample.$n" || { echo "STAR exited 0 and wrote no $sample.$n"; exit 1; }
done
mv "$sample.Aligned.sortedByCoord.out.bam" "$sample.bam"
tm samtools_index "$BIN/samtools" index "$sample.bam"
tm samtools_flagstat "$BIN/samtools" flagstat "$sample.bam" > "$sample.flagstat"
tm samtools_idxstats "$BIN/samtools" idxstats "$sample.bam" > "$sample.idxstats"
du -sb "$P/trimmed" "$P/star/$sample.bam" > "$P/time/du_before_delete.txt"
pe=()
if [ "$mates" = 2 ]; then pe=(-p --countReadPairs); fi
fc() {
    local name=$1; shift
    tm "fc_$name" "$BIN/featureCounts" -T "$T" -a "$GTF" "$@" -o "$P/featurecounts/$name.featureCounts.txt" "$sample.bam" \
        > "$P/featurecounts/$name.featureCounts.log" 2>&1
}
# count.FcOptions.args order: -s, then -p --countReadPairs, -O, -M, --fraction.
for s in 0 1 2; do
    fc "s$s" -s "$s" "${pe[@]}"
    fc "O_s$s" -s "$s" "${pe[@]}" -O
    fc "O_fraction_s$s" -s "$s" "${pe[@]}" -O --fraction
    fc "M_s$s" -s "$s" "${pe[@]}" -M
done
# One BAM's matrix is 133 MB on GENCODE 50 (pilot, 2026-10-22), nearly all annotation columns.
# ctx.process.fc_compact keeps Geneid and counts gzipped, the columns once per project, reads
# each compact file back against its original and only then deletes the original.
python3 "$CTX_REPO/src/ctx/process/fc_compact.py" "$OUT/featurecounts_annotation.tsv.gz" \
    "$P"/featurecounts/*.featureCounts.txt > "$P/featurecounts/compact.log" 2>&1
if [ "$KEEP" = 0 ]; then rm -f "$sample.bam" "$sample.bam.bai"; rm -rf "$P/trimmed"; fi
cd "$OUT"
kill "$SAMPLER" 2>/dev/null || true
du -sb "$P" > "$P/time/du_kept.txt"
mv "$P" "$D"
echo "done: $sample"
EOF

cat > ed6_tasks.py <<'PY'
# One row per sample, runs by accession, largest sample first. "partial" keeps only the samples
# whose every file and md5 file is present (a download still running); otherwise one absent file
# refuses the whole list.
import csv, os, sys
r, out = sys.argv[1], sys.argv[2]
partial = sys.argv[3:] == ["partial"]
rows = list(csv.DictReader(open(os.path.join(r, "samplesheet", "samplesheet.csv"))))
samples, absent = {}, {}
for x in sorted(rows, key=lambda x: x["run_accession"]):
    names = [os.path.basename(x[k]) for k in ("fastq_1", "fastq_2") if x[k]]
    for n in names:
        if not os.path.isfile(os.path.join(r, "fastq", n)) or not os.path.isfile(os.path.join(r, "fastq", "md5", n + ".md5")):
            absent.setdefault(x["sample"], []).append(n)
    samples.setdefault(x["sample"], []).append((x["run_accession"], names))
mixed = [s for s, rs in samples.items() if len({len(n) for _, n in rs}) != 1]
if mixed or (absent and not partial):
    sys.exit(f"NOT WRITTEN: {sum(map(len, absent.values()))} files or md5 files absent; samples mixing layouts {mixed}")
lines = []
for s, rs in samples.items():
    if s in absent:
        continue
    mates = len(rs[0][1])
    size = sum(os.path.getsize(os.path.join(r, "fastq", n)) for _, ns in rs for n in ns)
    r1 = ",".join(ns[0] for _, ns in rs)
    r2 = ",".join(ns[1] for _, ns in rs) if mates == 2 else "-"
    lines.append((size, s, mates, ",".join(a for a, _ in rs), r1, r2))
if not lines:
    sys.exit(f"NOT WRITTEN: no sample of {len(samples)} has every file")
lines.sort(key=lambda l: (-l[0], l[1]))
with open(out, "w") as f:
    f.write("sample\tmates\truns\tr1\tr2\tbytes\n")
    for size, s, mates, runs, r1, r2 in lines:
        f.write(f"{s}\t{mates}\t{runs}\t{r1}\t{r2}\t{size}\n")
print(f"{out}: {len(lines)} of {len(samples)} samples, {len(rows)} runs in the sheet, {sum(l[0] for l in lines)} bytes; first {lines[0][1]}, {lines[0][0]} bytes")
PY

for p in "$@"; do
    R="$S/$p/nfcore_fetchngs/results"
    mkdir -p "$E/$p"
    # fetchngs writes the samplesheet only once every run has downloaded (first run, 2026-10-22:
    # lung stopped this block with no samplesheet), so a project still downloading waits.
    if [[ ! -f "$S/$p/nfcore_fetchngs/results.complete" ]]; then
        echo "$p: its fetchngs run has not ended, so it has no samplesheet yet; pilot it once it has"
        continue
    fi
    T=tasks.tsv
    if [[ ! -s "$E/$p/$T" ]]; then python3 "$E/ed6_tasks.py" "$R" "$E/$p/$T"; fi
    echo "$T" > "$E/$p/pilot_tasks.txt"
    awk -F'\t' 'NR == 2 {print $1}' "$E/$p/$T" > "$E/$p/pilot.txt"
    echo "$p: pilot $(cat "$E/$p/pilot.txt") from $T"
done

dep=()
if [[ ! -s "$IDX/SA" ]]; then
    old=$(awk '$1 == "index" {print $2}' jobs.txt | tail -n 1)
    if [[ -n "$old" ]] && squeue -h -j "$old" 2>/dev/null | grep -q .; then
        echo "index job $old still queued or running"
    else
        old=$(sbatch --parsable -A "$A" --chdir="$E" --export=ALL,E="$E",REFS="$REFS",IDX="$IDX",BIN="$BIN" ed6_index.sbatch)
        echo "index $old" | tee -a jobs.txt
    fi
    dep=(--dependency=afterok:"$old")
fi
for p in "$@"; do
    # A project the loop above skipped has no pilot task list (2026-10-22: lung and pancreas got empty pilots).
    [[ -s "$E/$p/pilot_tasks.txt" ]] || continue
    if grep -q "^pilot_$p " jobs.txt; then echo "$p: a pilot was submitted before: $(grep "^pilot_$p " jobs.txt)"; continue; fi
    j=$(sbatch --parsable -A "$A" -c 8 --mem=56G -t 0-12:00:00 --array=1 "${dep[@]}" --chdir="$E/$p" \
        --export=ALL,E="$E",PROJ="$p",IDX="$IDX",GTF="$GTF",BIN="$BIN",CTX_REPO="$CTX_REPO",KEEP=1,PILOT=1,TASKS="$(cat "$E/$p/pilot_tasks.txt")" "$E/ed6_align.sbatch")
    echo "pilot_$p $j" | tee -a jobs.txt
done
squeue -u "$USER"
SCRIPT
bash ed6b_1_pilot.sh tier1 breast_tissue breast_mac breast_mono lung_gse229705 pancreas_gse192647 2>&1 | tee ed6b_1_pilot.log
```

Ready to submit now: the four complete projects pilot their largest sample. Lung and pancreas pilot after
their downloads end, with `bash ed6b_1_pilot.sh lung_gse229705 pancreas_gse192647`: a project still downloading is skipped,
since fetchngs writes its samplesheet last.

Paste back `ed6b_1_pilot.log`: the free space, the commit, `versions.txt`, `refs.txt`, the two MD5
lines "OK", each project's tasks line (samples written of the sheet's, bytes, the pilot sample), the job ids.

## Block B2 (ED.6b). Kebnekaise login node, once the pilots have ended: the measurements

```
cd /proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/ed6_align
cat > ed6b_2_report.sh <<'SCRIPT'
set -euo pipefail
E=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/ed6_align
cd "$E"
ids=$(awk '$1 == "index" || $1 ~ /^pilot_/ {print $2}' jobs.txt | paste -sd, -)
# Without -X, so each batch step's MaxRSS shows.
sacct -n -P -j "$ids" --format=JobID,JobName%30,Elapsed,ElapsedRaw,MaxRSS,ReqMem,AllocCPUS,TotalCPU,State
if [[ -f logs/star_index.time ]]; then grep -e 'Elapsed (wall' -e 'Maximum resident' logs/star_index.time; fi
. ./refs.txt
if [[ -s "$IDX/SA" ]]; then du -sb "$IDX" "$REFS/unzipped"; wc -l "$IDX/geneInfo.tab"; fi   # geneInfo.tab: expect 78942, a count line and 78,941 genes
for p in tier1 breast_tissue breast_mac breast_mono lung_gse229705 pancreas_gse192647; do
    [[ -s "$p/pilot.txt" ]] || continue
    sample=$(cat "$p/pilot.txt")
    D="$p/runs/$sample"
    echo "== $p $sample"
    if [[ ! -d "$D" ]]; then
        echo "$p: no $D: the pilot has not ended or has failed"
        ls -d "$p"/runs/* 2>/dev/null || true
        tail -n 30 "$p"/slurm-ed6_align-*.out 2>/dev/null || true
        continue
    fi
    grep -H -e 'Elapsed (wall' -e 'Maximum resident' "$D"/time/*.time | sed "s|$D/time/||"
    cat "$D/check/fastp_equal.txt"
    python3 - "$p" "$D" "$sample" "$p/$(cat "$p/pilot_tasks.txt")" <<'PY'
import csv, glob, json, os, sys
p, d, s, tasks = sys.argv[1:]
task = next(t for t in csv.DictReader(open(tasks), delimiter="\t") if t["sample"] == s)
log = {k.strip(): v.strip() for k, _, v in (l.partition("|") for l in open(f"{d}/star/{s}.Log.final.out")) if _}
inp = int(log["Number of input reads"])
passed = sum(json.load(open(f))["filtering_result"]["passed_filter_reads"] for f in glob.glob(f"{d}/fastp/*.fastp.json"))
pairs = passed // 2 if task["mates"] == "2" else passed
cols = [0, 0, 0]
for line in open(f"{d}/star/{s}.ReadsPerGene.out.tab"):
    f = line.rstrip("\n").split("\t")
    cols = [c + int(x) for c, x in zip(cols, f[1:4])]
assigned = {}
for f in sorted(glob.glob(f"{d}/featurecounts/*.featureCounts.txt.summary")):
    rows = dict(l.rstrip("\n").split("\t")[:2] for l in open(f) if not l.startswith("Status"))
    assigned[os.path.basename(f).split(".")[0]] = int(rows["Assigned"])
samples = [int(l) for l in open(f"{d}/time/du_samples.txt") if l.strip()]
before = {l.split("\t")[1].strip(): int(l.split("\t")[0]) for l in open(f"{d}/time/du_before_delete.txt")}
kept = int(open(f"{d}/time/du_kept.txt").read().split()[0])
print(f"{p}\t{s}\tfastq_bytes {task['bytes']}\tSTAR input {inp}\tfastp passed (reads or pairs) {pairs}\tI1 {inp == pairs}")
print(f"{p}\tReadsPerGene column sums {cols}, each equal to the input: {all(c == inp for c in cols)}")
print(f"{p}\tunique {log['Uniquely mapped reads %']}, multi {log['% of reads mapped to multiple loci']}")
print(f"{p}\tAssigned {assigned}")
print(f"{p}\tpeak task bytes {max(samples) if samples else 'none sampled'} over {len(samples)} samples a minute; before delete {before}; kept {kept}")
PY
    ls -la "$D/star" "$D/featurecounts" | head -n 60
done
df -h /proj/nobackup/hpc2n2024-150
SCRIPT
bash ed6b_2_report.sh 2>&1 | tee ed6b_2_report.log
```

**Paste back `ed6b_2_report.log` whole.** What it must show, per project:
- the sacct rows: the batch step's `ElapsedRaw`, `MaxRSS`, `AllocCPUS` and `TotalCPU`, `State` COMPLETED;
- per step (fastp, FastQC, seqkit, STAR, samtools, the 12 featureCounts runs) the wall time and max RSS;
- `fastp JSON equal with and without -o, command aside: True` for every run;
- `I1 True`: STAR's input equals fastp's passed reads (pairs for a paired sample);
- each ReadsPerGene column summing to STAR's input;
- Assigned for all 12 settings;
- the peak task bytes, the bytes before the BAM and trimmed reads would be deleted, and the kept bytes;
- for the index, if it was built: its wall time, max RSS, size, and 78942 lines in `geneInfo.tab`.

**Stop** on any False, a non-COMPLETED state, or a time near 12 h: the limits are then wrong, and
the full runs wait.

## Block C1 (ED.6c). After the pilot's numbers: the full arrays

Written now, run only after block B2's log is read. Every resource is a variable set from that
log; the script refuses to run with one unset. Rules for setting them (choices, not measurements;
change them if the pilot says otherwise):

| variable | rule |
|---|---|
| `C` | the pilot's 8, unless `TotalCPU / ElapsedRaw` shows the tasks use far fewer cores |
| `MEM` | the largest batch-step MaxRSS of the five pilots times 1.25, rounded up to whole G |
| `TIME[<slug>]` | that project's pilot ElapsedRaw times 1.5, rounded up, written as sbatch's `D-HH:MM:SS`; the pilot is the project's largest sample |
| `N` | concurrent tasks: at most floor((free bytes - kept_total - 200 GB margin) / peak task bytes), and no more than the queue grants |

The arrays run one project after another (`afterany` on the previous array), so at most `N`
tasks run at once and the disk bound is one formula. The pilot's sample is already done, and its
task exits at once wherever it sits in the list.

```
cd /proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/ed6_align
cat > ed6c_1_full.sh <<'SCRIPT'
set -euo pipefail
S=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial
E="$S/ed6_align"
A=hpc2n2026-161
# Set from ed6b_2_report.log (ED.6c's rules) on 2026-10-22; N=20 is a choice, the disk allowed 167.
C=8
MEM=56G
N=20
declare -A TIME=([tier1]=0-01:32:00 [breast_tissue]=0-00:27:00 [breast_mac]=0-01:22:00 [breast_mono]=0-00:58:00 [lung_gse229705]= [pancreas_gse192647]=)
[[ -n "$C" && -n "$MEM" && -n "$N" ]] || { echo "set C, MEM and N from the pilot first"; exit 1; }
cd "$E"
BIN=$(cat bin.txt)
. ./refs.txt
CTX_REPO="$S/RNAseq_tutorial"
[[ -s "$IDX/SA" ]] || { echo "no index at $IDX"; exit 1; }
# The pilots ran an align script that kept the 133 MB matrices; block B1 run again rewrites it
# with fc_compact and submits nothing new (its pilots are in jobs.txt).
grep -q fc_compact.py ed6_align.sbatch || { echo "ed6_align.sbatch predates fc_compact: run block B1 again first"; exit 1; }
for p in "$@"; do
    s=$(cat "$p/pilot.txt")
    if compgen -G "$p/runs/$s/featurecounts/*.featureCounts.txt" > /dev/null; then
        python3 "$CTX_REPO/src/ctx/process/fc_compact.py" "$p/featurecounts_annotation.tsv.gz" "$p/runs/$s/featurecounts/"*.featureCounts.txt
    fi
done
df -h /proj/nobackup/hpc2n2024-150
prev=()
for p in "$@"; do
    [[ -n "${TIME[$p]:-}" ]] || { echo "set TIME[$p] from its pilot first"; exit 1; }
    [[ -d "$p/runs/$(cat "$p/pilot.txt")" ]] || { echo "$p: its pilot has not passed"; exit 1; }
    if [[ ! -s "$p/tasks.tsv" ]]; then
        # Lung or pancreas, piloted mid-download: the full list once fetchngs has written its marker.
        [[ -f "$S/$p/nfcore_fetchngs/results.complete" ]] || { echo "$p: download not complete"; exit 1; }
        python3 ed6_tasks.py "$S/$p/nfcore_fetchngs/results" "$p/tasks.tsv"
    fi
    if grep -q "^full_$p " jobs.txt; then echo "$p: submitted before: $(grep "^full_$p " jobs.txt)"; continue; fi
    n=$(($(wc -l < "$p/tasks.tsv") - 1))
    j=$(sbatch --parsable -A "$A" -c "$C" --mem="$MEM" -t "${TIME[$p]}" --array=1-"$n"%"$N" "${prev[@]}" --chdir="$E/$p" \
        --export=ALL,E="$E",PROJ="$p",IDX="$IDX",GTF="$GTF",BIN="$BIN",CTX_REPO="$CTX_REPO",KEEP=0,PILOT=0 "$E/ed6_align.sbatch")
    echo "full_$p $j $n" | tee -a jobs.txt
    prev=(--dependency=afterany:"$j")
done
squeue -u "$USER"
SCRIPT
bash ed6c_1_full.sh tier1 breast_tissue breast_mac breast_mono lung_gse229705 pancreas_gse192647 2>&1 | tee ed6c_1_full.log
```

**M2's fastp grid**, one job per project, on the run `ctx qc --reads` grids by default: the first
run of the dataset's gather sheet. Read it on the Mac first, for each project, and set `GRID`:

```
python3 -c 'import csv, glob, sys; f = sorted(glob.glob(sys.argv[1]))[-1]; print(f, next(csv.DictReader(open(f)))["run_accession"])' "$HOME/ctx_tier1/outputs/01_gather/*_gather_tier1_samplesheet_v*.csv"
```

Run on the Mac 2026-10-22, it printed ERR7425257 (tier1, `2026-10-20_gather_tier1_samplesheet_v1.csv`),
SRR5811633 (breast_tissue), SRR7631637 (breast_mac) and SRR7631576 (breast_mono), each from its
2026-10-22 sheet in `~/ctx_breast`; they are filled in below. Lung's waits for `~/ctx_lung_gse229705`
(ED.4b block 6), pancreas's for `~/ctx_pancreas_gse192647`. ED.8a's metadata-only gather sheet,
in a scratch root, put SRR17335068 first; it is not filled in, because the rule is the real root's
sheet. Then on Kebnekaise:

```
cd /proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/ed6_align
cat > ed6c_2_grid.sh <<'SCRIPT'
set -euo pipefail
E=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/ed6_align
A=hpc2n2026-161
# The first run of each gather sheet, printed on the Mac 2026-10-22; lung's and pancreas's once their roots exist.
declare -A GRID=([tier1]=ERR7425257 [breast_tissue]=SRR5811633 [breast_mac]=SRR7631637 [breast_mono]=SRR7631576 [lung_gse229705]= [pancreas_gse192647]=)
cd "$E"
BIN=$(cat bin.txt)
cat > ed6_grid.sbatch <<'EOF'
#!/usr/bin/env bash
#SBATCH -J ed6_grid
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 0-06:00:00
#SBATCH -o slurm-%x-%j.out
# ctx.qc.reads.fastp_grid: fastp without -o at every GRID_Q x GRID_L, -u 40 fixed.
set -euo pipefail
: "${E:?}" "${PROJ:?}" "${RUN:?}" "${BIN:?}"
FQ="/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/$PROJ/nfcore_fetchngs/results/fastq"
G="$E/$PROJ/grid"
[ ! -e "$G" ] || { echo "$G exists: move it aside first"; exit 1; }
mkdir -p "$G.part"
row=$(awk -F'\t' -v r="$RUN" 'NR > 1 { n = split($3, a, ","); split($4, b, ","); split($5, c, ","); for (i = 1; i <= n; i++) if (a[i] == r) print $2 "\t" b[i] "\t" c[i] }' "$E/$PROJ/tasks.tsv")
[ -n "$row" ] || { echo "$RUN is in no task of $PROJ"; exit 1; }
IFS=$'\t' read -r mates f1 f2 <<< "$row"
given=(-i "$FQ/$f1")
if [ "$mates" = 2 ]; then given+=(-I "$FQ/$f2"); fi
for q in 10 15 20 25 30 35; do
    for l in 15 25 36 50 75 102; do
        t="grid_q${q}_l${l}"
        "$BIN/fastp" "${given[@]}" -R "$t" -j "$G.part/$t.fastp.json" -h /dev/null -w "$SLURM_CPUS_PER_TASK" -q "$q" -u 40 -l "$l" > "$G.part/$t.fastp.log" 2>&1
    done
done
mv "$G.part" "$G"
EOF
for p in "$@"; do
    [[ -n "${GRID[$p]:-}" ]] || { echo "set GRID[$p] from the Mac first"; exit 1; }
    j=$(sbatch --parsable -A "$A" --chdir="$E/$p" --export=ALL,E="$E",PROJ="$p",RUN="${GRID[$p]}",BIN="$BIN" ed6_grid.sbatch)
    echo "grid_$p $j ${GRID[$p]}" | tee -a jobs.txt
done
SCRIPT
bash ed6c_2_grid.sh tier1 breast_tissue breast_mac breast_mono lung_gse229705 pancreas_gse192647 2>&1 | tee ed6c_2_grid.log
```

The grid's 4 cores, 8G and 6 h are estimates: 36 fastp passes over one run, unmeasured.

**After the arrays: the check.**

```
cd /proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/ed6_align
cat > ed6c_3_check.sh <<'SCRIPT'
set -euo pipefail
E=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/ed6_align
cd "$E"
ids=$(awk '$1 ~ /^full_|^grid_/ {print $2}' jobs.txt | paste -sd, -)
sacct -n -P -j "$ids" --format=JobID,ElapsedRaw,MaxRSS,AllocCPUS,State | python3 -c '
import sys, collections
states, cpu, rss = collections.Counter(), 0, 0
for l in sys.stdin:
    j, el, m, c, st = l.rstrip("\n").split("|")
    if ".batch" in j and m:
        unit = {"K": 1e3, "M": 1e6, "G": 1e9}.get(m[-1], 1)
        rss = max(rss, float(m[:-1] if m[-1].isalpha() else m) * unit)
    elif "." not in j:
        states[st.split()[0]] += 1; cpu += int(el or 0) * int(c or 0)
print("tasks by state", dict(states), "; core-hours", round(cpu / 3600, 1), "; largest MaxRSS bytes", int(rss))'
for p in tier1 breast_tissue breast_mac breast_mono lung_gse229705 pancreas_gse192647; do
    [[ -s "$p/tasks.tsv" ]] || continue
    want=$(($(wc -l < "$p/tasks.tsv") - 1))
    done_=$(find "$p/runs" -mindepth 1 -maxdepth 1 -type d ! -name '*.part*' | wc -l)
    printf '%s: %s of %s samples done; failed tries: %s; grid JSONs: %s; bytes %s\n' "$p" "$done_" "$want" \
        "$(find "$p/runs" -maxdepth 1 -name '*.part*' | wc -l)" "$(find "$p/grid" -name '*.fastp.json' 2>/dev/null | wc -l)" "$(du -sb "$p" | cut -f1)"
    python3 - "$p" <<'PY'
import glob, os, sys
p = sys.argv[1]
bad = []
for d in sorted(glob.glob(f"{p}/runs/*")):
    if ".part" in d:
        continue
    s = os.path.basename(d)
    log = {k.strip(): v.strip() for k, _, v in (l.partition("|") for l in open(f"{d}/star/{s}.Log.final.out")) if _}
    n = int(log["Number of input reads"])
    cols = [0, 0, 0]
    for line in open(f"{d}/star/{s}.ReadsPerGene.out.tab"):
        cols = [c + int(x) for c, x in zip(cols, line.split("\t")[1:4])]
    if any(c != n for c in cols) or len(glob.glob(f"{d}/featurecounts/*.summary")) != 12 or glob.glob(f"{d}/star/*.bam"):
        bad.append(s)
print(f"{p}: samples failing a ReadsPerGene sum, lacking 12 featureCounts summaries, or still holding a BAM: {len(bad)} {bad[:10]}")
PY
done
df -h /proj/nobackup/hpc2n2024-150
SCRIPT
bash ed6c_3_check.sh 2>&1 | tee ed6c_3_check.log
```

Paste back `ed6c_3_check.log`. The pilot samples (`<slug>/pilot.txt`) are expected under "still holding a BAM" until
the last line of ED.6d. A failed try (`runs/<sample>.part.failed_<job>`) is resubmitted by the same
`sbatch` line of block C1 with `--array=<its task line>`, after Rytis has read its log and removed it.

**Expected cost, as formulas** (all estimates until ed6c_3_check.log measures them; `E_p` the pilot's
ElapsedRaw in s, `b_p` the pilot sample's FASTQ bytes, `B_p` the project's FASTQ bytes from the
table above, `n_p` its samples, `k_p` the pilot's kept bytes, `peak` the largest pilot peak task bytes):

- core-hours of the arrays: `C x sum_p (E_p x B_p / b_p) / 3600`. It assumes time linear in bytes; the
  index load is a fixed cost per task that this undercounts for the small samples;
- core-hours of the grids: 4 x each grid job's ElapsedRaw / 3600;
- kept on the cluster: `kept_total = sum_p n_p x k_p`, plus the pilots' BAMs and trimmed reads until removed;
- peak disk while running: `kept_total + N x peak`;
- copied to the Mac (ED.6d): `kept_total` less the FastQC HTML, fastp HTML, STAR's `Log.out`,
  `Log.progress.out` and `SJ.out.tab` and the logs; measured by the copy itself.

## Block C4 (ED.6c). Run 2026-10-22 16:17: lung and pancreas unattended, block E for four

Block B1's rerun at 16:10 skipped lung and pancreas as still downloading, but its submit loop sent
each a pilot with no task list (48051523, 48051524). This block cancelled them, moved their lines
to `jobs_cancelled.txt`, patched the cluster's `ed6b_1_pilot.sh` as block B1 above now reads, and
ran block E with the values above. For lung and pancreas it queued one 1-core job per project
after its download (48041315, 48047926): B1 for that project, then, after the pilot, C1 with
`TIME` from the pilot's sacct by C1's rule (1.5 x ElapsedRaw, rounded up), refusing if the
pilot's MaxRSS x 1.25 exceeds 56G. Their grids and block E wait for their Mac roots.

```
cd /proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/ed6_align
cat > ed6_overnight.sh <<'NIGHT'
set -euo pipefail
S=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial
E="$S/ed6_align"
A=hpc2n2026-161
D="$S/RNAseq_tutorial/docs/plans/2026-10-22_ed6_m2_m3_core_cluster_blocks.md"
cd "$E"
command -v flock > /dev/null || { echo "no flock on this node"; exit 1; }
# 1. B1 skipped lung and pancreas as still downloading, but its submit loop still sent each a
# pilot with no task list. Cancel them and move their lines aside, so B1 sends them again.
scancel 48051523 48051524 || true
cp jobs.txt jobs.txt.before_overnight
grep -E '^pilot_(lung_gse229705|pancreas_gse192647) ' jobs.txt >> jobs_cancelled.txt || true
grep -vE '^pilot_(lung_gse229705|pancreas_gse192647) ' jobs.txt > jobs.txt.new
mv jobs.txt.new jobs.txt
for p in lung_gse229705 pancreas_gse192647; do
    # A task list written mid-download would be reused by B1 and C1 as if it were complete.
    if [[ -e "$p/tasks.tsv" ]]; then mv "$p/tasks.tsv" "$p/tasks.tsv.before_download_$(date +%Y%m%d%H%M)"; fi
done
# 2. B1 fixed: the submit loop skips a project with no pilot task list; pilots get C1's 56G.
python3 - ed6b_1_pilot.sh <<'PY'
import sys
path = sys.argv[1]
text = open(path).read()
loop = '    if grep -q "^pilot_$p " jobs.txt; then'
guard = '    [[ -s "$E/$p/pilot_tasks.txt" ]] || continue\n'
mem_old, mem_new = "-c 8 --mem=48G -t 0-12:00:00", "-c 8 --mem=56G -t 0-12:00:00"
if guard + loop not in text:
    assert text.count(loop) == 1 and text.count(mem_old) == 1, "ed6b_1_pilot.sh is not the version this patch expects"
    text = text.replace(loop, guard + loop).replace(mem_old, mem_new)
    open(path, "w").write(text)
print("ed6b_1_pilot.sh: skip guard and 56G in place")
PY
# 3. Block E for the four projects: C1's C and MEM, TIME 3 x each pilot's ElapsedRaw.
awk -v s="cat > ed6e_featured.sh <<'SCRIPT'" '$0 == s {on=1} on {print} on && /^SCRIPT$/ {exit}' "$D" > write_ed6e_featured.sh
bash write_ed6e_featured.sh
sed -i -e 's/^C=$/C=8/' -e 's/^MEM=$/MEM=56G/' \
    -e 's/\[tier1\]= \[breast_tissue\]= \[breast_mac\]= \[breast_mono\]=/[tier1]=0-03:03:00 [breast_tissue]=0-00:54:00 [breast_mac]=0-02:44:00 [breast_mono]=0-01:56:00/' ed6e_featured.sh
grep -nE '^(C|MEM)=|declare -A TIME' ed6e_featured.sh
bash ed6e_featured.sh tier1 breast_tissue breast_mac breast_mono 2>&1 | tee ed6e_featured.log
# 4. Lung and pancreas unattended: after the download, B1 (task list, pilot); after the pilot,
# C1 with TIME from the pilot's sacct by C1's rule. Each step is a 1-core job.
cat > ed6_after_download.sh <<'EOF'
set -euo pipefail
p=$1
S=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial
E="$S/ed6_align"
A=hpc2n2026-161
cd "$E"
for v in ${!SLURM_@}; do unset "$v"; done
[[ -f "$S/$p/nfcore_fetchngs/results.complete" ]] || { echo "$p: no results.complete"; exit 1; }
flock "$E/b1.lock" bash ed6b_1_pilot.sh "$p" 2>&1 | tee "logs/ed6b_1_pilot_$p.log"
pid=$(awk -v k="pilot_$p" '$1 == k {print $2}' jobs.txt | tail -n 1)
[[ -n "$pid" ]] || { echo "$p: B1 submitted no pilot"; exit 1; }
j=$(sbatch --parsable -A "$A" -c 1 --mem=4G -t 0-01:00:00 -J ed6_full_plan --dependency=afterok:"$pid" \
    --chdir="$E" -o "$E/logs/ed6_full_plan_${p}_%j.out" --wrap "bash $E/ed6_full_after_pilot.sh $p")
echo "plan_$p $j" | tee -a jobs.txt
EOF
cat > ed6_full_after_pilot.sh <<'EOF'
set -euo pipefail
p=$1
E=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/ed6_align
cd "$E"
for v in ${!SLURM_@}; do unset "$v"; done
pid=$(awk -v k="pilot_$p" '$1 == k {print $2}' jobs.txt | tail -n 1)
el=$(sacct -n -X -P -j "$pid" --format=ElapsedRaw | head -n 1)
rss=$(sacct -n -P -j "$pid" --format=MaxRSS | sed -n 's/K$//p' | sort -n | tail -n 1)
[[ -n "$el" && -n "$rss" ]] || { echo "$p: pilot $pid has no ElapsedRaw or no MaxRSS in K"; exit 1; }
# C1's MEM rule over this pilot too: its MaxRSS x 1.25 must fit C1's 56G.
(( rss * 5 / 4 <= 56 * 1024 * 1024 )) || { echo "$p: pilot MaxRSS ${rss}K x 1.25 is over 56G; set MEM by hand"; exit 1; }
m=$(( (el * 3 + 119) / 120 ))
T=$(printf '%d-%02d:%02d:00' $((m / 1440)) $((m % 1440 / 60)) $((m % 60)))
echo "$p: pilot $pid ElapsedRaw $el s, MaxRSS ${rss}K; TIME $T (1.5 x ElapsedRaw, rounded up)"
sed "s/\[$p\]=\([ )]\)/[$p]=$T\1/" ed6c_1_full.sh > "ed6c_1_full_$p.sh"
grep -n 'declare -A TIME' "ed6c_1_full_$p.sh"
flock "$E/b1.lock" bash "ed6c_1_full_$p.sh" "$p" 2>&1 | tee "logs/ed6c_1_full_$p.log"
EOF
declare -A DL=([lung_gse229705]=48041315 [pancreas_gse192647]=48047926)
for p in lung_gse229705 pancreas_gse192647; do
    dep=()
    if squeue -h -j "${DL[$p]}" 2>/dev/null | grep -q .; then
        dep=(--dependency=afterok:"${DL[$p]}")
    elif [[ ! -f "$S/$p/nfcore_fetchngs/results.complete" ]]; then
        echo "$p: download job ${DL[$p]} has left the queue without results.complete; read its log"; continue
    fi
    j=$(sbatch --parsable -A "$A" -c 1 --mem=4G -t 0-02:00:00 -J ed6_after_dl "${dep[@]}" --chdir="$E" \
        -o "$E/logs/ed6_after_dl_${p}_%j.out" --wrap "bash $E/ed6_after_download.sh $p")
    echo "after_dl_$p $j" | tee -a jobs.txt
done
squeue -u "$USER"
NIGHT
bash ed6_overnight.sh 2>&1 | tee ed6_overnight.log
```

## Block E (ED.6e-6). The featured sample's settings and the references, after C1

Rytis decided on 2026-10-22: the featured sample's M3 teaching variants run on the cluster too,
for each project. The local route (`ctx.process.run.run` and `align_all`) runs two of them on the
featured sample, and block E runs both exactly as it does:

- `cut`: `reference.prepare(genome, gtf, parent, half=True)` (each contig's first floor(length / 2)
  bases, the genes that end inside), `reference.star_index`'s genomeGenerate on that half
  (`--runThreadN <n> --genomeDir <dir> --genomeFastaFiles genome.fa --sjdbGTFfile annotation.gtf
  --sjdbOverhang 100 --genomeSAindexNbases <sa_index_bases(half length)>`, nothing else), STAR with
  `STAR_OPTIONS` on the trimmed reads, `align.bam_stats`, then featureCounts with the half's
  `annotation.gtf` as `count_all` counts it, `-s <k> [-p --countReadPairs]`, at all three strands,
  since the Mac calls the strand;
- `untrimmed`: STAR with `STAR_OPTIONS` on the raw FASTQs (`Trimmed.raw` without a cap), and
  `align.bam_stats`.

`--multimap-max` is not here. The local route applies it to every sample's main alignment and to
both featured alignments (`AlignOptions.multimap`), so on the cluster it would be block C1 again
with `--outFilterMultimapNmax` added, for every sample; one featured run cannot stand for it.
`ctx process --align-dir` refuses it, as it refuses any setting the cluster did not run. Whether
to run that second array is Rytis's call.

ctx prepares both references here, the whole one too: its `reference.json` gives the Mac page the
genome's facts and start codons (D-039's check) without the 3.1 GB genome, and the Mac rebuilds the
annotation lines with the same code and checks their SHA-256. ctx needs its analysis environment for
this, so the block installs `laptop` from the committed lock beside `reads`.

The featured sample is the Mac's default `--featured-sample`: the sample of the gather sheet's first
run, printed by the same one-liner as the grid's run (block C1) with `["sample"]` in place of
`["run_accession"]`. Read on the Mac 2026-10-22 from the sheets named in block C1: ERX6994339
(tier1), SRX2990193 (breast_tissue), SRX4495392 (breast_mac), SRX4495339 (breast_mono); lung's once
its root exists.

Resources, none measured. The reference job: 8 cores, 48G, 4 h, estimates (M3.11's index guess,
the half smaller; ctx's Python passes over the 3.1 GB genome unmeasured). Each featured job takes
C1's `C` and `MEM` and `3 x` the pilot's ElapsedRaw (two alignments where the pilot ran one, with
C1's 1.5 margin), a choice like C1's; the pilot is the project's largest sample, so it bounds the
featured one.

```
cd /proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/ed6_align
cat > ed6e_featured.sh <<'SCRIPT'
set -euo pipefail
S=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial
E="$S/ed6_align"
A=hpc2n2026-161
CTX_REPO="$S/RNAseq_tutorial"
PIXI="${PIXI_HOME:-$HOME/.pixi}/bin/pixi"
# C1's C and MEM; TIME is 3 x each pilot's ElapsedRaw, as sbatch's D-HH:MM:SS. None is guessed.
C=8
MEM=56G
declare -A TIME=([tier1]=0-03:03:00 [breast_tissue]=0-00:54:00 [breast_mac]=0-02:44:00 [breast_mono]=0-01:56:00 [lung_gse229705]=)
# The sample of each gather sheet's first run, read on the Mac 2026-10-22; lung's once its root exists.
declare -A FEATURED=([tier1]=ERX6994339 [breast_tissue]=SRX2990193 [breast_mac]=SRX4495392 [breast_mono]=SRX4495339 [lung_gse229705]=)
[[ $# -ge 1 ]] || { echo "usage: bash ed6e_featured.sh <slug> [<slug> ...]"; exit 1; }
[[ -n "$C" && -n "$MEM" ]] || { echo "set C and MEM from the pilot first"; exit 1; }
cd "$E" && mkdir -p logs reference && touch jobs.txt
BIN=$(cat bin.txt)
. ./refs.txt
git -C "$CTX_REPO" log -1 --format='%h %s'
PIXI_CACHE_DIR="$S/pixi_cache" "$PIXI" install --manifest-path "$CTX_REPO/pixi.toml" -e laptop --frozen
PY="$("$PIXI" run --manifest-path "$CTX_REPO/pixi.toml" -e laptop --frozen printenv CONDA_PREFIX)/bin/python"
"$PY" -c 'import ctx.process.reference as r; print(r.__file__)'

cat > ed6_reference.sbatch <<'EOF'
#!/usr/bin/env bash
#SBATCH -J ed6_reference
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 0-04:00:00
#SBATCH -o slurm-%x-%j.out
# ctx's reference.prepare, whole and half, on the MD5-checked GENCODE 50 files; then
# reference.star_index's genomeGenerate on the half, its options and nothing else.
set -euo pipefail
: "${E:?}" "${REFS:?}" "${BIN:?}" "${PY:?}"
FA="$REFS/GRCh38.primary_assembly.genome.fa.gz"
GZ="$REFS/gencode.v50.primary_assembly.annotation.gtf.gz"
(cd "$REFS" && printf '%s  %s\n' da1a11258be075cfa7af718162c894e7 "$(basename "$FA")" 289b91e5e95e8b0450d223246f10a12e "$(basename "$GZ")" | md5sum -c -)
/usr/bin/time -v -o "$E/logs/reference_prepare.time" "$PY" - "$FA" "$GZ" "$E/reference" <<'PY'
import shutil, sys
from pathlib import Path
from ctx.process import reference
fa, gtf, out = map(Path, sys.argv[1:])
for name, half in (("main", False), ("cut", True)):
    ref = reference.prepare(fa, gtf, out / "prepared", half=half)
    (out / name).mkdir(exist_ok=True)
    for f in ("reference.json", "ctx_manifest.json"):
        shutil.copyfile(ref.folder / f, out / name / f)
    bases = reference.sa_index_bases(ref.genome.length)
    print(name, ref.folder, ref.genome.length, "bases", bases, "ATG share", ref.start_codons.share)
    if half:
        (out / "half.txt").write_text(f"FOLDER={ref.folder}\nBASES={bases}\n")
PY
. "$E/reference/half.txt"
IDX_CUT="$E/reference/star_2.7.11b_gencode50_half"
if [ ! -s "$IDX_CUT/SA" ]; then
    [ ! -e "$IDX_CUT.part" ] || { echo "$IDX_CUT.part exists, from a failed build: move it aside first"; exit 1; }
    mkdir -p "$IDX_CUT.part"
    /usr/bin/time -v -o "$E/logs/star_index_cut.time" "$BIN/STAR" --runMode genomeGenerate --runThreadN 8 \
        --genomeDir "$IDX_CUT.part" --genomeFastaFiles "$FOLDER/genome.fa" --sjdbGTFfile "$FOLDER/annotation.gtf" \
        --sjdbOverhang 100 --genomeSAindexNbases "$BASES" --outFileNamePrefix "$IDX_CUT.part/"
    mv "$IDX_CUT.part" "$IDX_CUT"
fi
printf 'GENOME=%s\nCUT_INDEX=%s\nCUT_GTF=%s\n' "$FA" "$IDX_CUT" "$FOLDER/annotation.gtf" > "$E/reference/inputs.txt"
du -sb "$E/reference"/* "$IDX_CUT"
EOF

cat > ed6_featured.sbatch <<'EOF'
#!/usr/bin/env bash
#SBATCH -J ed6_featured
#SBATCH -n 1
#SBATCH -o slurm-%x-%j.out
# The featured sample: fastp again for the trimmed reads (its JSON must equal C1's, command aside),
# STAR untrimmed on the raw reads, STAR on the half, samtools, featureCounts on the half at -s 0, 1, 2.
set -euo pipefail
: "${E:?}" "${PROJ:?}" "${SAMPLE:?}" "${IDX:?}" "${BIN:?}" "${KEEP:?}" "${CTX_REPO:?}"
. "$E/reference/inputs.txt"
T="${SLURM_CPUS_PER_TASK:?submit with sbatch -c}"
FQ="/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/$PROJ/nfcore_fetchngs/results/fastq"
OUT="$E/$PROJ"
line=$(awk -F'\t' -v s="$SAMPLE" 'NR > 1 && $1 == s' "$OUT/tasks.tsv")
[ -n "$line" ] || { echo "$SAMPLE is in no task of $OUT/tasks.tsv"; exit 1; }
IFS=$'\t' read -r sample mates runs r1 r2 bytes <<< "$line"
[ -d "$OUT/runs/$sample" ] || { echo "$sample: C1 has not finished it"; exit 1; }
D="$OUT/featured/$sample"
P="$D.part"
if [ -e "$D" ]; then echo "$sample: $D exists; nothing run"; exit 0; fi
if [ -e "$P" ]; then mv "$P" "$P.failed_${SLURM_JOB_ID}"; fi
mkdir -p "$P"/{fastp,trimmed,check} "$P"/untrimmed/{star,time} "$P"/cut/{star,featurecounts,time}
IFS=, read -ra RUN <<< "$runs"
IFS=, read -ra F1 <<< "$r1"
F2=()
if [ "$mates" = 2 ]; then IFS=, read -ra F2 <<< "$r2"; fi
(cd "$FQ" && for f in "${F1[@]}" "${F2[@]}"; do md5sum -c "md5/$f.md5"; done) | tee "$P/check/md5.txt"
t1=(); t2=(); u1=(); u2=()
for i in "${!RUN[@]}"; do
    run=${RUN[$i]}
    given=(-i "$FQ/${F1[$i]}"); u1+=("$FQ/${F1[$i]}")
    if [ "$mates" = 2 ]; then
        given+=(-I "$FQ/${F2[$i]}"); u2+=("$FQ/${F2[$i]}")
        trimmed=(-o "$P/trimmed/${run}_1.fastq.gz" -O "$P/trimmed/${run}_2.fastq.gz")
        t1+=("$P/trimmed/${run}_1.fastq.gz"); t2+=("$P/trimmed/${run}_2.fastq.gz")
    else
        trimmed=(-o "$P/trimmed/$run.fastq.gz"); t1+=("$P/trimmed/$run.fastq.gz")
    fi
    /usr/bin/time -v -o "$P/cut/time/fastp_$run.time" "$BIN/fastp" "${given[@]}" -R "$run" -j "$P/fastp/$run.fastp.json" \
        -h /dev/null -w "$T" -q 15 -u 40 -l 15 "${trimmed[@]}" > "$P/fastp/$run.fastp.log" 2>&1
    python3 - "$OUT/runs/$sample/fastp/$run.fastp.json" "$P/fastp/$run.fastp.json" <<'PY' | tee -a "$P/check/fastp_equal.txt"
import json, sys
a, b = (json.load(open(p)) for p in sys.argv[1:])
a.pop("command", None); b.pop("command", None)
diff = sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))
print(f"{sys.argv[2]}: fastp JSON equal to C1's, command aside: {not diff}; differing keys {diff}")
sys.exit(1 if diff else 0)
PY
done
join() { local IFS=,; echo "$*"; }
star() {
    local dir=$1 index=$2; shift 2
    (cd "$dir" && /usr/bin/time -v -o "$dir/../time/star.time" "$BIN/STAR" --runThreadN "$T" --genomeDir "$index" --readFilesIn "$@" \
        --readFilesCommand gzip -dcf --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
        --outFileNamePrefix "$dir/$sample." > "$dir/$sample.star.log" 2>&1)
    for n in Log.final.out ReadsPerGene.out.tab Aligned.sortedByCoord.out.bam; do
        test -s "$dir/$sample.$n" || { echo "STAR exited 0 and wrote no $dir/$sample.$n"; exit 1; }
    done
    mv "$dir/$sample.Aligned.sortedByCoord.out.bam" "$dir/$sample.bam"
    /usr/bin/time -v -o "$dir/../time/samtools_index.time" "$BIN/samtools" index "$dir/$sample.bam"
    "$BIN/samtools" flagstat "$dir/$sample.bam" > "$dir/$sample.flagstat"
    "$BIN/samtools" idxstats "$dir/$sample.bam" > "$dir/$sample.idxstats"
}
raw=("$(join "${u1[@]}")"); cut=("$(join "${t1[@]}")")
if [ "$mates" = 2 ]; then raw+=("$(join "${u2[@]}")"); cut+=("$(join "${t2[@]}")"); fi
star "$P/untrimmed/star" "$IDX" "${raw[@]}"
star "$P/cut/star" "$CUT_INDEX" "${cut[@]}"
pe=()
if [ "$mates" = 2 ]; then pe=(-p --countReadPairs); fi
for s in 0 1 2; do
    (cd "$P/cut/star" && /usr/bin/time -v -o "$P/cut/time/fc_cut_s$s.time" "$BIN/featureCounts" -T "$T" -a "$CUT_GTF" -s "$s" "${pe[@]}" \
        -o "$P/cut/featurecounts/cut_s$s.featureCounts.txt" "$sample.bam" > "$P/cut/featurecounts/cut_s$s.featureCounts.log" 2>&1)
done
# The half's annotation is another one, so its columns are kept apart from the project's.
python3 "$CTX_REPO/src/ctx/process/fc_compact.py" "$P/cut/featurecounts_annotation.tsv.gz" \
    "$P"/cut/featurecounts/*.featureCounts.txt > "$P/cut/featurecounts/compact.log" 2>&1
if [ "$KEEP" = 0 ]; then rm -f "$P"/{untrimmed,cut}/star/"$sample".bam*; rm -rf "$P/trimmed"; fi
mv "$P" "$D"
echo "done: $sample"
EOF

ref=$(awk '$1 == "reference" {print $2}' jobs.txt | tail -n 1)
if [[ -s reference/inputs.txt ]]; then
    ref=
elif [[ -z "$ref" ]] || ! squeue -h -j "$ref" 2>/dev/null | grep -q .; then
    ref=$(sbatch --parsable -A "$A" --chdir="$E" --export=ALL,E="$E",REFS="$REFS",BIN="$BIN",PY="$PY" ed6_reference.sbatch)
    echo "reference $ref" | tee -a jobs.txt
fi
for p in "$@"; do
    [[ -n "${FEATURED[$p]:-}" && -n "${TIME[$p]:-}" ]] || { echo "set FEATURED[$p] and TIME[$p] first"; exit 1; }
    if grep -q "^featured_$p " jobs.txt; then echo "$p: submitted before: $(grep "^featured_$p " jobs.txt)"; continue; fi
    # After the references and C1's array for the project, whose fastp JSON the job compares.
    full=$(awk -v k="full_$p" '$1 == k {print $2}' jobs.txt | tail -n 1)
    conds=()
    if [[ -n "$ref" ]]; then conds+=("afterok:$ref"); fi
    if [[ -n "$full" ]]; then conds+=("afterany:$full"); fi
    after=()
    if [[ ${#conds[@]} -gt 0 ]]; then after=(--dependency="$(IFS=,; echo "${conds[*]}")"); fi
    j=$(sbatch --parsable -A "$A" -c "$C" --mem="$MEM" -t "${TIME[$p]}" "${after[@]}" --chdir="$E/$p" \
        --export=ALL,E="$E",PROJ="$p",SAMPLE="${FEATURED[$p]}",IDX="$IDX",BIN="$BIN",CTX_REPO="$CTX_REPO",KEEP=0 "$E/ed6_featured.sbatch")
    echo "featured_$p $j ${FEATURED[$p]}" | tee -a jobs.txt
done
squeue -u "$USER"
SCRIPT
bash ed6e_featured.sh tier1 breast_tissue breast_mac breast_mono lung_gse229705 2>&1 | tee ed6e_featured.log
```

Paste back `ed6e_featured.log`, then, once the jobs end, the reference job's printed lines (the two
references' lengths, SAindexNbases and ATG shares) from `slurm-ed6_reference-<job>.out`, each
`featured/<sample>/check/fastp_equal.txt`, and `logs/*.time` and `featured/*/*/time/*.time`'s
wall time and maximum resident set size lines. A featured job that finds a fastp JSON unlike C1's
stops before any alignment: paste its log.

## Block D (ED.6d). On the Mac: the copy back

Into each ctx root's `data/raw/<slug>_align/`, beside the `<slug>_meta/` and `<slug>_rnaseq/` folders
ED.4b's block 6 and the earlier copies use. `ctx qc --reads --align-dir <folder>` reads its
`tasks.tsv`, fastp JSON, FastQC zips, seqkit tables and `grid/` (ED.6e-1 to ED.6e-3), and
`ctx process --align-dir` its STAR and featureCounts tables (ED.6e-4, ED.6e-5); block E's featured
sample and references come with it. No command reads it unless `--align-dir` names it. The includes are the files M2's read stage and
M3 read (fastp JSON, FastQC zips, seqkit tables, the grid's JSON, STAR's two tables, flagstat,
idxstats, featureCounts matrices and summaries), plus `tasks.tsv`, `versions.txt` and the `time/`
files as the run's record, which the pages' Reproduce section needs (ED.6e-7).

```
cat > ~/ed6d_copy.sh <<'SCRIPT'
set -euo pipefail
K=rstak@kebnekaise.hpc2n.umu.se:/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/ed6_align
copy() {
    mkdir -p "$2/data/raw/$1_align"
    /usr/bin/rsync -avm --exclude='*.part*/' --include='*/' \
        --include='/tasks.tsv' --include='runs/*/fastp/*.fastp.json' --include='runs/*/fastqc/*_fastqc.zip' \
        --include='runs/*/seqkit/*.seqkit_stats.tsv' --include='runs/*/star/*.Log.final.out' \
        --include='runs/*/star/*.ReadsPerGene.out.tab' --include='runs/*/star/*.flagstat' --include='runs/*/star/*.idxstats' \
        --include='/featurecounts_annotation.tsv.gz' \
        --include='runs/*/featurecounts/*.featureCounts.tsv.gz' --include='runs/*/featurecounts/*.featureCounts.txt.summary' \
        --include='grid/*.fastp.json' --include='runs/*/time/*' \
        --include='featured/*/*/star/*.Log.final.out' --include='featured/*/*/star/*.ReadsPerGene.out.tab' \
        --include='featured/*/*/star/*.flagstat' --include='featured/*/*/star/*.idxstats' \
        --include='featured/*/cut/featurecounts_annotation.tsv.gz' \
        --include='featured/*/cut/featurecounts/*.featureCounts.tsv.gz' --include='featured/*/cut/featurecounts/*.featureCounts.txt.summary' \
        --include='featured/*/*/time/*' --include='featured/*/check/*' --exclude='*' "$K/$1/" "$2/data/raw/$1_align/"
    /usr/bin/rsync -av "$K/versions.txt" "$K/refs.txt" "$2/data/raw/$1_align/"
    # Block E's references: ctx's facts of the whole and the half, never the genome copies.
    /usr/bin/rsync -avm --include='/inputs.txt' --include='/main/' --include='/main/*.json' --include='/cut/' \
        --include='/cut/*.json' --exclude='*' "$K/reference/" "$2/data/raw/$1_align/reference/"
    python3 - "$2/data/raw/$1_align" <<'PY'
import csv, glob, os, sys
d = sys.argv[1]
tasks = list(csv.DictReader(open(f"{d}/tasks.tsv"), delimiter="\t"))
runs = sum(len(t["runs"].split(",")) for t in tasks)
files = sum(len(t["r1"].split(",")) * int(t["mates"]) for t in tasks)
count = lambda p: len(glob.glob(f"{d}/{p}"))
print(d, f"{len(tasks)} samples, {runs} runs, {files} FASTQ files expected")
for p, want in (("runs/*/fastp/*.fastp.json", runs), ("runs/*/fastqc/*_fastqc.zip", files),
                ("runs/*/seqkit/*.tsv", len(tasks)), ("runs/*/star/*.Log.final.out", len(tasks)),
                ("runs/*/star/*.ReadsPerGene.out.tab", len(tasks)), ("runs/*/star/*.flagstat", len(tasks)),
                ("runs/*/featurecounts/*.featureCounts.tsv.gz", 12 * len(tasks)), ("grid/*.fastp.json", 36),
                ("featurecounts_annotation.tsv.gz", 1), ("featured/*/*/star/*.Log.final.out", 2),
                ("featured/*/cut/featurecounts/*.featureCounts.tsv.gz", 3),
                ("reference/*/reference.json", 2), ("reference/inputs.txt", 1)):
    print(f"  {p}: {count(p)} of {want}")
print(f"  {sum(os.path.getsize(f) for f in glob.glob(f'{d}/**', recursive=True) if os.path.isfile(f))} bytes")
PY
}
copy tier1 ~/ctx_tier1
copy breast_tissue ~/ctx_breast
copy breast_mac ~/ctx_breast
copy breast_mono ~/ctx_breast
copy lung_gse229705 ~/ctx_lung_gse229705
copy pancreas_gse192647 ~/ctx_pancreas_gse192647
SCRIPT
bash ~/ed6d_copy.sh 2>&1 | tee ~/ed6d_copy.log
```

Paste back `~/ed6d_copy.log`: every "N of N" line equal. Then on Kebnekaise, and only then, the
pilots' BAMs and trimmed reads go:

```
cd /proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/ed6_align
for p in tier1 breast_tissue breast_mac breast_mono lung_gse229705 pancreas_gse192647; do s=$(cat "$p/pilot.txt"); ls -la "$p/runs/$s/star/$s.bam" "$p/runs/$s/trimmed"; done
# after reading that list:
for p in tier1 breast_tissue breast_mac breast_mono lung_gse229705 pancreas_gse192647; do s=$(cat "$p/pilot.txt"); rm -f "$p/runs/$s/star/$s.bam" "$p/runs/$s/star/$s.bam.bai"; rm -rf "$p/runs/$s/trimmed"; done
```

**The ctx commands on the Mac.** `ctx qc --reads --align-dir` reads the cluster's reports in
place of running FastQC, fastp and seqkit (ED.6e-1 to ED.6e-3). It refuses a folder whose
`tasks.tsv` and gather sheet disagree on the runs or their layout, a missing or unreadable report,
and a fastp report run on other files or at a filter other than the command's `--fastp-*`
settings. The grid is drawn on the sheet's first run unless `--grid-run` names another, and a
`grid/` run on other files is refused. `ctx process --align-dir` reads `<slug>_align/` (ED.6e-4,
ED.6e-5): no `--genome`, no FASTQ, the featured sample's `cut` and `untrimmed` from block E
(`--settings`, default `cut`). `ctx run --align-dir` passes the folder to process, which then runs
before qc (ED.6e-8).

```
cd /Users/rytis/Scripts/Pipeline/RNAseq_tutorial
# After ED.6e-1 to ED.6e-3: M2 with its read stage, from the cluster's reports.
pixi run -e laptop ctx --root ~/ctx_breast qc --dataset breast_tissue --dataset-file conf/datasets/breast_tissue.toml --rnaseq-dir ~/ctx_breast/data/raw/breast_tissue_rnaseq --reads --align-dir ~/ctx_breast/data/raw/breast_tissue_align
# After ED.6e-4 and ED.6e-5: M3 from the cluster's alignment, never the default.
pixi run -e laptop ctx --root ~/ctx_breast process --dataset breast_tissue --dataset-file conf/datasets/breast_tissue.toml --align-dir ~/ctx_breast/data/raw/breast_tissue_align --settings cut,untrimmed
# Then M2 again, for M3.12's align stage, which reads ctx process's tables.
pixi run -e laptop ctx --root ~/ctx_breast qc --dataset breast_tissue --dataset-file conf/datasets/breast_tissue.toml --rnaseq-dir ~/ctx_breast/data/raw/breast_tissue_rnaseq --reads --align-dir ~/ctx_breast/data/raw/breast_tissue_align
```

The same three lines for breast_mac and breast_mono in `~/ctx_breast`, tier1 in `~/ctx_tier1`, and
lung_gse229705 in `~/ctx_lung_gse229705`, pancreas_gse192647 in `~/ctx_pancreas_gse192647`, each
with its own dataset file and folders.

## What is not done here

- Only block A has run (Rytis, 2026-10-22). Nothing in the cost table is measured; block B2 measures it.
- The lung and pancreas pilots are the largest samples downloaded when block B1 ran, not the
  projects' largest.
- pancreas_gse192647 was not in block A; its folder and FASTQs are unchecked by this plan.
- `/usr/bin/time` on Kebnekaise's compute nodes is assumed, as M3.11 assumed it; the align script
  refuses to start without it.
- The GENCODE `MD5SUMS` line format (`<md5>  <file>`) is assumed from GENCODE's release folders and
  M3.11's sums; block B1 stops if either grep finds no line.
- `--multimap-max` is not run on the cluster: locally it changes every sample's alignment, so it
  would be block C1 again with `--outFilterMultimapNmax`; block E says why, and `ctx process
  --align-dir` refuses it. `cut` and `untrimmed` run in block E (ED.6e-6).
- `ctx process --align-dir` is tested on a fixture folder laid out as block D copies it
  (tests/process_cluster.py), single-end only, with STAR's two tables written in STAR's layout from
  HISAT2 BAMs, since the pinned osx-arm64 STAR maps nothing (D-036). No paired cluster folder and no
  real STAR `Log.final.out` has been read by it yet.
- nf-core's own FastQC zips on the raw reads may already sit in each rnaseq outdir (block A counts
  them); M2's read stage runs its own FastQC, so this plan runs it again rather than map nf-core's names.
