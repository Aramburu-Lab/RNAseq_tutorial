# M3.11 with M4's Q2: Kebnekaise blocks (written 2026-10-09, not run)

Approved by Rytis 2026-10-09, work-order item 13 answers 1 and 5. Seven blocks, one per wait. The
login node blocks are pasted whole, so every check sits inside one `if` on the file it reads, and a
job is submitted only when the files it needs pass their check. No block uses `set -e`; the batch
scripts written by the heredocs do.

## What is fixed before any block runs

**The reads.** M3.9's trim folder `~/ctx_tier1/outputs/03_process/tools/tier1_teach/trim_b72d1253088d/reads/`,
fastp 1.1.0 with `-q 15 -u 40 -l 15`, no cap, as its `ctx_manifest.json` records. Single-end, one run per
sample. The sample names are the ones M3's align folder `align_hisat2_a06356a7b142` keys, matched to the
runs by SHA-256. Bytes and reads were measured on the Mac on 2026-10-09 with `stat` and `gzip -dc | wc -l`
divided by 4.

| sample | run | file | bytes | reads | SHA-256 |
|---|---|---|---|---|---|
| ERX6994352 | ERR7425270 | ERR7425270.fastq.gz | 60,153,995 | 989,688 | dfa1b32700e5a5cacd9ff766926cf2756c3915defdeb5fdfc8dc6f84f02abc2c |
| ERX6994353 | ERR7425271 | ERR7425271.fastq.gz | 58,379,816 | 990,782 | bf4f48dea604560ecadb551251841b2a4d0b4ee0c47cb949d2f27e309a95d280 |
| ERX6994370 | ERR7425288 | ERR7425288.fastq.gz | 55,300,977 | 998,086 | 4611beb5cd569eb6ae24e198b7acb791cc7e21dc554ebaa5ca8b2816d8142c0f |
| ERX6994371 | ERR7425289 | ERR7425289.fastq.gz | 62,224,333 | 997,712 | 5a040cf1003035c7d97ebcb61d1249f3cb2dbd7d71001958c5d5eb8f86ce813a |
| ERX6994463 | ERR7425381 | ERR7425381.fastq.gz | 54,128,718 | 994,661 | 237b9d594d0e33791a511ec384ea7fdbf7388c240169d13a33340649e1a68ec3 |
| ERX6994464 | ERR7425382 | ERR7425382.fastq.gz | 53,617,743 | 994,728 | 04aa7717e70fec64311f33587fcf2c1f46c6ed6ac8a84644be796a0d8aba7aa5 |
| ERX6994487 | ERR7425405 | ERR7425405.fastq.gz | 60,674,339 | 998,321 | 6a0427b5f6c6528bb655b6c59fcc43a4c29dd6ed3c8497347f6e8438cca84b10 |
| ERX6994488 | ERR7425406 | ERR7425406.fastq.gz | 61,447,418 | 998,616 | 7894fb5339aad2305571689a1457ae724af685e2e8784deb0906da8b5c6febec |
| total | | 8 files | 465,927,339 | 7,962,594 | |

**STAR, as M3 runs it.** `src/ctx/process/reference.py` builds the index with `--runMode genomeGenerate
--sjdbGTFfile --sjdbOverhang 100 --genomeSAindexNbases <min(14, floor(log2(L)/2 - 1))>`, into a `.part`
folder renamed when done; for the 3.1 Gb genome that rule gives 14. `src/ctx/process/align.py`
`STAR_OPTIONS` and `star_align` map each sample with `--readFilesCommand gzip -dcf --quantMode GeneCounts
--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outFileNamePrefix <out>/<sample>.` and no
`--outFilterMultimapNmax`, so STAR's 10 holds (D-036). STAR is bioconda's linux-64 2.7.11b build
`h5ca1c30_8`, the build pixi.lock pins, as the depot.galaxyproject.org image of that name: 31,514,624 B by
HEAD on 2026-10-09. M3 feeds STAR uncompressed files because STAR refuses gzipped FASTA, so the job unzips
both. One departure, which changes no index byte: `--limitGenomeGenerateRAM 42000000000`, because STAR's
default of 31,000,000,000 B sits under the 32 GB its README asks for a mammal genome.

**Salmon, as nf-core ran it.** Step-1's run of nf-core/rnaseq 3.27.0 records, for ERX6994352
(`~/ctx_tier1/data/raw/tier1_salmon/ERX6994352`): salmon 2.7.0, "rust port, reads mode"; libType U;
bias correction off; optimiser vb; 644,914 targets and 194 decoys; `index_seq_hash`
`b05bb28ed408deeeabeb73b8714fbb87618eba5ae4895e1efe53bb1074ddf3d8`; and a gene map, since
`aux_info/genemap_unmatched_txps.json` and `quant.genes.sf` exist although its `cmd_info.json` lists no
`geneMap`. That `cmd_info.json` is therefore not a full record of the options, so block 0 prints nf-core's
own `.command.sh` for that task (trace hash `aa/693654`) if `work/` still holds it. The jobs run
`salmon quant --geneMap <nf-core's filtered GTF> --threads 4 --libType=U --index <nf-core's index> -r <reads>`,
plus `--dumpEq` in set A and `--gcBias --seqBias` in set B. **Stop and paste** if block 0's `.command.sh`
shows any other option. Whether 2.7.0 still takes `--dumpEq`, `--gcBias` and `--seqBias` is UNVERIFIED;
block 0 greps its help, and Salmon 1.x documented `--gcBias` as beta for single-end reads.

**Costs.** Measured values have a source; the rest are marked.

| Step | Request | Basis |
|---|---|---|
| Upload | 465,927,339 B | measured on the Mac |
| Unzipped genome and GTF, `$P/refs/unzipped/` | about 3.1 GB and 4,688,772,094 B | the GTF measured on the Mac; the genome from nf-core's `GUNZIP_FASTA` trace, 2.9 GB written (1024-based) |
| STAR index | 8 cores, 48G, 3 h | UNMEASURED: no GRCh38 STAR index has been built in this project. RAM: the STAR README's 32 GB. Time: a guess, about 8 times nf-core's measured `SALMON_INDEX` of the same genome (21 min 56 s). Disk: unmeasured; block 4 measures it |
| STAR alignment, 8 samples | 8 cores, 48G, 2 h | UNMEASURED; the index is read from disk once per sample. BAMs about 31 MB per 1M reads on chr22 (M3 plan), unmeasured on the full genome |
| Q2, 16 Salmon runs | 4 cores, 16G, 1 h | the request approved in the M4 plan; nf-core's subsampled strandedness runs took 15 to 37 s at 10.0 to 10.3 GB each |
| Q2 outputs | about 0.55 GB | `quant.sf` 29,878,586 B and `quant.genes.sf` 3,931,602 B per sample in nf-core's run, times 16; `--dumpEq`'s files unmeasured |

Nextflow does not run in these blocks, so `monochrome_logs` has nothing to switch off. Every job writes
`slurm-%x-%j.out` in `$X`, as P3.0's do. The time limits fall to the measured values once block 4 has
run; for M3.13 nothing here is a basis until then.

---

## Block 0. Kebnekaise login node: inputs, space and the two images

Produces: `$X/{reads,containers,logs}`, the STAR image (31,514,624 B), `$X/salmon_image.txt`.
Paste back: the whole output.

```
# Block 0. Kebnekaise login node. Paste back the whole output.
export P=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/tier1
export X=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/tier1_teach/full_reference
C=/proj/nobackup/hpc2n2024-150/Rytis/Cache/Apptainer/nf-core_rnaseq_3.27.0
G="$P/nfcore_rnaseq/results/genome"
mkdir -p "$X/reads" "$X/containers" "$X/logs"
df -h /proj/nobackup/hpc2n2024-150          # stop if under about 100G free: the index size is unmeasured
if test -s "$P/refs/GRCh38.primary_assembly.genome.fa.gz" && test -s "$P/refs/gencode.v50.primary_assembly.annotation.gtf.gz" && test -d "$G/index/salmon"; then
    echo "refs and nf-core's salmon index present"
    (cd "$P/refs" && printf '%s  %s\n' da1a11258be075cfa7af718162c894e7 GRCh38.primary_assembly.genome.fa.gz 289b91e5e95e8b0450d223246f10a12e gencode.v50.primary_assembly.annotation.gtf.gz | md5sum -c -)   # expect two OK
    ls -la "$P/refs"                                                  # expect no star_2.7.11b_gencode50_sjdb100 and no unzipped/ yet
    find "$G" -maxdepth 1 -printf '%y %s %f\n' | sort -k3             # nf-core's genome folder: expect one *.filtered.gtf
    du -sh "$G/index/salmon"; ls "$G/index/salmon"
    W=$(ls -d "$P"/nfcore_rnaseq/work/aa/693654* 2>/dev/null | head -n 1)
    if [ -n "$W" ]; then
        echo "nf-core's SALMON_QUANT (ERX6994352) command:"; grep -n 'salmon' "$W/.command.sh"
        IMG=$(grep -o '[^ ]*\.img' "$W/.command.run" | sort -u)
    else
        echo "no work dir aa/693654*: work/ was removed; the image comes from the cache"
        IMG=$(ls "$C"/*salmon*.img 2>/dev/null)
    fi
    ls "$C" | grep -i salmon
    if [ "$(printf '%s\n' "$IMG" | grep -c .)" = 1 ] && test -s "$IMG"; then
        echo "$IMG" | tee "$X/salmon_image.txt"
        apptainer exec "$IMG" salmon --version                        # expect 2.7.0
        apptainer exec "$IMG" salmon quant --help 2>&1 | grep -o -e '--dumpEq' -e '--gcBias' -e '--seqBias' -e '--geneMap' -e '--libType' | sort | uniq -c   # expect all five
    else
        rm -f "$X/salmon_image.txt"; echo "NOT ONE salmon image: [$IMG]"
    fi
    test -s "$X/containers/star_2.7.11b.sif" || wget -nv -O "$X/containers/star_2.7.11b.sif" https://depot.galaxyproject.org/singularity/star:2.7.11b--h5ca1c30_8
    stat -c '%s %n' "$X/containers/star_2.7.11b.sif"                   # expect 31514624
    sha256sum "$X/containers/star_2.7.11b.sif"
    apptainer exec "$X/containers/star_2.7.11b.sif" STAR --version     # expect 2.7.11b
else
    echo "NO refs or no salmon index under $P: nothing checked"
    ls -la "$P/refs" "$G" 2>&1 | head -n 30
fi
```

## Block 1. On the Mac: upload the 8 trimmed runs

Produces: 8 files, 465,927,339 B, in `$X/reads/` on the cluster, and the trim folder's manifest
beside them. Paste back: the whole output.

```
# Block 1. On the Mac, after block 0. Paste back the whole output.
T=~/ctx_tier1/outputs/03_process/tools/tier1_teach/trim_b72d1253088d
X=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/tier1_teach/full_reference
ssh rstak@kebnekaise.hpc2n.umu.se "test -d $X/reads" && \
rsync -av "$T"/reads/*.fastq.gz rstak@kebnekaise.hpc2n.umu.se:"$X"/reads/ && \
rsync -av "$T"/ctx_manifest.json rstak@kebnekaise.hpc2n.umu.se:"$X"/reads/trim_ctx_manifest.json || echo "NOT COPIED: no $X/reads on the cluster (run block 0), or the copy failed"
ls "$T"/reads/*.fastq.gz | wc -l                                      # expect 8
```

## Block 2. Kebnekaise: check the upload, submit the STAR index and Q2

Produces: `$X/reads_check.txt`, `$X/samples.tsv`, `$X/star_index.sbatch`, `$X/q2_salmon.sbatch`,
`$X/jobs.txt`, and two jobs. Paste back: the whole output.

```
# Block 2. Kebnekaise login node, after block 1. Paste back the whole output.
export P=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/tier1
export X=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/tier1_teach/full_reference
cd "$X"
cat > "$X/reads/SHA256SUMS" <<'EOF'
dfa1b32700e5a5cacd9ff766926cf2756c3915defdeb5fdfc8dc6f84f02abc2c  ERR7425270.fastq.gz
bf4f48dea604560ecadb551251841b2a4d0b4ee0c47cb949d2f27e309a95d280  ERR7425271.fastq.gz
4611beb5cd569eb6ae24e198b7acb791cc7e21dc554ebaa5ca8b2816d8142c0f  ERR7425288.fastq.gz
5a040cf1003035c7d97ebcb61d1249f3cb2dbd7d71001958c5d5eb8f86ce813a  ERR7425289.fastq.gz
237b9d594d0e33791a511ec384ea7fdbf7388c240169d13a33340649e1a68ec3  ERR7425381.fastq.gz
04aa7717e70fec64311f33587fcf2c1f46c6ed6ac8a84644be796a0d8aba7aa5  ERR7425382.fastq.gz
6a0427b5f6c6528bb655b6c59fcc43a4c29dd6ed3c8497347f6e8438cca84b10  ERR7425405.fastq.gz
7894fb5339aad2305571689a1457ae724af685e2e8784deb0906da8b5c6febec  ERR7425406.fastq.gz
EOF
printf 'sample\trun\treads\nERX6994352\tERR7425270\t989688\nERX6994353\tERR7425271\t990782\nERX6994370\tERR7425288\t998086\nERX6994371\tERR7425289\t997712\nERX6994463\tERR7425381\t994661\nERX6994464\tERR7425382\t994728\nERX6994487\tERR7425405\t998321\nERX6994488\tERR7425406\t998616\n' > "$X/samples.tsv"
# Counted both ways, so a missing file is a failure and not a silent 0.
(cd "$X/reads" && sha256sum -c SHA256SUMS) 2>&1 | awk '/: OK$/{ok++; next} {bad++} END{print ok+0, "OK;", bad+0, "not OK"}' | tee "$X/reads_check.txt"   # expect 8 OK; 0 not OK
find "$X/reads" -name '*.fastq.gz' -printf '%s\n' | awk '{n++; s+=$1} END{printf "%d files, %.0f bytes\n", n, s}'   # expect 8 files, 465927339 bytes

cat > "$X/star_index.sbatch" <<'EOF'
#!/usr/bin/env bash
#SBATCH -J m311_star_index
#SBATCH -A hpc2n2026-161
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 0-03:00:00
#SBATCH -o slurm-%x-%j.out
set -euo pipefail
cd "$X"
IMG="$X/containers/star_2.7.11b.sif"
IDX="$P/refs/star_2.7.11b_gencode50_sjdb100"
U="$P/refs/unzipped"
FA=GRCh38.primary_assembly.genome.fa
GTF=gencode.v50.primary_assembly.annotation.gtf
[ ! -e "$IDX" ] || { echo "$IDX exists: move it aside first"; exit 1; }
[ ! -e "$IDX.part" ] || { echo "$IDX.part exists, from a failed build: move it aside first"; exit 1; }
(cd "$P/refs" && printf '%s  %s\n' da1a11258be075cfa7af718162c894e7 "$FA.gz" 289b91e5e95e8b0450d223246f10a12e "$GTF.gz" | md5sum -c -)
apptainer exec "$IMG" STAR --version
mkdir -p "$U"
for f in "$FA" "$GTF"; do
    if [ ! -s "$U/$f" ]; then gzip -dc "$P/refs/$f.gz" > "$U/$f.part"; mv "$U/$f.part" "$U/$f"; fi
done
ls -l "$U"
mkdir -p "$IDX.part"
/usr/bin/time -v -o "$X/logs/star_index.time" apptainer exec -B /proj/nobackup/hpc2n2024-150 "$IMG" STAR \
    --runMode genomeGenerate --runThreadN 8 --genomeDir "$IDX.part" \
    --genomeFastaFiles "$U/$FA" --sjdbGTFfile "$U/$GTF" \
    --sjdbOverhang 100 --genomeSAindexNbases 14 \
    --limitGenomeGenerateRAM 42000000000 --outFileNamePrefix "$IDX.part/"
mv "$IDX.part" "$IDX"
du -sb "$IDX"
EOF

cat > "$X/q2_salmon.sbatch" <<'EOF'
#!/usr/bin/env bash
#SBATCH -J q2_salmon
#SBATCH -A hpc2n2026-161
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 0-01:00:00
#SBATCH -o slurm-%x-%j.out
set -euo pipefail
cd "$X"
IMG=$(cat "$X/salmon_image.txt")
SI="$P/nfcore_rnaseq/results/genome/index/salmon"
GTF=$(find "$P/nfcore_rnaseq/results/genome" -maxdepth 1 -name '*.filtered.gtf')
[ "$(printf '%s\n' "$GTF" | grep -c .)" = 1 ] || { echo "not one *.filtered.gtf: [$GTF]"; exit 1; }
for d in salmon_nfcore_dumpeq salmon_gcseqbias; do
    [ ! -e "$X/$d" ] && [ ! -e "$X/$d.part" ] || { echo "$X/$d or its .part exists: move it aside first"; exit 1; }
done
apptainer exec "$IMG" salmon --version
for pair in ERX6994352:ERR7425270 ERX6994353:ERR7425271 ERX6994370:ERR7425288 ERX6994371:ERR7425289 ERX6994463:ERR7425381 ERX6994464:ERR7425382 ERX6994487:ERR7425405 ERX6994488:ERR7425406; do
    s=${pair%%:*}; r=${pair#*:}
    /usr/bin/time -v -o "$X/logs/salmon_nfcore_dumpeq_$s.time" apptainer exec -B /proj/nobackup/hpc2n2024-150 "$IMG" salmon quant \
        --geneMap "$GTF" --threads 4 --libType=U --index "$SI" -r "$X/reads/$r.fastq.gz" \
        --dumpEq -o "$X/salmon_nfcore_dumpeq.part/$s" < /dev/null
    /usr/bin/time -v -o "$X/logs/salmon_gcseqbias_$s.time" apptainer exec -B /proj/nobackup/hpc2n2024-150 "$IMG" salmon quant \
        --geneMap "$GTF" --threads 4 --libType=U --index "$SI" -r "$X/reads/$r.fastq.gz" \
        --gcBias --seqBias -o "$X/salmon_gcseqbias.part/$s" < /dev/null
done
mv "$X/salmon_nfcore_dumpeq.part" "$X/salmon_nfcore_dumpeq"
mv "$X/salmon_gcseqbias.part" "$X/salmon_gcseqbias"
EOF

# Submitted only when the upload, both images and the STOP rule of block 0 hold.
if grep -qx '8 OK; 0 not OK' "$X/reads_check.txt" && test -s "$X/containers/star_2.7.11b.sif" && test -s "$X/salmon_image.txt"; then
    sbatch --parsable --export=ALL,P="$P",X="$X" "$X/star_index.sbatch" | sed 's/^/star_index /' | tee -a "$X/jobs.txt"
    sbatch --parsable --export=ALL,P="$P",X="$X" "$X/q2_salmon.sbatch" | sed 's/^/q2_salmon /' | tee -a "$X/jobs.txt"
else
    echo "NOT SUBMITTED: the upload is not 8 OK, or an image is missing (block 0)"
fi
squeue -u "$USER"
```

## Block 3. Kebnekaise, once `m311_star_index` has ended: check the index, submit the alignment

Produces: `$X/star_align.sbatch` and one job. Paste back: the whole output.

```
# Block 3. Kebnekaise login node, once squeue shows no m311_star_index job. Paste back the whole output.
export P=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/tier1
export X=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/tier1_teach/full_reference
IDX="$P/refs/star_2.7.11b_gencode50_sjdb100"
cd "$X"
if test -s "$IDX/SA" && test -s "$IDX/genomeParameters.txt"; then
    echo "index present"
    tail -n 15 "$X"/slurm-m311_star_index-*.out
    grep -e 'Elapsed (wall' -e 'Maximum resident' "$X/logs/star_index.time"
    sacct -n -j "$(awk '$1 == "star_index" {print $2}' "$X/jobs.txt")" --format=JobID,JobName%20,Elapsed,MaxRSS,ReqMem,AllocCPUS,State   # without -X, so the batch step's memory shows
    du -sh "$IDX" "$P/refs/unzipped"
    cat "$IDX/genomeParameters.txt"
    wc -l "$IDX/chrName.txt" "$IDX/sjdbList.out.tab" "$IDX/exonGeTrInfo.tab" "$IDX/geneInfo.tab"   # geneInfo.tab: expect 78942, a count line and 78,941 genes
    cat > "$X/star_align.sbatch" <<'EOF'
#!/usr/bin/env bash
#SBATCH -J m311_star_align
#SBATCH -A hpc2n2026-161
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 0-02:00:00
#SBATCH -o slurm-%x-%j.out
set -euo pipefail
cd "$X"
IMG="$X/containers/star_2.7.11b.sif"
IDX="$P/refs/star_2.7.11b_gencode50_sjdb100"
[ ! -e "$X/star" ] && [ ! -e "$X/star.part" ] || { echo "$X/star or star.part exists: move it aside first"; exit 1; }
mkdir -p "$X/star.part"
# M3's star_align: src/ctx/process/align.py STAR_OPTIONS, no --outFilterMultimapNmax (D-036).
for pair in ERX6994352:ERR7425270 ERX6994353:ERR7425271 ERX6994370:ERR7425288 ERX6994371:ERR7425289 ERX6994463:ERR7425381 ERX6994464:ERR7425382 ERX6994487:ERR7425405 ERX6994488:ERR7425406; do
    s=${pair%%:*}; r=${pair#*:}
    /usr/bin/time -v -o "$X/logs/star_$s.time" apptainer exec -B /proj/nobackup/hpc2n2024-150 "$IMG" STAR \
        --runThreadN 8 --genomeDir "$IDX" --readFilesIn "$X/reads/$r.fastq.gz" \
        --readFilesCommand gzip -dcf --quantMode GeneCounts \
        --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
        --outFileNamePrefix "$X/star.part/$s." < /dev/null
done
mv "$X/star.part" "$X/star"
EOF
    sbatch --parsable --export=ALL,P="$P",X="$X" "$X/star_align.sbatch" | sed 's/^/star_align /' | tee -a "$X/jobs.txt"
else
    echo "NO INDEX at $IDX: the build has not finished or has failed; nothing submitted"
    tail -n 30 "$X"/slurm-m311_star_index-*.out
    squeue -u "$USER"
fi
```

## Block 4. Kebnekaise, once `m311_star_align` and `q2_salmon` have ended: the checks

Produces nothing new. Paste back: the whole output.

```
# Block 4. Kebnekaise login node, once squeue shows no m311_star_align and no q2_salmon job. Paste back the whole output.
export P=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/tier1
export X=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/tier1_teach/full_reference
cd "$X"
if test -d "$X/star" && test -d "$X/salmon_nfcore_dumpeq" && test -d "$X/salmon_gcseqbias"; then
    echo "star and both salmon sets present"
    tail -n 5 "$X"/slurm-m311_star_align-*.out "$X"/slurm-q2_salmon-*.out
    sacct -n -j "$(awk '{print $2}' "$X/jobs.txt" | paste -sd,)" --format=JobID,JobName%20,Elapsed,MaxRSS,ReqMem,AllocCPUS,State
    grep -H -e 'Elapsed (wall' -e 'Maximum resident' "$X"/logs/*.time
    # STAR: every read read, and each ReadsPerGene column summing to the input, as M3's probe asks.
    printf 'sample\texpected\tinput\tunique_pct\tmulti_pct\tcol2\tcol3\tcol4\n'
    tail -n +2 "$X/samples.tsv" | while IFS=$'\t' read -r s r n; do
        L="$X/star/$s.Log.final.out"; RPG="$X/star/$s.ReadsPerGene.out.tab"
        inp=$(awk -F'|' '/Number of input reads/ {gsub(/[ \t]/, "", $2); print $2}' "$L")
        uni=$(awk -F'|' '/Uniquely mapped reads %/ {gsub(/[ \t]/, "", $2); print $2}' "$L")
        mul=$(awk -F'|' '/% of reads mapped to multiple loci/ {gsub(/[ \t]/, "", $2); print $2}' "$L")
        sums=$(awk -F'\t' '{a += $2; b += $3; c += $4} END {print a "\t" b "\t" c}' "$RPG")
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$s" "$n" "$inp" "$uni" "$mul" "$sums"
    done                                                              # expect input = expected = col2 = col3 = col4 on every row
    wc -l "$X"/star/*.ReadsPerGene.out.tab | sort -u -k1,1 | head -n 3   # expect 78945 lines each: the GTF's 78,941 gene lines (counted on the Mac) and STAR's 4 N_ lines
    ls "$X"/star | sed 's/^ERX[0-9]*\.//' | sort | uniq -c            # expect 8 of each STAR file
    # Salmon: nf-core's index (its seq hash), every read processed, bias as asked, --dumpEq's files.
    python3 - "$X" <<'PY'
import json, os, sys
x = sys.argv[1]
want = {l.split("\t")[0]: int(l.split("\t")[2]) for l in open(os.path.join(x, "samples.tsv")).read().splitlines()[1:]}
print("set\tsample\tprocessed\texpected\tmapped_pct\tdecoy\tlibtype\tseq_bias\tgc_bias\tsame_index\taux_info")
for d in ("salmon_nfcore_dumpeq", "salmon_gcseqbias"):
    for s, n in want.items():
        m = json.load(open(os.path.join(x, d, s, "aux_info", "meta_info.json")))
        aux = ",".join(sorted(os.listdir(os.path.join(x, d, s, "aux_info"))))
        print("\t".join(map(str, [d, s, m["num_processed"], n, round(m["percent_mapped"], 2), m["num_decoy_fragments"], ",".join(m["library_types"]),
              m["seq_bias_correct"], m["gc_bias_correct"], m["index_seq_hash"] == "b05bb28ed408deeeabeb73b8714fbb87618eba5ae4895e1efe53bb1074ddf3d8", aux])))
PY
    # expect processed = expected, U, False False in the first set and True True in the second, same_index True,
    # and an eq_classes file in the first set's aux_info only
    grep -il -e warn -e error "$X"/salmon_*/ERX*/logs/salmon_quant.log | sed 's/^/log with a warning or error: /'
    du -sh "$X/star" "$X/salmon_nfcore_dumpeq" "$X/salmon_gcseqbias" "$X"/star/*.bam | tail -n 12
else
    echo "NOT ALL OUTPUTS: a job has not finished or has failed; nothing checked"
    ls -d "$X"/star* "$X"/salmon_* 2>&1
    tail -n 30 "$X"/slurm-m311_star_align-*.out "$X"/slurm-q2_salmon-*.out 2>&1
    squeue -u "$USER"
fi
```

## Block 5. On the Mac: copy back

Produces `~/ctx_tier1/data/raw/tier1_teach_full_reference/` (a new folder; its name is a proposal):
`star/` with 8 `Log.final.out`, 8 `ReadsPerGene.out.tab` and 8 `Log.out`; `star_index/` with the
index's `genomeParameters.txt` and `Log.out`; both Salmon sets whole, about 0.55 GB plus `--dumpEq`'s
files; `samples.tsv`, `jobs.txt`, `logs/*.time` and the jobs' `slurm-*.out`. The BAMs, the index and
the reads stay on the cluster. Paste back: the whole output.

```
# Block 5. On the Mac, after block 4 printed "star and both salmon sets present". Paste back the whole output.
X=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/tier1_teach/full_reference
IDX=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/tier1/refs/star_2.7.11b_gencode50_sjdb100
D=~/ctx_tier1/data/raw/tier1_teach_full_reference
ssh rstak@kebnekaise.hpc2n.umu.se "test -s $X/star/ERX6994488.ReadsPerGene.out.tab && test -s $X/salmon_nfcore_dumpeq/ERX6994488/quant.sf && test -s $X/salmon_gcseqbias/ERX6994488/quant.sf" && mkdir -p "$D/star_index" && \
rsync -avm --include='*/' --include='star/*.Log.final.out' --include='star/*.ReadsPerGene.out.tab' --include='star/*.Log.out' --include='salmon_nfcore_dumpeq/**' --include='salmon_gcseqbias/**' --include='samples.tsv' --include='jobs.txt' --include='logs/*.time' --include='slurm-*.out' --exclude='*' rstak@kebnekaise.hpc2n.umu.se:"$X"/ "$D"/ && \
rsync -av rstak@kebnekaise.hpc2n.umu.se:"$IDX"/genomeParameters.txt rstak@kebnekaise.hpc2n.umu.se:"$IDX"/Log.out "$D"/star_index/ || echo "NOT COPIED: the outputs are not all on the cluster, or the copy failed"
if test -d "$D/star"; then
    for p in 'star/*.Log.final.out' 'star/*.ReadsPerGene.out.tab' 'salmon_nfcore_dumpeq/*/quant.sf' 'salmon_gcseqbias/*/quant.sf' 'salmon_nfcore_dumpeq/*/aux_info/eq_classes*'; do
        printf '%s: %s\n' "$p" "$(cd "$D" && ls $p 2>/dev/null | wc -l | tr -d ' ')"
    done                                                              # expect 8 on every line
    echo "$(find "$D" -type f | wc -l | tr -d ' ') files, $(find "$D" -type l | wc -l | tr -d ' ') symlinks"   # expect 0 symlinks
    du -sh "$D"
else
    echo "NO star folder in $D: nothing checked"
fi
```

## What is not done here

- No block was run: this session cannot reach Kebnekaise.
- The STAR index's RAM, time and disk, and the alignment's, are unmeasured; their limits are flagged
  estimates, and block 3 and block 4 print the measured values that replace them for M3.13.
- featureCounts on the full-genome BAMs (the M4 plan's "M3.11 adds featureCounts") is not in these
  blocks: the task named STAR's two files and the Salmon outputs. The BAMs stay in `$X/star/` for it.
- `ctx process` reading these files as the full-reference side is M3.11's code step, not a cluster step.
