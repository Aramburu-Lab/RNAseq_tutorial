# NF.2: nf-core/rnaseq 3.27.0 with STAR on the six core projects (D-113)

Run by Rytis on Kebnekaise 2026-10-23 at the clone's 0e54360. Each project's Salmon-only
`results/` was moved to `results_salmon_only_20261023`; its references came from its own
`params.yml` (lung, pancreas) or nf-core's `results/pipeline_info/params_*.json` (tier1 from
2026-09-30_09-26-16, breast_tissue 2026-10-01_09-48-27, breast_mac 2026-09-30_14-05-35,
breast_mono 2026-10-02_08-57-08). A first submission (48070934, 48070935) stopped in 2 s:
`kebnekaise_env.sh` needs `CTX_STORE`, which that script did not export; this one exports it and
activates the environment once before submitting.

| project | head job |
|---|---|
| tier1 | 48071030 |
| breast_tissue | 48071031 |
| breast_mac | 48071032 |
| breast_mono | 48071033 |
| lung_gse229705 | 48071034 |
| pancreas_gse192647 | 48071035 |

```
cd /proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial
cat > nf_launch.sh <<'LAUNCH'
set -euo pipefail
S=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial
A=hpc2n2026-161
# Both are read by kebnekaise_env.sh inside every job; without CTX_STORE the 10:16 jobs stopped in 2 s.
export CTX_STORE="$S"
export CTX_REPO="$CTX_STORE/RNAseq_tutorial"
git -C "$CTX_REPO" pull --ff-only -q
git -C "$CTX_REPO" log -1 --format='%h %s'
[[ -f "$CTX_REPO/scripts/bash/job_name.sh" ]] || { echo "the clone predates the job names"; exit 1; }
grep -q '^SITE_RNASEQ_ALIGN=true$' "$CTX_REPO/scripts/bash/site_kebnekaise.sh" || { echo "the clone predates D-113"; exit 1; }
# The activation the dry run skips, tried once here, so a job never dies on it in 2 s.
( source "$CTX_REPO/scripts/bash/kebnekaise_env.sh" && nextflow -version | grep -m1 version ) || { echo "the fetch environment does not activate"; exit 1; }
stamp=$(date +%Y%m%d)
for p in tier1 breast_tissue breast_mac breast_mono lung_gse229705 pancreas_gse192647; do
    RUN="$S/$p"
    L="$RUN/nfcore_rnaseq"
    [[ -f "$RUN/nfcore_fetchngs/results.complete" ]] || { echo "$p: no results.complete; skipped"; continue; }
    if squeue -h -u "$USER" -o "%j %Z" | awk -v r="$RUN" '$1 ~ /_nf-head_|ctx_rnaseq/ && $2 == r' | grep -q .; then
        echo "$p: a head job already runs in $RUN; skipped"; continue
    fi
    # The references of the Salmon-only run: our launcher's params.yml, else nf-core's own record
    # of the run (results/pipeline_info/params_*.json), for runs started by the lab's launcher.
    if [[ -s "$L/params.yml" ]]; then
        fasta=$(sed -n 's/^fasta: //p' "$L/params.yml"); gtf=$(sed -n 's/^gtf: //p' "$L/params.yml")
        gencode=$(sed -n 's/^gencode: //p' "$L/params.yml")
    else
        last=$(ls "$L"/results/pipeline_info/params_*.json 2>/dev/null | sort | tail -n 1)
        [[ -n "$last" ]] || { echo "$p: neither $L/params.yml nor a params_*.json under $L/results; skipped"; continue; }
        read -r fasta gtf gencode < <(python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); print(d["fasta"], d["gtf"], str(d.get("gencode", False)).lower())' "$last")
        echo "$p: references from $last"
    fi
    [[ -s "$fasta" && -s "$gtf" && -n "$gencode" ]] || { echo "$p: references not found: $fasta $gtf"; continue; }
    # The Salmon-only results move aside, never removed; ctx gather refuses an outdir whose params disagree.
    if [[ -d "$L/results" ]] && ! grep -qs '^aligner: star_salmon$' "$L/params.yml"; then
        mv "$L/results" "$L/results_salmon_only_$stamp"
        echo "$p: Salmon-only results moved to $L/results_salmon_only_$stamp"
    fi
    (cd "$RUN" && DRY_RUN=1 bash "$CTX_REPO/scripts/bash/rnaseq.sbatch" "$RUN" "$fasta" "$gtf" "$gencode" | tail -n 1)
    grep -q '^aligner: star_salmon$' "$L/params.yml" || { echo "$p: the dry run wrote no star_salmon route"; exit 1; }
    grep -q '_nf-' "$L/site.config" || { echo "$p: the dry run wrote no task names"; exit 1; }
    j=$(cd "$RUN" && sbatch --parsable -A "$A" "$CTX_REPO/scripts/bash/rnaseq.sbatch" "$RUN" "$fasta" "$gtf" "$gencode")
    echo "$p: rnaseq $j" | tee -a "$S/nf_rnaseq_jobs.txt"
done
squeue -u "$USER" -o "%.10i %.60j %.2t %.10M %R"
df -h /proj/nobackup/hpc2n2024-150
LAUNCH
bash nf_launch.sh 2>&1 | tee nf_launch.log
```
