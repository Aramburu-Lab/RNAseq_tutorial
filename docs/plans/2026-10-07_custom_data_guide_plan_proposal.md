# Your own data, from an accession to your own course pages: plan, 2026-10-07

Saved as [PLAN.md](/private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/plan_guide/PLAN.md). Paths start at `/Users/rytis/Scripts/Pipeline/RNAseq_tutorial`, lines as at b9254be; I changed nothing there, unlike another writer's uncommitted build_kit.py, naming.py and test_kit.py.

**Outcome.** One generated page, "Your own data", holds Rytis's three parts for a laptop and a cluster; each command comes from a Usage block, the CLI parser or the code's constants, and a test runs each in page order. Missing today: own-table curation, `ctx run --check`, site files for one launcher per pipeline, a `Code` block.

## Facts measured or read today

| Fact | Value | Source |
|---|---|---|
| Salmon index, GENCODE 50, 4 cluster runs | peak RSS 9.1, 9.1, 9.3, 16.8 GiB; 21.8 to 25.9 min; the genome indexed as decoys | `SALMON_INDEX` rows, `execution_trace_*.txt` of `~/ctx_tier1/data/raw/tier1_rnaseq`, `~/ctx_breast/data/raw/breast_*_rnaseq` (GB 1024-based, MemoryUnit.groovy:137-146); [salmon/index/main.nf:24-28](https://github.com/nf-core/rnaseq/blob/3.27.0/modules/nf-core/salmon/index/main.nf) |
| `SALMON_QUANT` | 10.2 to 10.3 GiB, all 858 tasks | same traces |
| GENCODE 50 genome, GTF | 845,635,028 and 124,650,284 B | HTTP HEAD, ftp.ebi.ac.uk release_50; tier1.toml:39 |
| G.0's four smallest Step-1 runs | 1,268,992,130 B, all with FTP links | ENA filereport; TODO.md:577 |
| Their Step-1 `num_processed` | 5,544,123; 5,723,703; 6,377,572; 6,333,288, each under ENA's `read_count` | `tier1_rnaseq/salmon/<ERX>/aux_info/meta_info.json` |
| `md5sum -c -` on empty input | exit 0 with macOS 15.7.9's /sbin/md5sum; exit 1 with GNU 9.12 from the `fetch` environment | run in scratch |

## 1. The guide

**Where.** One page, `outputs/09_report/<date>_report_guide_v1.html`, written by a new `ctx guide [--dataset-file F]`, linked from START.html as the glossary is (course.py:98-101), so in every kit. README.md gains a pointer: the kit README sends readers to commands (course.toml:61) the template-stub README lacks (README.md:7-8). No `docs/new_dataset.md` (TODO.md:634), a second copy. Words in `content/guide.toml`, prose-linted (CONVENTIONS.md:66-67).

**A dataset file drives it.** `[dataset] slug` and a new `accessions` list fill Parts 1 and 2. `ctx guide --dataset-file ~/ctx_mine/mine.toml` writes a student's own commands; the course's page uses `conf/datasets/four_runs.toml`.

**Part 1, download with fetchngs 1.12.0,** after `pixi install -e fetch --frozen` (pixi.toml:182-187):

```
mkdir -p "$RUN" && cd "$RUN"
printf '%s\n' ERR7425276 ERR7425397 ERR7425403 ERR7425436 > ids.txt
CTX_SITE=laptop SKIP_DOWNLOAD=1 bash "$CTX_REPO/scripts/bash/fetchngs.sbatch" "$RUN" ids.txt
CTX_SITE=laptop bash "$CTX_REPO/scripts/bash/fetchngs.sbatch" "$RUN" ids.txt
```

The cluster shows fetchngs.sbatch:19-30 after kebnekaise_env.sh:12-19; the laptop lines join that Usage block.

**Part 2, count with rnaseq 3.27.0, Salmon only:**

```
mkdir -p "$RUN/refs" && cd "$RUN/refs"
B=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50
curl -fL --remote-name-all "$B/MD5SUMS" "$B/GRCh38.primary_assembly.genome.fa.gz" "$B/gencode.v50.primary_assembly.annotation.gtf.gz"
cd "$RUN" && CTX_SITE=laptop bash "$CTX_REPO/scripts/bash/rnaseq.sbatch" "$RUN" refs/GRCh38.primary_assembly.genome.fa.gz refs/gencode.v50.primary_assembly.annotation.gtf.gz true
```

The cluster line starts `sbatch -A <compute project>` (rnaseq.sbatch:23-24). U.4 moves the MD5SUMS check into rnaseq.sbatch, run with the environment's GNU md5sum: on a Mac a shell check passes an empty match. GENCODE 50 is Step-1's release (TODO.md:2191).

**Part 3, into the course.** The cluster route copies back with stock macOS rsync, openrsync; on `tests/fixtures/rnaseq` it copied the 22 files rsync 3.4.3 copied under TODO.md:2326's includes:

```
/usr/bin/rsync -av <user>@kebnekaise.hpc2n.umu.se:$RUN/nfcore_fetchngs_meta/results/ ~/ctx_mine/data/raw/mine_meta/
/usr/bin/rsync -avmL --include='*/' <an --include per COPIED_FILES pattern> --exclude='*' <user>@kebnekaise.hpc2n.umu.se:$RUN/nfcore_rnaseq/results/ ~/ctx_mine/data/raw/mine_rnaseq/
```

Both routes:

```
pixi install -e laptop --frozen && pixi install -e reads --frozen
pixi run -e laptop ctx --root ~/ctx_mine gather --dataset mine --fetchngs-dir <meta outdir>
cp conf/dataset_template.toml ~/ctx_mine/mine.toml
pixi run -e laptop ctx --root ~/ctx_mine setup --dataset-file ~/ctx_mine/mine.toml
pixi run -e laptop ctx --root ~/ctx_mine run --dataset mine --dataset-file ~/ctx_mine/mine.toml --fetchngs-dir <meta outdir> --rnaseq-dir <rnaseq outdir>   # first with --check
pixi run -e laptop python scripts/python/build_kit.py --out ~/ctx_mine_kit --real mine=~/ctx_mine --zip
```

The uncurated gather writes the sheet a student's table starts from (harmonise.py:491); the kit builder takes any root with `outputs/` (build_kit.py:182-191).

**Another SLURM site** copies the Kebnekaise site file, changing the profile (`apptainer`, or its nf-core/configs profile, loaded by both pipelines: fetchngs nextflow.config:37,65, rnaseq nextflow.config:154,331; `uppmax` and `pdc_kth` exist, `hpc2n` not, api.github.com/repos/nf-core/configs/contents/conf), queueSize (rnaseq.sbatch:90), home rule, scratch variable (kebnekaise_env.sh:94) and `#SBATCH` header. Whether its compute nodes reach ENA, where downloads run (kebnekaise_env.sh:89), is unverified.

**Tests.** T1: each `ctx` line parses with the real click parser. T2: each other line is a Usage block's line; the blocks run in page order from a fresh folder under G.1's stubs (test_cluster_scripts.py:124-157), the recorded `nextflow` argv carrying pin, profile and `NXF_SYNTAX_PARSER`. T3: rsync includes come from COPIED_FILES (rnaseq.py:29-35), checked as test_repo.py:262-281 does. T4: every number fills a placeholder (CONVENTIONS.md:66). T5: `ctx run` writes the module commands' tables.

## 2. Tooling, smallest first

- **U.1 own table (S).** `[sdrf] file = "<table>"` in place of `resource`, relative to the dataset file, CSV or TSV, BOM and CRLF stripped, SHA-256 recorded; gather accepts it (harmonise.py:510-511). Mutations: path read from the working folder; tab forced; BOM kept; SHA-256 dropped; a run missing from the table accepted.
- **U.2 template (S, G.3).** `conf/dataset_template.toml`, three routes; MSigDB and GTF blocks. Mutations: resource blocks unlike tier1.toml's; the filled template failing to curate the fetchngs fixture.
- **U.3 trace (S, G.4).** `ctx.gather.trace`: peak RSS and real time per process, sample from the parentheses in `name`; a cost table on the gather page. Mutations: GB read as MB; CACHED rows (40 of the fixture's 41) counted as cost; sample from the tag.
- **U.4 sites (M).** Launchers read `CTX_SITE` (default kebnekaise) before the dry run; `laptop` on arm64: `docker,arm` for fetchngs, `docker,arm64` for rnaseq (fetchngs nextflow.config:110-112; rnaseq nextflow.config:224-228), fetchngs capped `--max_cpus 4 --max_memory 8.GB` (decisions.md:1622); the MD5SUMS check; `tutorial_test_run.py:80-99`, a second fetchngs command, retired. Mutations: Kebnekaise env on a laptop; profiles swapped; home refused on a laptop; an empty MD5SUMS match accepted; G.1's 92 planted mutations, rerun.
- **U.5 Code (S, G.2).** File text, SHA-256, Copy (builder.py:261; module.html.j2:55). Mutations: text reflowed; another file's hash; two buttons copying one block.
- **U.6 run (M, 2D).** `--check` lists every missing table with the lines to add, and unfetched resources; then gather, qc, de, function, report as subprocesses, each page keeping its Reproduce and peak RAM; after gather, each contrast's realised n via `select_samples` (de/deseq.py:120); stops at the first failure, naming it. Mutations: modules out of order, the test asserting order, not only the set; the check stopping at its first problem; a suggested line still failing once pasted; n counted before filters; tables unlike the module-by-module run.
- **U.7 kit (S).** A kit's dataset folder as a `--real` root. Mutations: that folder refused.
- **U.9 page (M).** `ctx guide`, `content/guide.toml`, START link. Tests: T1 to T4.

## 3. The proof

- **P1, every commit:** T1 to T5, the U tests.
- **P2, this Mac, laptop route**, once Rytis approves about 2.2 GB of reads and references plus images: a fresh clone in scratch, the page only, `/usr/bin/time -l` per step, current Docker memory first; on exit 137 Rytis raises it and `-resume` continues. Pass: the four `num_processed` above, 78,470 genes as Step-1, `ctx run` exit 0, kit checks pass; the measured floor goes on the page; Apptainer laptops and WSL2 stay "not run".
- **P3, Kebnekaise (Rytis, G.5)**, after U.4; one block per wait (lessons.md:36). Now, for the laptop floor:

```
du -sh /proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/tier1/nfcore_rnaseq/results/genome/index/salmon   # PASTE
du -sh /proj/nobackup/hpc2n2024-150/Rytis/Cache/Apptainer/nf-core_rnaseq_3.27.0                                   # PASTE
```

On the Mac: `git archive --format=tar.gz -o ~/ctx_guide_src.tar.gz HEAD`, `git rev-parse --short HEAD` (PASTE), scp to the Tasks folder. Then a new login:

```
export CTX_STORE=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/guide_rehearsal
export CTX_REPO="$CTX_STORE/RNAseq_tutorial" RUN="$CTX_STORE/four_runs"
mkdir -p "$CTX_REPO" "$RUN" && tar -xzf /proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial/ctx_guide_src.tar.gz -C "$CTX_REPO"
PIXI_CACHE_DIR="$CTX_STORE/pixi_cache" "${PIXI_HOME:-$HOME/.pixi}/bin/pixi" install --manifest-path "$CTX_REPO/pixi.toml" -e fetch --frozen
source "$CTX_REPO/scripts/bash/kebnekaise_env.sh" && nextflow -version 2>&1 | grep -i version     # PASTE, expect 26.04.6
cd "$RUN" && printf '%s\n' ERR7425276 ERR7425397 ERR7425403 ERR7425436 > ids.txt
SKIP_DOWNLOAD=1 sbatch -A hpc2n2026-161 "$CTX_REPO/scripts/bash/fetchngs.sbatch" "$RUN" ids.txt
```

When it ends:

```
cd "$RUN" && grep -h 'bytes of FASTQ at ENA' slurm-ctx_fetchngs-*.out          # PASTE, expect 4 accessions; 4 runs; 1,268,992,130 bytes
grep -q ' 1,268,992,130 bytes' slurm-ctx_fetchngs-*.out && sbatch -A hpc2n2026-161 "$CTX_REPO/scripts/bash/fetchngs.sbatch" "$RUN" ids.txt || echo "NOT SUBMITTED"
```

When the download ends, relaunch; every download should come back CACHED; a mid-run cancel is not rehearsed:

```
cd "$RUN" && test -f nfcore_fetchngs/results.complete && mv nfcore_fetchngs/results nfcore_fetchngs/results_try1 && sbatch -A hpc2n2026-161 "$CTX_REPO/scripts/bash/fetchngs.sbatch" "$RUN" ids.txt || echo "NOT SUBMITTED"
```

Then Part 2 with `sbatch -A hpc2n2026-161`; after it, the checks in one `if` (lessons.md:49):

```
S="$RUN/nfcore_rnaseq/results/salmon"
if test -s "$S/salmon.merged.gene_counts_length_scaled.tsv"; then
    for s in ERX6994358 ERX6994479 ERX6994485 ERX6994518; do printf '%s ' "$s"; grep -o '"num_processed": [0-9]*' "$S/$s/aux_info/meta_info.json"; done
    grep -hE 'nf-core/rnaseq|salmon:|Nextflow' "$RUN"/nfcore_rnaseq/results/pipeline_info/*versions*.yml
    sacct -n -S <first job's day> --name=ctx_fetchngs,ctx_rnaseq --format=JobID,JobName,Elapsed,MaxRSS,State
    du -sh "$RUN/nfcore_fetchngs" "$RUN/nfcore_rnaseq/results" "$RUN/nfcore_rnaseq/work"
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="status") s=i; next} {n[$s]++} END{for(k in n) print n[k], k}' "$RUN"/nfcore_fetchngs/results/pipeline_info/execution_trace_*.txt
else
    echo "NO MATRIX: nothing checked"; squeue -u "$USER"
fi                                                                 # PASTE all; expect 5544123, 5723703, 6377572, 6333288; every download CACHED
```

Then Part 3 on the Mac in `~/ctx_guide`, never `~/ctx_tier1` (TODO.md:582).
- **P4, newcomer (G.8):** someone new to the repository, given only the kit, runs the worked example, then a 2 to 4 run study of their own, logging every edit. Pass: none on the worked example; their pages reach a kit.

## 4. Order and parallel work

U.3, U.4, U.5 touch no module planner's files: start now. U.1, U.2 and U.6 touch `dataset.py`, `cli.py` and `cards.toml`, which M7 to M10 extend: small merges; `ctx run`'s order is one tuple each new module joins, report last (`MODULES`, naming.py:12-26, is not run order). U.9 needs U.4 to U.6 and U.8; P2 needs U.9; P3 needs U.4; P4 needs P2 and P3.

Absorbed: G.2 to G.8; Phase 2D's `ctx run` and new-dataset items (TODO.md:633-634), the doc replaced by the page; K.5's rerun pointer. Superseded: `tutorial_test_run.py`; TODO.md:158's typed command (test_repo.py:210 reads it); G.6's D-046 and D-047, numbers now the Step-1 record's and the kit's. Left: cohort bundles, containers, `ctx figure` (TODO.md:635-639), docs/accessions.md (TODO.md:659).

## 5. Decisions for Rytis; 1 to 4 block the build

1. **Where the guide lives.** Recommended: the generated page, READMEs pointing there; cost M. `docs/new_dataset.md`: S, a second copy, regex-checked at best (test_repo.py:199-259). Kit README: no Copy buttons, untested.
2. **Laptop counting.** Recommended: above a floor measured in P2, the cluster below it; cost: Rytis raises this Mac's Docker memory, 11.7 of 24 GiB by `docker info` and `sysctl`, under the 16.8 GiB index peak seen once. Laptops only fetch: no counts without a cluster account. Transcript-only index (prepare_genome_indices/main.nf:273-277): fits small laptops but changes the method; counts not comparable with Step-1's.
3. **One launcher per pipeline, site files.** Recommended; cost: G.1's tests rerun. Separate laptop scripts: second copies of the list checks, D-035's retry config and end checks.
4. **Worked example.** Recommended: G.0's four smallest runs throughout, curated from an own table as 2 neutrophils (LDN, NDN, both PDAC) against 2 monocytes (G-CSF, IPMN), confounding stated on the page. Alternative: UPN13 and UPN16, NDN and monocyte each, steady state, 101 nt, paired: 3,786,695,065 B, 3.0 times the bytes.
5. **G.0's two D-entries**, the student route and own-table curation: the next free numbers when U.1 lands, D-053 and D-054 unless taken first.
6. **G.8's newcomer:** recommended now, one colleague with read access, before release (TODO.md:591).
7. **Head memory (P4.5, TODO.md:627):** recommended 12G in rnaseq.sbatch before P3, from peaks of 6.16, 6.61 and 7.49 GiB of 8.

## 6. TODO for a Phase 2U, in build order

- [ ] U.1 (S) Own-table route; its D-entry.
- [ ] U.2 (S) Template (G.3).
- [ ] U.3 (S) Trace reader, gather cost table (G.4).
- [ ] U.4 (M) Site files; laptop lines and the MD5SUMS check in the launchers; head memory; `tutorial_test_run.py` retired; TODO.md:158 points at the launcher.
- [ ] U.5 (S) `Code` result (G.2).
- [ ] U.6 (M) `ctx run` and `--check`; Track B budget measured.
- [ ] U.7 (S) Kit round trip; README.md section; course.toml:61.
- [ ] U.8 (S) `conf/datasets/four_runs.toml`, its own table, `[guard]` and `ena_study` (test_repo.py:73,126 read them).
- [ ] U.9 (M) `ctx guide`, its content, START link.
- [ ] U.10 (M) P2 on this Mac; the floor on the page.
- [ ] U.11 (M, Rytis) P3 (G.5).
- [ ] U.12 (S) G.6: CONVENTIONS.md:31 gains both routes; items closed; this plan into `docs/plans/` (lessons.md:73).
- [ ] U.13 (S) G.7: 1280 and 375 px, both themes, every Copy clicked.
- [ ] U.14 (L, calendar) P4 (G.8).
