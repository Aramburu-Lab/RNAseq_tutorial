# CONVENTIONS.md — Cancer Transcriptomics Teaching Suite

Project rules for agents and humans in this repository, imported by `CLAUDE.md`. The repo is scaffolded from `codebase_template`, whose rules are in `AGENTS.md`. Where `CONVENTIONS.md` conflicts with `AGENTS.md`, `CONVENTIONS.md` wins; each divergence is a D-entry in `knowledge/decisions.md`.

## Read first, every session

1. `knowledge/lessons.md` — rules written after corrections; apply them.
2. `TODO.md` — the plan; work only on unchecked items of the current phase.
3. This file.

## Workflow

- Plan mode for any task with 3+ steps or an architectural choice. Write the plan (checkable items) into `TODO.md`, check in before implementing, mark items done as you go, add a short review section when a phase completes.
- Verify before done: run tests, execute the CLI on fixture data, open the HTML offline, measure wall time and peak RAM (the ctx reproduction block; for other commands `/usr/bin/time -l` on macOS, `-v` on Linux). Never report a module complete without this evidence.
- After any correction from the user, add a rule to `knowledge/lessons.md` in the same session.
- Use subagents for research, dataset scouting, and parallel exploration; one task per subagent.
- If something goes sideways, stop and re-plan; do not push through.
- Doc edits: targeted string replacements with grep verification, never whole-file rewrites.
- Never fabricate accessions, versions, or URLs: search, then record with the source, or flag as unverified.
- No access tokens or PATs in the repo or in chat.

## Stack

- Python 3.12; environment via `pixi` (channels: conda-forge, bioconda); containers via multi-stage Docker from the template, Apptainer SIF for HPC.
- CLI package `ctx` (typer). Layout and module list: see `TODO.md` → Repository. One subcommand per module: `ctx setup | gather | qc | process | quant | de | function | proteo | ncrna | report | lite-build`.
- Runtime profiles `laptop` (default) / `hpc` / `lite`; thread and memory budgets come from the profile, command syntax does not change.
- Budgets: Track A < 15 min and < 6 GB peak RAM on 4 cores; Track B < 10 min and < 4 GB. Every report states the measured values.
- Tables: Polars, lazy where possible, chunked concat, downcast to Int32/Float32 only after checking the value range, print `estimated_size()` for every materialised frame. Pandas only inside JupyterLite notebooks.
- Raw reads enter only through nf-core/fetchngs, pinned `-r 1.12.0` (TODO.md D18) and run with `--nf_core_pipeline rnaseq --monochrome_logs true` into a fresh outdir (on Kebnekaise, `monochrome_logs = true` goes in `<project_dir>/nfcore_fetchngs/nextflow.config`: the lab launcher takes no extra flags and starts Nextflow in that directory, whose `nextflow.config` Nextflow reads). Inside a ctx root that outdir is `data/raw/<slug>/`, never `./results/`; the Kebnekaise lab script's `nfcore_fetchngs/results/` is copied there unchanged. `ctx gather` is the only ctx module that reads that layout (`samplesheet/`, `fastq/`, `fastq/md5/`, `metadata/`, `pipeline_info/`), and it never writes into it; every later module reads the harmonised sample sheet in `outputs/01_gather/`. Never write custom ENA/SRA download code (D-006).
- Nextflow comes from the pixi `fetch` environment on a laptop (pinned `==`, both platforms; `-profile docker,arm` on Apple silicon) and from the shared HPC2N install through the lab's `SHARED_run_fetchngs.sbatch` on Kebnekaise. fetchngs 1.12.0 on Nextflow 26 needs `NXF_SYNTAX_PARSER=v1`.
- Bioinformatics tools (bioconda): fastp, seqkit, MultiQC, STAR, HISAT2 (optional), samtools, subread/featureCounts, Salmon. Python: pydeseq2, pytximport, gseapy, decoupler (optional), pysradb, GEOparse, cptac (verify), scikit-bio, scipy, statsmodels, plotly, seaborn, matplotlib, jinja2, typer, pytest, ruff, mypy, interrogate.
- Pin every version in `pixi.toml`; record data-release identifiers (GDC release, GENCODE, MSigDB, TargetScan, miRTarBase, PDC) in the `_README.md` of each output.

## Code rules

- Google-style docstrings with type hints on every function and class; `interrogate` at 100 % for `src/`.
- `ruff` clean, `mypy` clean on `src/`.
- Tests in `tests/` run on chr22-sized fixtures in < 5 min total; every module has an end-to-end test.
- Elegance check before presenting non-trivial code: "is there a simpler way; am I reinventing a wheel?"

## Outputs

- File naming: `YYYY-MM-DD_<module>_<description>_v<X>.<ext>`; underscores, no spaces. A lint (`ctx.io.naming`) enforces this.
- Every data output gets a companion `_README.md`: creation timestamp, source script and commit, tool and data versions, and a data dictionary (one row per column: name, type, unit, meaning).
- `data/` and `resources/` are never committed; `ctx setup` fetches resources with checksums. External inputs (fetchngs outdirs, GDC and Xena downloads) sit under `data/raw/` as fetched; ctx writes only under `<root>/outputs/NN_<module>/`.

## Figures

- Every figure exists in four variants: PDF light and PDF dark (Seaborn/Matplotlib, vector) and HTML with a light/dark toggle (Plotly). One data contract per figure type; `static()` and `interactive()` render from the same object.
- Style: Okabe–Ito palette in this order — #E69F00, #56B4E9, #009E73, #F0E442, #0072B2, #D55E00, #CC79A7, #000000; Helvetica with a sans-serif fallback; 14 pt base; `sns.set_context("paper")`, `sns.despine()`, `plt.tight_layout()`; Nature-style bold lowercase panel letters; no chart junk.
- Plotly: `plotly_white` base plus a custom dark template; hover shows every variable a reader would ask about; exact p-values and n in the caption printed once under each interactive figure, not inside the plot (D-011); the PDFs keep theirs in the figure.
- Statistical annotation is mandatory: test name, n per group, adjusted p (BH) or FDR, effect size, and the caption block written by `ctx.viz.annotate`. A figure without its test is not finished.

## HTML reports

- Self-contained single files: `plotly.min.js` embedded, no CDN, no external requests except optional gene links opened on click. Size < 15 MB; heatmaps capped at 100 genes; embedded tables capped to tested genes with ≤ 8 columns.
- Theme toggle via CSS variables and `Plotly.relayout`; Helvetica 18 px base, the figures' 14 pt plus 3 pt for reading on screen, the interactive plots raised by the same 3 pt (`ctx.viz.style.PAGE_BOOST_PX`, Rytis 2026-09-28); Okabe–Ito.
- Report skeleton, in this order: big question → concept in plain English → parameter cards (parameter, plain-English meaning, default, what changes if you move it) → annotated results → what it means in practice → think-about-it prompts and `<details>` self-check answers → common mistakes → reproduction block with exact command, wall time, peak RAM, versions.
- Every parameter exposed in the CLI must have a parameter card; a test reads CLI defaults and fails if a card is missing or stale.

## Statistics rules

- Multiple testing: BH everywhere, stated in the figure. Report effect sizes next to p-values.
- DE: pre-filter by expression (edgeR-style rule) before pydeseq2; paired or nested designs when the data have them; show the p-value histogram and dispersion plot as health checks; shrunken and raw LFC both available.
- ORA uses an explicit background of tested genes; GSEA uses the Wald statistic, 1000 permutations, seed 42, min size 15, max size 500.
- Proteomics: missing values are not zeros; restrict tests to proteins quantified in ≥ 70 % of samples; state TMT batch design.
- Non-coding: stratify by biotype; never claim miRNA quantities from standard RNA-seq; predicted and validated targets shown separately.
- Cross-study comparisons are labelled as such with batch diagnostics beside them.

## Data and scope boundaries

- TCGA raw reads are controlled-access; only processed expression from GDC is used. Do not attempt dbGaP data.
- Tier-1 dataset: E-MTAB-11190, ENA PRJEB49000 (confirmed 2026-09-27, `TODO.md` D5; read length is confounded with the cancer groups and is taught, never hidden); fetchngs takes it as its ENA study accession (ERP or PRJEB), because it rejects `E-MTAB-*`. Cohort identifiers never go in `src/`, `tests/` (except the D20 fixture of nf-core's own test IDs), `scripts/`, `pyproject.toml` or `pixi.toml` (`tests/test_repo.py::test_no_cohort_identifiers`). Teaching subset and indices are prepared once on HPC and hosted on Zenodo; students never download full FASTQs.
- Licences: code MIT, materials CC BY 4.0; respect KEGG, COSMIC, DisGeNET and miRTarBase terms (documented in `docs/`).

## Definition of done (per module)

- `ctx <module>` runs end to end on fixtures in CI; HTML opens offline; all figure variants present with annotations; every output has `_README.md` and passes the naming lint; notebook executes headless; budgets measured and within limits; biology sanity checks from `TODO.md` pass; `TODO.md` updated.
