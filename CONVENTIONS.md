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
- CLI package `ctx` (typer). Layout and module list: see `TODO.md` → Repository. One subcommand per module: `ctx setup | gather | qc | process | quant | de | function | proteo | ncrna | genecard | clinical | report | lite-build`, plus, from Phase 2D, `ctx run | figure | source | verify-figure` and the global `--cohort`.
- Cancer type is a first-class parameter (Phase 2D): `--cohort` selects a bundle in `conf/cohorts/<cancer>.toml`, which lists the cohort's dataset files per Step and module with its laptop caps and teaching notes; contrasts, covariates, sanity rules and worked-example genes stay in the dataset files. No cohort-specific logic in module code.
- Runtime profiles `laptop` (default) / `hpc`, and `lite` once the Phase 5 JupyterLite site exists; thread and memory budgets come from the profile, command syntax does not change.
- Budgets: Track A < 15 min and < 6 GB peak RAM on 4 cores; Track B < 10 min and < 4 GB. Every report states the measured values.
- Tables: Polars, lazy where possible, chunked concat, downcast to Int32/Float32 only after checking the value range, print `estimated_size()` for every materialised frame. Pandas only inside JupyterLite notebooks at the pydeseq2 boundary in `ctx.de.deseq`, which pydeseq2 requires, and at the lifelines boundary in `ctx.clinical` (M10), whose fitters take a pandas DataFrame.
- Raw reads enter only through nf-core/fetchngs, pinned `-r 1.12.0` (TODO.md D18) and run with `--nf_core_pipeline rnaseq --monochrome_logs true` into a fresh outdir (on Kebnekaise, `monochrome_logs = true` goes in `<project_dir>/nfcore_fetchngs/nextflow.config`: the lab launcher takes no extra flags and starts Nextflow in that directory, whose `nextflow.config` Nextflow reads). Inside a ctx root that outdir is `data/raw/<slug>/`, never `./results/`; the Kebnekaise lab script's `nfcore_fetchngs/results/` is copied there unchanged. `ctx gather` is the only ctx module that reads that layout (`samplesheet/`, `fastq/`, `fastq/md5/`, `metadata/`, `pipeline_info/`), and it never writes into it; `ctx qc` and Track A open the FASTQ files the harmonised sheet names, relative to `--fetchngs-dir`, and nothing else of it (D-031); gather parses the samplesheet and metadata tables by column name, never by position or row order, because fetchngs fixes neither the `metadata/*.runinfo_ftp.tsv` row order nor the samplesheet sort key (1.12.0 sorts per-sample files by name, 1.13.0 sorts rows as text; D-006); every later module reads the harmonised sample sheet in `outputs/01_gather/`. Never write custom ENA/SRA download code (D-006).
- Nextflow comes from the pixi `fetch` environment on a laptop (pinned `==`, both platforms; `-profile docker,arm` on Apple silicon) and from the shared HPC2N install through the lab's `SHARED_run_fetchngs.sbatch` on Kebnekaise. fetchngs 1.12.0 on Nextflow 26 needs `NXF_SYNTAX_PARSER=v1`.
- Step-1 counts come from nf-core/rnaseq -r 3.27.0 (`ctx.gather.rnaseq.RNASEQ_VERSION`; TODO.md D21), Salmon only (`pseudo_aligner: salmon`, `skip_alignment: true`), fed by the fetchngs samplesheet cut to `sample, fastq_1, fastq_2, strandedness`, under `NXF_SYNTAX_PARSER=v2`; `ctx gather` refuses any other release, and the run's MultiQC stays on the cluster.
- Bioinformatics tools (bioconda), in the `reads` pixi environment outside the analysis solve and called by subprocess through `ctx.tools`, which refuses any version other than the pin (D-027): FastQC, fastp, seqkit, MultiQC, STAR, HISAT2 (optional), samtools, subread/featureCounts, Salmon, Trim Galore; seqtk for the instructor-side teaching subset. Python: pydeseq2, pytximport, gseapy, decoupler (optional), pysradb, GEOparse, cptac (verify), scikit-bio, scipy, statsmodels, lifelines (M10), plotly, seaborn, matplotlib, jinja2, typer, pytest, ruff, mypy, interrogate.
- Pin every version in `pixi.toml`; record data-release identifiers (GDC release, GENCODE, MSigDB, TargetScan, miRTarBase, PDC) in the `_README.md` of each output.

## Code rules

- Google-style docstrings with type hints on every function and class; `interrogate` at 100 % for `src/`.
- `ruff` clean, `mypy` clean on `src/`.
- Tests in `tests/` run on chr22-sized fixtures in < 5 min total; every module has an end-to-end test.
- Elegance check before presenting non-trivial code: "is there a simpler way; am I reinventing a wheel?"

## Outputs

- File naming: `YYYY-MM-DD_<module>_<description>_v<X>.<ext>`; underscores, no spaces. A lint (`ctx.io.naming`) enforces this.
- Every data output gets a companion `_README.md`: creation timestamp, source script and commit, tool and data versions, and a data dictionary (one row per column: name, type, unit, meaning).
- `data/` and `resources/` are never committed; `ctx setup` fetches resources with checksums; `ctx genecard` alone writes `resources/cache/`, with its own manifest (TODO.md M9). External inputs (fetchngs outdirs, GDC and Xena downloads) sit under `data/raw/` as fetched. Writers: modules write only under `<root>/outputs/NN_<module>/`, `ctx setup` under `resources/`, `ctx genecard` under `resources/cache/`; from Phase 2D, `ctx figure` run from a Reproduce panel writes into the reader's own working folder, outside any ctx root.

## Figures

- Every figure exists as `_white` and `_black` PDF and SVG from one render (Seaborn/Matplotlib, vector), one HTML with a light/dark toggle (Plotly), and its numbers as CSV with `_README.md` (TODO.md D16). One data contract per figure type; `static()` and `interactive()` render from the same object.
- Style: Okabe–Ito palette in this order — #E69F00, #56B4E9, #009E73, #F0E442, #0072B2, #D55E00, #CC79A7, #000000. Sizes and parameters of the Candida Ribo-seq figures (Rytis 2026-09-29, D-024), set in `ctx.viz.style` alone: Arial, Liberation Sans, DejaVu Sans; 12 pt panel letters (bold capitals), 10 pt axis titles, 9 pt ticks, legends, annotations and captions, nothing below 8 pt; 1 pt axes and lines, 1.5 pt boxes and error bars; 3 pt ticks, outward; 6 pt points at alpha 0.5; no seaborn context; `sns.despine()`, `plt.tight_layout()`; no chart junk.
- Plotly: `plotly_white` base plus a custom dark template; hover shows every variable a reader would ask about; exact p-values and n in the caption printed once under each interactive figure, not inside the plot (D-011); the PDFs keep theirs in the figure.
- Continuous values (filled contours, colour bars) use cividis in both themes, with isolines and their labels in ink; Okabe–Ito stays for categories, including per-group contour lines. Contour plots only for a continuous value over two continuous axes (computed surfaces, or densities from ≥ 30 points per group), never for categorical axes; density levels are 50 % / 90 % highest-density regions from one grid shared by both renderers.
- Statistical annotation is mandatory: test name, n per group, adjusted p (BH) or FDR, effect size, and the caption block written by `ctx.viz.annotate`. A figure without its test is not finished.
- External summary data (GTEx, HPA) get descriptive annotation only (n, median, IQR, tau, release); never inferential tests on external medians. Descriptive figures are the exception to the mandatory test above: M1's counts, the QC bands and external summaries say "no test" in their caption, and an external panel's caption also names the release and retrieval date.
- From Phase 2D, every figure has a stable ID (`<cohort>_<module>_<figure>`) and a parameter model as a frozen dataclass with one loader that reports every error together, and is produced only through `ctx figure <id>`. Its Reproduce panel is generated from the `FigureSpec` (environment, inputs with checksums, params heredoc, generate naming the plotting file, check; upstream collapsed); never hand-write panel commands. Pasted blocks contain no `set -e`; the downloadable script uses `set -euo pipefail`.
- Fonts: PDFs embed fonts as Type 42; from Phase 2D the image ships Liberation Sans, Arial's redistributable metric twin, and figure files carry no creation timestamp.

## HTML reports

- Self-contained single files: `plotly.min.js` embedded, no CDN, no external requests except optional gene links opened on click. Size < 15 MB; heatmaps capped at 100 genes; embedded tables capped to tested genes with ≤ 8 columns.
- Design of the ALI report (D-024): a frosted sticky top bar, an eyebrow above the title, numbered sections under hairlines, white cards with a soft shadow on a near-white page, stat tiles, rounded table wraps, tinted callouts for practice and mistakes; a skip link and focus rings, which ALI lacks; every text colour at 4.5:1 or more on each surface it sits on (`tests/test_report.py::test_text_contrast`). Theme toggle via CSS variables and `Plotly.relayout`; an 18 px page base; the interactive plots at the figure sizes plus 3 pt, converted to px at 96/72 (`ctx.viz.style.on_page`, Rytis 2026-09-28); Okabe–Ito.
- Report skeleton, in this order: big question → concept in plain English → parameter cards (parameter, plain-English meaning, default, what changes if you move it) → annotated results → what it means in practice → think-about-it prompts and `<details>` self-check answers → common mistakes → reproduction block with exact command, wall time, peak RAM, versions.
- Page text lives in `src/ctx/report/content/*.toml`, one file per page plus `ui.toml` and `cards.toml`, read through `ctx.report.text`; code names a key and fills its `{placeholders}`, never writes a sentence (D-024, Rytis 2026-09-29). Setup and gather do so now; M5, M6, the hub and the captions move after the Step-1 record (P2S.4). An unfilled placeholder or an unread key fails a test.
- Every parameter exposed in the CLI must have a parameter card; a test reads CLI defaults and fails if a card is missing or stale.

## Statistics rules

- Multiple testing: BH everywhere, stated in the figure. Report effect sizes next to p-values.
- DE: pre-filter by expression (edgeR-style rule) before pydeseq2; paired or nested designs when the data have them; show the p-value histogram and dispersion plot as health checks; shrunken and raw LFC both available, shrinkage on unless `ctx de --no-lfc-shrink` turns it off (D-019).
- ORA uses an explicit background of tested genes; GSEA uses the Wald statistic, 1000 permutations, seed 42, min size 15, max size 500.
- Proteomics: missing values are not zeros; restrict tests to proteins quantified in ≥ 70 % of samples; state TMT batch design.
- Non-coding: stratify by biotype; never claim miRNA quantities from standard RNA-seq; predicted and validated targets shown separately.
- Cross-study comparisons are labelled as such with batch diagnostics beside them.
- Survival: the endpoint is chosen per TCGA-CDR's recommendations (Liu et al. 2018, Table 3) and named in every figure; continuous Cox is the primary model; Kaplan–Meier splits are pre-specified, and any data-driven cut-off carries the corrected p of the maximally selected statistic; the PH assumption is tested and reported; every KM plot has a risk table and censor marks.
- Treatment: TCGA tumours were sampled before systemic therapy for them (exclusionary on TCGA's case quality-control forms), with exceptions GDC notifies: neoadjuvant therapy (42 patients in the CDR, 102 GDC notifications) and treatment for a prior malignancy, so expression by treatment is never presented as a treatment effect; treatment effects come only from pre/on-treatment datasets (TODO.md D28).
- QC bands live only in `conf/qc_metrics.toml`; every band has a source and a source type, course conventions are labelled as such, and no QC threshold is hard-coded elsewhere, except the sex check in `ctx.de.run.sex_check`, frozen with M5 until the Step-1 record; after it the check moves to `ctx.qc.identity` and M5 imports it from there. MultiQC, pinned in the `reads` environment (D-027), is the only parser of tool QC logs, read through `ctx.qc.multiqc_io`.

## Data and scope boundaries

- TCGA raw reads are controlled-access; only processed expression from GDC is used. Do not attempt dbGaP data.
- Step-1 dataset: E-MTAB-11190, ENA PRJEB49000 (confirmed 2026-09-27, `TODO.md` D5; read length is confounded with the cancer groups and is taught, never hidden); fetchngs takes it as its ENA study accession (ERP or PRJEB), because it rejects `E-MTAB-*`. Cohort identifiers never go in `src/`, `tests/` (except the D20 fixture of nf-core's own test IDs), `scripts/`, `pyproject.toml` or `pixi.toml` (`tests/test_repo.py::test_no_cohort_identifiers`). Teaching subset and indices are prepared once on HPC and hosted on Zenodo; students never download full FASTQs.
- External knowledge bases (GTEx, HPA, Ensembl): documented APIs and bulk files only, cached under `resources/cache/` per release with the retrieval date and the SHA-256 measured on arrival; never scrape, mirror or iframe their pages; cite HPA per its licence page: a primary publication, proteinatlas.org, and the versioned URL `v<NN>.proteinatlas.org`.
- Clinical data: open-access, de-identified only; harmonised vocabulary with explicit missingness; GDC annotations are read and logged, and any exclusion made from them is the course's own rule, stated on the page (the neoadjuvant and prior-malignancy annotations are notifications; only a redaction removes a case).
- Licences: code MIT, materials CC BY 4.0; respect KEGG, COSMIC, DisGeNET and miRTarBase terms (documented in `docs/`). CC BY 4.0 content (HPA, AlphaFold 2) may be embedded with its citation; link, never embed, what HPA lists under user terms (Ensembl, GTEx, TCGA, Tabula Sapiens, IntAct, BioPlex, OpenCell, TFClass, DepMap) when taken from HPA, AlphaFold 3 (CC BY-NC-SA 4.0) and CC BY-SA sources such as IUPHAR, TCDB and Metabolic Atlas; GTEx and Ensembl data fetched from their own APIs follow their own terms: draw GTEx data only as our own figures with the GTEx Portal acknowledgement and retrieval date, its terms of use UNVERIFIED until the external-data spike reads them, and link instead if they forbid it.

## Definition of done (per module)

- `ctx <module>` runs end to end on fixtures in CI; HTML opens offline; all figure variants present with annotations; from Phase 2D every figure has a Reproduce panel that CI executes successfully; every output has `_README.md` and passes the naming lint; notebook executes headless; budgets measured and within limits; biology sanity checks from `TODO.md` pass; `TODO.md` updated.
