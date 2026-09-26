# Cancer Transcriptomics Teaching Suite — plan

Status: **draft v3.1 (2026-09-26) — D1, D2, D4, D11–D13 resolved; remaining Phase 0 items: D5 dataset verification and version pinning.**
Outputs of this plan: a Python package + CLI, one Jupyter notebook per module, and one self-contained interactive HTML report per module, linked from a hub page.

---

## 🎯 Goal and learning outcomes

Students analyse public human cancer RNA-seq (ENA / GEO / TCGA) end to end on a laptop, connect it to public proteomics, and leave able to:

1. Choose the right data tier and repository for a question (raw reads vs processed counts; open vs controlled access).
2. Judge whether reads are usable from QC metrics, and know what trimming/filtering parameters do.
3. Explain what an aligner, a pseudo-aligner and a read counter each produce, and why counts differ between them.
4. Set up a differential-expression (DE) design correctly (pairing, covariates, filtering, multiple testing, shrinkage) and read the diagnostics.
5. Translate a DE gene list into biology (pathways, cellular roles, disease associations) without over-interpreting.
6. Read every figure type in the syllabus (box, violin, heatmap, volcano, MA, PCA, Venn/UpSet, enrichment) with its statistical annotation.
7. Test whether an RNA change is carried through to protein, and name the reasons it may not be (translational control, protein stability, timing of synthesis, complex stoichiometry, non-coding function).
8. State what standard RNA-seq can and cannot say about lncRNAs and miRNAs, and choose the right library type, annotation and filters for each.

Pedagogical device used throughout: **"turn the knob"** — every module exposes 2–4 parameters as interactive controls or precomputed alternatives, so students see how a choice changes the result.

---

## 🏗️ Architecture: two tracks, one downstream

A full human alignment does not fit the RAM/speed constraint (STAR full-genome index needs ~32 GB). The design therefore separates the *read-level concepts* (taught on a laptop-sized subset of raw reads) from the *cohort-level statistics* (taught on count tables). Cohorts follow the agreed ladder: innate immune cells first, other immune cells second, large solid-tumour cohorts with large healthy groups third. For the immune-cell tiers one public dataset serves both tracks (raw FASTQs on ENA for Track A, the full cohort for Track B), so students follow one biology from reads to pathways. TCGA enters at tier 3 as processed counts only, since TCGA FASTQ/BAM are dbGaP-controlled.

```mermaid
flowchart LR
    accTitle: Two-Track Teaching Architecture
    accDescr: Track A runs the complete read-level pipeline on subsampled FASTQ files so it fits a laptop. Track B uses count tables from immune-cell cohorts and, at tier 3, from TCGA and GTEx, for realistic cohort statistics. Both tracks feed the same differential expression and functional annotation modules; two integration modules compare RNA changes with proteomics and handle non-coding RNAs before reporting.

    subgraph track_a["Track A — raw reads, laptop scale"]
        fastq["ENA FASTQ, tier-1 subset<br/>subsampled ~1M pairs/sample"] --> qc["fastp<br/>QC + trimming"]
        qc --> align["STAR on chr22 or HISAT2<br/>BAM + featureCounts"]
        qc --> salmon["Salmon<br/>full transcriptome"]
        salmon --> tximport["pytximport<br/>gene-level counts"]
    end

    subgraph track_b["Track B — cohorts, count tables"]
        gdc["ArrayExpress/GEO counts, tiers 1–2<br/>GDC + GTEx/TOIL, tier 3"] --> matrix["Polars count matrix<br/>+ harmonised sample sheet"]
    end

    tximport --> de["pydeseq2 DE<br/>PCA, volcano, MA, heatmap,<br/>box/violin, Venn/UpSet"]
    matrix --> de
    align -.->|"compare count methods"| de
    de --> func["gseapy + Open Targets<br/>ORA, GSEA, disease association"]
    de --> prot["CPTAC / PDC proteomics<br/>RNA–protein concordance"]
    de --> ncrna["Biotype-aware DE<br/>lncRNA, miRNA-seq targets"]
    func --> report["Self-contained HTML reports<br/>light/dark, sliders, quizzes"]
    prot --> report
    ncrna --> report

    classDef reads fill:#dbeafe,stroke:#2563eb,stroke-width:2px,color:#1e3a5f
    classDef cohort fill:#fef9c3,stroke:#ca8a04,stroke-width:2px,color:#713f12
    classDef shared fill:#dcfce7,stroke:#16a34a,stroke-width:2px,color:#14532d

    class fastq,qc,align,salmon,tximport reads
    class gdc,matrix cohort
    class de,func,prot,ncrna,report shared
```

| Concern | Track A (raw reads) | Track B (count tables) |
| --- | --- | --- |
| Teaches | Modules 1–4 (fetch, QC, trim/align, quantify) | Modules 1, 5–8 (fetch, DE, function, RNA→protein, non-coding RNA) at realistic n |
| Data size | ~6–8 samples × 1M read pairs (~1 GB total) | tier 1 ≈ 200 samples (tens of MB); tier 3 ≈ 150–1000 samples (< 1 GB) |
| Peak RAM target | < 6 GB | < 4 GB |
| Wall time target | < 15 min on 4 cores | < 10 min on 4 cores |
| Why | Real tools, real parameters, visible effects | Real biology, real statistics, nested and paired designs |

Bridges between tracks: for tiers 1–2 the dataset itself — Track A quantifies a subset of the very samples whose full cohort Track B analyses, and the two count vectors are compared directly; for tier 3, STAR `--quantMode GeneCounts` produces the same four-column table (unstranded / stranded-first / stranded-second) that GDC exposes in its STAR-Counts files, so students see where TCGA numbers come from.

---

## 📦 Datasets

Cohort ladder agreed at check-in: innate immune cells → other immune cells → large solid-tumour cohorts with large healthy groups. Criteria for every candidate: human, bulk RNA-seq, cancer vs healthy groups with n ≥ 10 per group, raw FASTQs open on ENA (needed for Track A), rich sample metadata, permissive licence.

| Tier | Candidate | Design | Why | Verify before use |
| --- | --- | --- | --- | --- |
| 1 — innate immune cells (primary) | ArrayExpress E-MTAB-11190 (released 2022-11): bulk RNA-seq of normal-density neutrophils, low-density neutrophils and monocytes from blood; healthy controls n = 19, G-CSF-treated donors n = 17, HSCT recipients n = 8, PDAC n = 15, IPMN n = 14; 210 samples from 73 individuals; reanalysed by Caronni et al. 2023 Nature | Cell type × condition, nested within donor | Neutrophils and monocytes in cancer vs healthy at usable n; stress-induced myelopoiesis story (G-CSF, LDN vs NDN); one dataset serves both tracks | ENA project accession and FASTQ availability, library layout and read length, SDRF fields (cell type, condition, donor, batch), processed counts provided or not, source publication and its positive-control genes, licence |
| 1 — alternative | Cassetta et al. 2019 Cancer Cell: circulating monocytes from breast / endometrial cancer patients vs healthy donors; TAMs vs tissue-resident macrophages | Cell type × cancer status | Systemic and local myeloid reprogramming; complements the neutrophil dataset | GEO accession, FASTQs on ENA, group sizes |
| 2 — other immune cells | TCGA-LAML (myeloid leukaemia, ~150 cases) or TCGA-DLBC (B-cell lymphoma, ~48 cases) plus a healthy sorted-cell comparator from GEO (bone-marrow CD34+ / B cells) | Cancer vs healthy across studies | Teaches the cross-study batch problem explicitly: same lineage, different labs | Comparator accession, protocol match (poly-A vs total RNA), current TCGA case counts |
| 3 — solid tumours, large cohorts | TCGA tumour vs solid-tissue-normal (KIRC ~70 pairs, BRCA ~110 pairs, LUAD ~55 pairs) and TCGA vs GTEx through the UCSC Xena TOIL recompute (uniformly processed, hundreds of healthy samples per tissue) | Paired (TCGA) and unpaired large-n (TOIL) | Large groups; VHL→HIF→hypoxia story in KIRC; harmonisation and batch lessons in the TOIL comparison | Pair counts in current GDC release; current TOIL file names and gene model |
| Proteomics, tier 3 | CPTAC-3 ccRCC (Clark et al. 2019 Cell): tumour and normal-adjacent tissue with matched RNA-seq (open at GDC) and TMT proteomics (open at the Proteomic Data Commons); CPTAC LUAD / PDAC as alternatives | Same samples measured at RNA and protein level | Direct RNA–protein concordance on the tier-3 cancer; PDAC option echoes tier 1 | Case counts with both layers; PDC / `cptac` access route; TMT batch layout |
| Proteomics, tier 1 | Rieckmann et al. 2017 Nat Immunol immune-cell proteome (PRIDE; healthy donors; protein copy numbers per cell type incl. neutrophils and monocytes) | Cross-dataset: RNA from tier 1 vs protein per cell type | Neutrophil vs monocyte RNA–protein discordance (granule proteins are transcribed early in granulopoiesis, so mature neutrophils carry high protein with low mRNA) | PRIDE accession; processed copy-number table availability; cell-type labels matching tier 1 |
| miRNA-seq, tier 3 | TCGA-KIRC miRNA-seq (mature-miRNA and isoform quantification, open at GDC) for the same cases as mRNA-seq | Matched mRNA + miRNA per sample | miRNA–target anticorrelation and its protein footprint in CPTAC | Number of normals with miRNA-seq; GDC file types and normalisation |
| ncRNA resources | GENCODE biotypes; TargetScan predicted targets (context++ scores); miRTarBase validated targets; LncBook / Lnc2Cancer for disease-linked lncRNAs; GENCODE Ribo-seq ORF annotation to flag "non-coding" genes with translation evidence | — | Annotation is a status, not a proof | Download terms (miRTarBase academic licence); current versions |
| Track A subset | 3–4 donors × 2–3 cell types from tier 1, subsampled (`seqtk sample -s42`, 1M pairs/sample) on Kebnekaise, hosted on Zenodo with DOI + md5 | — | Students never download tens of GB; same biology as Track B | Zenodo quota, redistribution terms |
| References | GENCODE primary assembly + GTF (current release), chr22 slice for STAR, transcript FASTA for Salmon, tx2gene from GTF | — | Versions pinned in every `_README.md` | Current GENCODE release; GENCODE version behind current GDC STAR-Counts |
| Gene sets | MSigDB Hallmark, GO BP, Reactome (GMTs cached by `ctx setup`); Enrichr disease libraries (DisGeNET, Jensen DISEASES, GWAS Catalog) as GMTs; Open Targets GraphQL for per-gene disease scores; NCG cancer-gene list; immune sets (MSigDB C7 ImmuneSigDB, neutrophil/monocyte state signatures from the source paper) | — | Offline-capable, permissive; immune sets make tier-1 results interpretable | Current MSigDB version string; Enrichr library names |

---

## 🧩 Module specification

Each module = one CLI subcommand (`ctx <module>`), one notebook (`notebooks/0N_<module>.ipynb`, student + solution variants), one HTML report. Every report follows the same skeleton: **big question → concept in plain English → parameter cards → annotated results → what it means in practice → think-about-it prompts / self-check quiz → exact reproduction command with runtime and peak RAM**.

### M0 — Setup

- `ctx setup`: creates `resources/` (references, indices, GMTs, teaching subset) with checksums; prints environment report (versions, CPU, RAM).
- Environment: pixi (conda-forge + bioconda) → Docker (multi-stage, from `codebase_template`) → Apptainer SIF.
- Knob: none. Output: `resources/manifest.json` + README.

### M1 — Gather public data

- Tools: BioStudies/ArrayExpress API (SDRF + processed files) and ENA Portal API `filereport` for FASTQ URLs/md5 (tier 1); `pysradb` for GEO↔SRA mapping and `GEOparse` for GSE metadata (tier 2); GDC REST API (`files`, `cases`, `data`) via `requests` and UCSC Xena TOIL downloads (tier 3).
- Teaches: repositories and their data tiers (raw vs processed; open vs controlled), how one study appears in ArrayExpress, ENA and GEO/SRA, sample-type codes (01 tumour / 11 normal), one aliquot per sample, pairing by donor or case ID, metadata hygiene, cross-study batch as a design constraint.
- Knob: cell-type and condition filters (tier 1), cancer project and sample type (tier 3), minimum group size.
- Figures: cohort composition bars; metadata completeness heatmap; library-size bars (from `N_*` summary rows).
- Stats: none inferential; descriptive tables.
- Output: `data/raw/<dataset>/` + harmonised sample sheet (CSV: sample, donor, cell_type, condition, batch) + `_README.md`; Polars matrix assembled lazily (`scan_csv`, comment rows dropped, count column selected, Int32 downcast, `estimated_size()` printed).

### M2 — Read QC

- Tools: `fastp` (JSON parsed into Plotly), `seqkit stats`, MultiQC for the classic aggregate.
- Teaches: per-position quality, GC content, adapter content, duplication, insert size, Q30 rate; what "unusable" looks like (planted bad sample: adapter-heavy or low-Q30 subset); why neutrophil libraries are special (low RNA content, granule-protein and ribosomal signal, degradation), read from the real tier-1 data.
- Knob: quality threshold `-q`, unqualified-base % `-u`, min length `-l`, `--cut_tail` window/mean.
- Figures: per-position quality curves with threshold line; GC/duplication distributions; before/after read-count Sankey (Mermaid in notebook, Plotly in report).
- Stats: descriptive; flag rules explained (why Q30 < 80 % is a warning, not a verdict).

### M3 — Trim, filter, align, count

- Tools: `fastp` (trim), `STAR` on chr22 index (`--genomeSAindexNbases 11`, `--quantMode GeneCounts`, `--outSAMtype BAM SortedByCoordinate`), optional `HISAT2` full-genome low-RAM path, `samtools flagstat/idxstats`, `featureCounts`.
- Teaches: splice-aware alignment, multimappers (`--outFilterMultimapNmax`), strandedness (`-s 0/1/2`; inferred from Salmon `lib_format_counts.json`), fragment vs read counting (`-p --countReadPairs`), overlap handling (`-O`, `--fraction`), incomplete reference as an error source (chr22-only vs full → mapping-rate drop).
- Knob: strandedness setting, multimapper limit, chr22 vs full reference.
- Figures: mapping-rate bars per setting; strandedness column comparison; STAR-GeneCounts vs featureCounts scatter with concordance stats.
- Stats: Spearman/Pearson concordance with CI; Bland–Altman style bias plot.

### M4 — Transcript/gene quantification

- Tools: `Salmon` (`-l A`, `--gcBias`, `--seqBias`, decoy-aware index if hosted, else no-decoy with the trade-off explained), `pytximport` (`counts_from_abundance` options).
- Teaches: TPM vs counts vs CPM; effective length; why DE wants counts; why pydeseq2 receives `length_scaled_tpm` counts (no per-sample offset matrix support — stated honestly); mapping rate vs assignment rate.
- Knob: `counts_from_abundance`, bias correction on/off.
- Figures: library-size and count-distribution box/violin; gene-detection curves; sample–sample correlation heatmap with dendrogram; Salmon vs STAR vs featureCounts Venn of detected genes.
- Stats: Spearman correlation matrix; hierarchical clustering (average linkage, optimal leaf ordering).

### M5 — Differential expression

- Tools: `pydeseq2` (design formulas `~ donor + condition` for paired tissue, `~ cell_type + condition` or `~ cell_type * condition` for tier 1; `DeseqStats` with BH, independent filtering, Cook's; `lfc_shrink`; `lfc_null` threshold test; `vst()` for PCA/heatmap), `scipy.stats`, `scikit-bio` PERMANOVA.
- Teaches: pre-filtering (edgeR-style "≥10 counts in ≥ smallest group"), size factors (median-of-ratios), dispersion trend/shrinkage, Wald test, FDR in plain English, NA `padj` meaning, effect size vs significance, shrinkage, paired vs unpaired design, cell type as the dominant axis with disease nested inside it (contrasts and interaction terms), p-value histogram as a health check, batch inspection (sequencing batch in tier 1; TCGA barcode TSS/plate in tier 3).
- Knob: design (paired / unpaired; with / without cell type), contrast selector, `padj` and |LFC| sliders, shrunken vs raw LFC toggle, `lfc_null` threshold, min-count filter (precomputed grid).
- Figures: PCA (VST, top-500 variable genes; PERMANOVA pseudo-F and p annotated; colour = condition, symbol = cell type or batch); dispersion plot; p-value histogram; volcano (sliders, top-gene labels, up/down counts); MA (raw vs shrunk); heatmap of top 50 (z-scored VST, clustered, group colour bar); box + violin of top genes with paired lines, Wilcoxon signed-rank (paired) or Mann–Whitney, BH across plotted genes, rank-biserial effect size, DESeq2 `padj` shown alongside with a note on why they differ; Venn/UpSet: paired vs unpaired, shrunk vs raw thresholds, neutrophil vs monocyte cancer responses (tier 1), KIRC vs BRCA vs LUAD (tier 3).
- Output: results table (Parquet + CSV) with data dictionary; `_README.md` records design, versions, filters.

### M6 — Function and disease association

- Tools: `gseapy` (`enrich` with explicit background = tested genes; `prerank` on Wald stat, 1000 permutations, seed 42, `min_size 15`, `max_size 500`), ImmuneSigDB and source-paper signatures as extra GMTs, optional `decoupler` PROGENy pathway activity, Open Targets GraphQL, NCG list; `ctx.function.lite` pure-NumPy ORA and prerank-GSEA (validated against gseapy in tests) for the browser.
- Teaches: ORA vs GSEA; the background trap (same list, three backgrounds → inflated results demo); ranking metric choice; term redundancy; "association ≠ causation" for disease sets; how tier-1 results map onto emergency granulopoiesis and immunosuppressive myeloid states, and how KIRC hypoxia/glycolysis output maps onto VHL biology and therapy.
- Knob: background choice, ranking metric, DE thresholds feeding ORA (live, client-side hypergeometric on embedded Hallmark + Reactome).
- Figures: dot/bar plots (−log10 FDR, gene ratio, count); GSEA running-score plots; enrichment-map network (Jaccard-linked terms); disease-association bars (Open Targets score, DisGeNET ORA); gene-to-term heatmap.
- Stats: hypergeometric/Fisher with BH; GSEA permutation FDR; term-similarity Jaccard.

### M7 — Is the transcript a template? RNA versus protein

- Tools: `cptac` package or PDC GraphQL API via `requests` for CPTAC protein tables; GDC RNA-seq for the same cases; `scipy.stats` (Spearman, Welch t, Kolmogorov–Smirnov); `statsmodels` (BH); optional `inmoose` limma port for moderated t-tests (verify).
- Teaches: what MS proteomics measures (peptides → protein groups; TMT ratios; missing values are not zeros; ratio compression); joining layers by case ID and gene symbol / UniProt mapping (isoforms, many-to-one); three questions — does protein track RNA across samples (per-gene Spearman ρ), does the tumour-vs-normal change carry to protein (log2FC concordance), which genes change at RNA only or protein only; reasons for discordance: translational control, protein half-life, complex stoichiometry (subunits buffered), secretion, timing of synthesis (neutrophil granule proteins), and genuinely non-coding function; Ribo-seq as the missing middle layer (translation efficiency; microproteins/smORFs) — concept only (D12).
- Knob: correlation method, minimum samples with protein quantified, log2FC concordance threshold, gene class (complex subunits / enzymes / secreted), predicted vs validated protein evidence.
- Figures: RNA log2FC vs protein log2FC scatter with quadrant counts and Spearman ρ (bootstrap 95 % CI); violin of per-gene ρ by gene class; concordance MA (mean abundance vs RNA−protein LFC difference); Venn/UpSet of DE at RNA vs protein; heatmap of paired RNA/protein z-scores for top genes; tier 1: neutrophil-vs-monocyte RNA fold change vs protein copy-number fold change with granule proteins highlighted.
- Stats: per-gene Spearman with BH; across-gene concordance ρ; class comparisons Kruskal–Wallis + Dunn with BH and effect sizes; protein DE by Welch t-test on log2 ratios (moderated t if available) with BH; missingness reported per protein and per sample.
- Output: merged RNA–protein table (Parquet + CSV) with data dictionary; per-gene class labels with source.

### M8 — Beyond mRNA: lncRNA and miRNA nuances

- Tools: GENCODE GTF biotypes via `pyranges` / Polars; GDC miRNA-seq quantification (reads per million miRNA-mapped; isoform tables); TargetScan / miRTarBase tables; `gseapy` for target-set enrichment; `pydeseq2` for biotype-stratified DE.
- Teaches: library type decides what is visible — poly-A selection vs rRNA depletion (non-polyadenylated lncRNAs, histone mRNAs, circRNAs) and why mature miRNAs are invisible in standard RNA-seq (size selection; small-RNA-seq with adapter trimming of inserts shorter than the read — explained, not run: D13); strandedness is needed to separate antisense lncRNAs from overlapping mRNAs (TCGA mRNA-seq is unstranded — stated as a limit); lncRNA pitfalls: low counts, length bias in TPM, multi-mapping, and how the min-count filter and independent filtering remove them; no GO terms for lncRNAs — guilt-by-association through correlated mRNAs and cis neighbours, with caveats; miRNA units, isomiRs, family-level multimapping; targets — predicted vs validated, anticorrelation across matched samples, target enrichment among downregulated mRNAs, and why repression is often clearer at protein level (M7); the ceRNA hypothesis presented as contested; "non-coding" is an annotation status — Ribo-seq ORF catalogues show translated smORFs inside annotated lncRNAs.
- Knob: biotype selection, min-count filter per biotype, target source (predicted / validated / both) and context-score threshold, library-type filter where both types exist.
- Figures: biotype composition stacked bars per sample and library type; count distributions per biotype (violin); volcano coloured by biotype; lncRNA–neighbour co-expression scatter; miRNA volcano; miRNA–target anticorrelation scatter; cumulative distribution of target vs non-target mRNA log2FC with KS statistic; miRNA→mRNA→protein tri-layer plot for selected miRNAs.
- Stats: biotype-stratified DE with BH; Spearman anticorrelation with BH; KS and Mann–Whitney for target shifts with median-shift effect size; Fisher test for target enrichment among DE genes with matched background.
- Output: biotype-annotated DE results; miRNA results; target tables with evidence source.

### Hub

- `reports/index.html`: module cards, glossary, "start here" path, teacher notes link; offline; no external requests except optional gene links (GeneCards/NCBI/Open Targets) opened on click.

---

## 📊 Figure and statistics conventions

- Palette: Okabe–Ito; font Helvetica 14 pt base; `sns.set_context("paper")`, `sns.despine()`, `plt.tight_layout()`; panel letters bold lowercase (Nature style); exact p-values, n per group, test name and effect size in-figure.
- Every figure in four variants: PDF light + dark (Seaborn, vector) and HTML with theme toggle (Plotly, `plotly_white` base + custom dark template; hover shows all relevant variables).
- Single implementation per figure type in `ctx.viz` with `static()` and `interactive()` methods sharing one data contract, so the two renderings cannot drift.
- Statistical annotation helper `ctx.viz.annotate` writes the same caption block to both renderers.

| Figure | Test / annotation |
| --- | --- |
| Box / violin | Wilcoxon signed-rank (paired) or Mann–Whitney U; BH across plotted genes; rank-biserial r; n |
| PCA | % variance per PC; PERMANOVA (Euclidean, 999 perms) pseudo-F, p |
| Volcano / MA | Wald `padj` (BH); thresholds shown; up/down counts |
| Heatmap | Linkage method, distance metric, z-score definition in caption |
| Venn / UpSet | Set sizes and intersections; Fisher test for overlap enrichment where two sets are compared |
| Enrichment | Test (hypergeometric / permutation), background size, FDR |
| Concordance scatter | Spearman ρ with 95 % CI (bootstrap) |
| RNA–protein scatter | Spearman ρ with 95 % CI; quadrant counts; n genes; missingness stated |
| Target CDF | Kolmogorov–Smirnov D and p; median log2FC shift; target and background sizes |

---

## 🖥️ Interactive HTML design

- One self-contained file per module (Jinja2 template; `plotly.min.js` embedded; no CDN). Light/dark toggle via CSS variables + `Plotly.relayout` with embedded light/dark layout templates. Helvetica 14 px, Okabe–Ito.
- Interactivity without a server: results tables embedded as JSON; vanilla-JS sliders/toggles call `Plotly.react`; client-side hypergeometric ORA on a small embedded gene-set library; precomputed grids for heavy knobs (design, filters).
- Parameter cards: table of parameter → plain-English meaning → default → what changes if you move it.
- Pedagogy elements: think-about-it prompts, `<details>` self-check answers, "common mistakes" callouts, reproduction block with measured runtime/RAM.
- Size guard: heatmaps capped at 50–100 genes; embedded tables capped to tested genes with ≤ 8 columns.
- In-browser Python (D4, agreed): a JupyterLite site (`jupyter lite build`, static, GitHub Pages) ships notebooks 05–08 with the count tables and DE results bundled as CSV/JSON (< 50 MB total), so students re-run M5–M8 in the browser with nothing installed. Pyodide rules: `pandas`/`numpy`/`scipy`/`plotly` only (no Polars, no compiled bioinformatics tools); DE re-runs use `pydeseq2` only if the Phase 1 spike shows it installs via micropip, otherwise precomputed grids plus live nonparametric tests; enrichment uses `ctx.function.lite`. Track A (M2–M4) stays CLI-only.

---

## 🗂️ Repository and engineering conventions

Scaffold from Aramburu-Lab `codebase_template` (Copier): pixi multi-stage containers, `digests.lock`/ORAS SIFs, HPC2N configs for the instructor-side data preparation.

```text
cancer-transcriptomics-teaching/
├── pixi.toml, pixi.lock, Dockerfile, apptainer.def
├── README.md, TODO.md, knowledge/lessons.md, CITATION.cff, LICENSE (code MIT, materials CC BY 4.0)
├── src/ctx/
│   ├── acquire/   gdc.py ena.py geo.py
│   ├── qc/        fastp_parse.py
│   ├── process/   trim.py align.py count.py strandedness.py
│   ├── quant/     salmon.py tximport.py
│   ├── de/        deseq.py filters.py stats.py
│   ├── function/  ora.py gsea.py disease.py lite.py
│   ├── proteo/    cptac.py pdc.py concordance.py
│   ├── ncrna/     biotypes.py mirna.py targets.py
│   ├── viz/       palette.py static.py interactive.py annotate.py
│   ├── report/    builder.py templates/*.html.j2 assets/
│   ├── io/        naming.py readme_companion.py polars_utils.py
│   └── cli.py     (typer)
├── notebooks/     01_gather … 08_ncrna (student + solutions; 05–08 also run in JupyterLite)
├── lite/          JupyterLite config + bundled data (built by `ctx lite-build`)
├── resources/     references, indices, GMTs (gitignored; fetched by ctx setup)
├── data/          raw/ interim/ processed/ (gitignored)
├── reports/       index.html + committed example reports
├── figures/       pdf light/dark
├── tests/         pytest, chr22-sized fixtures
├── docs/          glossary.md teacher_guide.md parameter_cards/
└── .github/workflows/ci.yml
```

- Runtime profiles (D2): one CLI, three profiles — `laptop` (pixi on Linux / WSL2 / Intel Mac; Docker image for Apple-silicon Macs), `hpc` (Apptainer SIF + SLURM submit script from the template; HPC2N configs), `lite` (browser). Thread and memory budgets come from the profile; commands do not change.
- Polars everywhere for tables (lazy scan, chunked concat, Int32/Float32 downcast after range check, `estimated_size()` printed).
- Google-style docstrings + type hints on every function; `ruff`, `mypy`, `interrogate` in CI.
- Output naming `YYYY-MM-DD_<module>_<description>_v<X>.<ext>`; companion `_README.md` (timestamp, source script, versions, data dictionary) written by `ctx.io.readme_companion`.
- Optional later: nf-core-style Nextflow wrapper for M2–M4; docs map each step to nf-core/rnaseq.

---

## ✅ Phased checklist

### Phase 0 — Decide and pin (this check-in)

- [x] D1 cohort ladder: innate immune cells → other immune cells → large solid-tumour cohorts
- [x] D2 runtime: student laptops first, HPC profile kept
- [x] D4 interactivity: client-side sliders plus JupyterLite
- [x] D11 tier-1 proteomics: Rieckmann 2017 immune-cell proteome, cross-dataset
- [x] D12 Ribo-seq: concept only inside M7
- [x] D13 miRNA: processed TCGA miRNA-seq only, no raw small-RNA-seq in Track A
- [ ] D5 confirm tier-1 dataset (E-MTAB-11190 recommended) after the verification items below
- [ ] Pin tool versions in `pixi.toml` (verify list at the end)

### Phase 1 — Foundations

- [ ] Scaffold repo from `codebase_template`; CI skeleton; `knowledge/lessons.md` in place
- [ ] Feasibility spike: `pydeseq2`, `plotly`, `scipy` under Pyodide via micropip; record outcome; choose the DE strategy for JupyterLite
- [ ] `ctx.viz` palette/static/interactive/annotate with one round-trip test (same data → PDF light/dark + HTML)
- [ ] `ctx.report.builder` with theme toggle, parameter card, quiz and reproduction-block components
- [ ] `ctx.io` naming lint + README companion; runtime profiles `laptop` / `hpc`
- [ ] Chr22 fixture data for tests (tiny FASTQ pair, tiny count matrix)

### Phase 2 — Tier 1 cohort (highest value first)

- [ ] M1 tier-1 acquisition: SDRF metadata → harmonised sample sheet; processed counts if provided, otherwise Salmon counts produced once on Kebnekaise and released with the teaching subset
- [ ] M5 pydeseq2 on tier 1: cell-type axis, cancer-vs-healthy contrasts within cell type, all M5 figures, precomputed knob grid
- [ ] M6 ORA/GSEA/disease association with immune gene sets; `ctx.function.lite` validated against gseapy
- [ ] Reports M1, M5, M6 + hub; sanity checks pass (see verification)

### Phase 3 — Track A on the tier-1 subset

- [ ] Instructor-side: subsample FASTQs on Kebnekaise, build chr22 STAR index, host subset + indices on Zenodo
- [ ] M2 fastp/seqkit/MultiQC parsing and report
- [ ] M3 STAR-chr22 / featureCounts / strandedness; HISAT2 optional path
- [ ] M4 Salmon + pytximport; cross-method Venn; per-sample comparison against the cohort counts of the same samples

### Phase 4 — Tiers 2 and 3, integration modules

- [ ] Tier 3: GDC acquisition, TCGA paired analysis (KIRC), TOIL tumour-vs-GTEx comparison with batch diagnostics
- [ ] Tier 2: TCGA-LAML or DLBC vs healthy comparator; cross-study batch lesson
- [ ] Pan-cancer / pan-cell-type Venn and UpSet views
- [ ] M7 tier 3: CPTAC-3 ccRCC RNA–protein join, concordance analyses, all M7 figures
- [ ] M7 tier 1: neutrophil vs monocyte RNA vs protein copy-number comparison (Rieckmann 2017)
- [ ] M8: biotype-stratified QC/DE on tiers 1 and 3; TCGA-KIRC miRNA-seq target analysis with its protein footprint

### Phase 5 — JupyterLite and pedagogy layer

- [ ] `ctx lite-build`: bundle notebooks 05–08 + data; deploy to GitHub Pages; offline fallback documented
- [ ] Parameter cards, prompts, quizzes, glossary, "common mistakes" for all modules
- [ ] Notebooks (student + solutions) executed headless in CI; JupyterLite smoke test in a headless browser
- [ ] Teacher guide with timing per module

### Phase 6 — Packaging and pilot

- [ ] Docker + Apptainer builds; Apple-silicon and Windows/WSL2 smoke tests; HPC profile test on Kebnekaise
- [ ] Runtime/RAM benchmarks recorded in reports
- [ ] Pilot with 2–3 students; collect friction points into `knowledge/lessons.md`; release v1.0 with Zenodo DOI

---

## 🔍 Verification (definition of done, per module)

- [ ] `ctx <module>` runs end to end on fixture data in CI in < 5 min
- [ ] HTML opens with network disabled; no `<script src="http…">`; size < 15 MB
- [ ] All figures exist in PDF light/dark + HTML; captions carry test, n, p, effect size
- [ ] Every output has `_README.md` with data dictionary; names pass the naming lint
- [ ] Biology sanity, tier 1: PCA separates cell types before disease; cancer-vs-healthy contrasts recover the source paper's reported signatures; p-value histograms show the expected anti-conservative shape. Tier 3: KIRC shows CA9/NDUFA4L2 up and hypoxia/glycolysis enriched
- [ ] JupyterLite notebooks reproduce the CLI results on the bundled data within tolerance (tests compare tables)
- [ ] M7 sanity: per-gene RNA–protein ρ distribution matches the source publication; complex subunits show lower ρ than metabolic enzymes; missingness plotted before any test
- [ ] M8 sanity: biotype fractions match the library type; known KIRC miRNA–target pairs recovered (e.g. hypoxia-induced miR-210); mature miRNAs absent from the mRNA-seq counts, as expected
- [ ] Paired vs unpaired designs produce the expected difference in DE counts (documented)
- [ ] Notebook executes headless; parameter cards match the code defaults (test reads defaults from CLI)
- [ ] Peak RAM and wall time measured (`/usr/bin/time -v`) and within budget

---

## ⚠️ Risks and mitigations

| Risk | Mitigation |
| --- | --- |
| Bioconda lacks osx-arm64 builds for STAR/Salmon/HISAT2 | Docker for Mac users; course JupyterHub/HPC as primary runtime (D2) |
| Salmon decoy-aware index too heavy to build/run on 8 GB laptops | Host prebuilt index on Zenodo or use no-decoy index; explain trade-off |
| GDC data-release or API changes | Pin release in `ctx setup`; cache manifests; tests use fixtures |
| Enrichr/Open Targets availability during class | GMTs cached at setup; Open Targets results cached with timestamp |
| Licence friction (KEGG, COSMIC CGC, DisGeNET API) | Use MSigDB Hallmark/GO/Reactome, NCG, Enrichr GMTs |
| pydeseq2 API drift between releases | Pin version; wrap in `ctx.de.deseq` adapter; test against fixture results |
| Report files large (embedded plotly.js ~3.5 MB each) | Accepted for self-containment; heatmap/table caps |
| Students over-interpret enrichment | Background-trap demo and "association ≠ causation" callouts in M6 |
| Pyodide lacks a needed package (`pydeseq2` deps, `gseapy` Rust backend) | Phase 1 spike decides; pure-NumPy `ctx.function.lite`; precomputed DE grids as fallback |
| Tier-1 processed counts absent or built with a different pipeline | Re-quantify all samples once on Kebnekaise with the course Salmon setup; versions documented |
| Cross-study comparisons (TCGA vs GTEx, TCGA-LAML vs GEO) confounded by batch | Presented as the lesson, not hidden: TOIL-harmonised data, batch diagnostics, limits stated in reports |
| Proteomics access changes (`cptac` package, PDC API) | Processed tables cached in `resources/` at setup; PDC GraphQL as fallback; versions recorded |
| Proteomics pitfalls mislead students (missing values, TMT ratio compression, batch) | Explicit QC section in M7; missingness plotted; tests restricted to proteins quantified in ≥ 70 % of samples |
| miRNA target lists are noisy | Predicted and validated targets shown separately; enrichment tests use matched backgrounds |

---

## ❓ Decisions (resolved, or recommended default in bold)

| ID | Decision | Status |
| --- | --- | --- |
| D1 | Cohort ladder | Resolved: innate immune cells → other immune cells → large solid-tumour cohorts |
| D2 | Student runtime | Resolved: laptops first (pixi / Docker), HPC profile kept |
| D3 | Alignment demo | **STAR on chr22 (mirrors GDC)** + optional HISAT2 full genome / HISAT2 only |
| D4 | Interactivity depth | Resolved: client-side sliders + JupyterLite |
| D5 | Tier-1 dataset | **E-MTAB-11190** / Cassetta 2019 monocytes / other |
| D6 | Hosting of teaching subset and indices | **Zenodo with DOI** / GitHub release assets |
| D7 | Orchestration | **Python CLI first**, Nextflow wrapper later / Nextflow from the start |
| D8 | Licences | **Code MIT, materials CC BY 4.0** |
| D9 | Audience and language | **MSc level, English** / BSc / other |
| D10 | Tier-3 primary solid tumour | **KIRC** / BRCA / LUAD |
| D11 | Tier-1 proteomics source | Resolved: Rieckmann 2017 immune-cell proteome, cross-dataset |
| D12 | Ribo-seq bridge | Resolved: concept only inside M7 |
| D13 | Small-RNA-seq raw reads in Track A | Resolved: processed TCGA miRNA-seq only |
| D14 | Placement of the plan files under the template | Resolved: `TODO.md`, `knowledge/lessons.md`, and `CONVENTIONS.md` imported by `CLAUDE.md`, overriding `AGENTS.md` where they conflict |
| D15 | labcore at runtime | Resolved: no; ctx is self-contained, labcore stays for labdocs |
| D16 | Files per figure | Resolved: `_white`/`_black` PDF and SVG from one render, one HTML with a theme toggle, CSV with `_README.md` |

### Items you could search or verify for me (saves tokens)

- E-MTAB-11190: ENA project accession, FASTQ availability and layout, read length, SDRF fields (cell type, condition, donor, batch), whether processed counts are provided, source publication and its positive-control genes, licence.
- Cassetta et al. 2019 (Cancer Cell): GEO accession and whether FASTQs are on ENA.
- Current TCGA-KIRC / BRCA / LUAD paired tumour–normal counts and TCGA-LAML / DLBC case counts in the current GDC release; GENCODE version behind current STAR-Counts; current UCSC Xena TOIL file names.
- Current release and API of `pydeseq2` (design formula, `vst()`, `lfc_shrink`, `lfc_null`) and of `pytximport` (`counts_from_abundance` options); whether either installs under Pyodide.
- Bioconda osx-arm64 availability for `star`, `salmon`, `hisat2`, `subread`, `fastp`.
- Whether a prebuilt decoy-aware GRCh38 Salmon index is publicly downloadable and its size / RAM at quant time.
- Current MSigDB version string accepted by `gseapy.Msigdb` (including C7 ImmuneSigDB) and current Enrichr disease library names.
- CPTAC-3 ccRCC: number of cases with both RNA-seq (GDC) and proteomics (PDC); current `cptac` package status and data route; TMT batch design.
- Rieckmann et al. 2017 (Nat Immunol): PRIDE accession and whether the per-cell-type copy-number table is downloadable.
- TCGA-KIRC miRNA-seq: number of solid-tissue-normal samples at GDC; file types (miRNA and isoform quantification).
- TargetScan and miRTarBase current versions and download terms; GENCODE Ribo-seq ORF annotation availability.
- `inmoose` (Python ports of limma/edgeR/DESeq2): current status for moderated t-tests.

---

## 📝 Review (filled at completion)

- What was built, what deviated from plan and why
- Measured budgets per module
- Lessons captured in `knowledge/lessons.md`

---

## 🧱 Phase 1 — implementation order (approved 2026-09-26)

Phase 1 depends on no dataset. Nothing reads E-MTAB-11190 or any cohort, and P1.2 has a test
that fails if a cohort identifier enters the code, the tests or the config.

Map from the six Phase 1 items above to the steps below:

| Phase 1 item | Steps |
|---|---|
| 1 | P1.1, P1.2 |
| 2 | P1.5 |
| 3 | P1.7, P1.8 |
| 4 | P1.8, P1.9 |
| 5 | P1.4, P1.6 |
| 6 | P1.3, P1.10 |

Decisions D14–D16 were resolved at this check-in. Engineering decisions are recorded in
`knowledge/decisions.md` as D-001 onwards. CLAUDE.md says to work only on the current phase,
and Phase 0 is still open (D5 and pinning). Phase 1 proceeds anyway because Rytis asked for
it. P1.2 pins every dependency the project adds.

- [x] P1.1 Adopt `codebase_template` in place (D14). Files are moved, never deleted or overwritten
  - c1: `git init`, then commit the three files unchanged.
  - c2, moves only:
    - `mkdir knowledge`
    - `git mv CLAUDE.md CONVENTIONS.md`
    - todo.md to `TODO.md` in two `git mv` steps (via a temporary name)
    - `git mv lessons.md knowledge/lessons.md`
  - c3: `copier copy --trust --defaults --skip-tasks --vcs-ref v0.9.10 --data-file <scratchpad answers> gh:Aramburu-Lab/codebase_template .`
    - Answers: flavour `analysis`; languages python and bash; licence MIT; both HPC2N answers blank.
    - Run it under `/usr/bin/time -l` and keep the log.
    - v0.9.11 arrives later through `copier update`, once the tag is pushed.
  - c4: one scripted pass of targeted edits, each asserting exactly one match:
    - Generated `CLAUDE.md`: append `@CONVENTIONS.md` plus one sentence saying CONVENTIONS.md overrides AGENTS.md where they conflict and each divergence is a D-entry in `knowledge/decisions.md`. Add the same sentence to `AGENTS.md` as one pointer line.
    - Fix every path reference that `grep -n 'lessons\.md\|todo\.md\|CLAUDE\.md'` finds in `CONVENTIONS.md`, and in `TODO.md` outside this implementation-order section, which describes the moves. That includes the old title and the "merge this file" instruction.
    - `.gitignore` (the project owns it after generation):
      - add anchored `/data/` and `/resources/`;
      - add `!tests/fixtures/**`;
      - after that negation, re-add `.DS_Store`, `.env`, `*.pem` and `*.key`.
    - Add `LICENSE` (MIT).
    - In `knowledge/decisions.md`, replace the placeholder D-001 with:
      - D-001: `src/ctx` is the installable product, an exception to "no shared helper module".
      - D-002: the moves and `CONVENTIONS.md`.
      - D-003: no runtime labcore (D15).
      - D-004: laptop root defaults to the working directory; hpc requires an explicit root.
  - c5: run by hand the tasks that `--skip-tasks` skipped:
    - `pixi install` (timed)
    - `pixi run prek install`
    - `pixi run labdocs render`
    - commit `pixi.lock` and every render output.
  - Proven by (output pasted under this item):
    - `git diff --name-status -M c1 c2` shows exactly three `R100` lines.
    - `git diff --name-status c2 c3` shows only `A` lines.
    - The copier log has exit 0 and no `overwrite`.
    - The `/usr/bin/time -l` output for copier and pixi.
    - In a fresh scratchpad clone, `pixi install --locked && pixi run preflight` exits 0. At conformance 1, labdocs lint checks only LD008.
  - Evidence, 2026-09-26, macOS 15 arm64:
    - Commits: fdd58a8 (c1), 45632ed (c2), 9800065 (c3), 9829ac2 (c4), 3b5d329 (c5).
    - c1 to c2 is exactly `R100 CLAUDE.md CONVENTIONS.md`, `R100 todo.md TODO.md` and `R100 lessons.md knowledge/lessons.md`. c2 to c3 is 31 `A` lines and nothing else.
    - copier log: exit 0, zero `overwrite`, and `conflict` + `skip` only for `knowledge/lessons.md` and `TODO.md`. The three files keep their c1 SHA-256.
    - `/usr/bin/time -l`: copier copy took 3.08 s wall at 59,047,936 bytes max RSS; pixi install took 25.41 s wall at 564,822,016 bytes max RSS.
    - `git check-ignore --no-index`:
      - ignored: `data/x.csv`, `resources/manifest.json`, `tests/fixtures/chr22/.DS_Store`, `tests/fixtures/x.key`;
      - not ignored: `tests/fixtures/x.fa`, `x.fq`, `x.gz`, `src/ctx/data/x.py`.
    - Fresh clone: `pixi install --locked` exits 0 and `pixi run preflight` exits 0. pytest collected no tests and still passed; P1.2 closes that exit-5 trap.
    - Deviations:
      - The source is `git@github.com:Aramburu-Lab/codebase_template.git` rather than `gh:`, because SSH is the access path that was verified. `_src_path` records it.
      - The D-entries are appended after the template's entry skeleton instead of replacing it, which leaves Copier's lines untouched.

- [x] P1.2 `ctx` package skeleton, gates, and the full lock (after P1.1)
  - `pyproject.toml`:
    - hatchling build, src layout;
    - distribution `cancer-transcriptomics-teaching`, import package and console script `ctx`;
    - mypy `strict` on `src`, interrogate `fail-under = 100`;
    - pytest marker `network`, deselected by default.
  - `src/ctx/__init__.py` takes its version from `importlib.metadata`. `src/ctx/cli.py` is a typer app with `--version`. Add `src/ctx/py.typed`.
  - `pixi.toml`:
    - Append a `[feature.ctx]` block, every dependency `==`-pinned to the resolved lock: typer, jinja2, seaborn, pydeseq2, pypdf, mypy, interrogate, hatchling, plus the editable `ctx` path dependency.
    - Add `ctx` to the `default` and `laptop` environments only. It stays out of `prod`, because the Dockerfile builds `prod` from pixi.toml and pixi.lock alone. Phase 6 installs the wheel there.
    - pydeseq2 goes into the lock now, so the template's typing_extensions `--locked` trap shows up before anything depends on the environment.
  - `.pre-commit-config.yaml`: a second `- repo: local` block at the end with mypy and interrogate on `^src/`. CI's prek step and preflight's `hooks` task then run them with no further edit.
  - `tests/conftest.py`: fail the session if zero tests are collected, because the template's `test` task maps exit 5 to success.
  - Every deliberate lazy import carries `# noqa: PLC0415` and a reason; the template's ruff selects `PL`.
  - Proven by:
    - `tests/test_gates.py`:
      - `::test_gates_fail_on_violation`: copy `src` to tmp and add one untyped, undocumented function. mypy and interrogate, run with the project config, exit non-zero. Both exit 0 on the clean copy.
      - `::test_zero_collected_fails`: `python -m pytest -k zz_no_match`, run as a subprocess, exits non-zero. Everything is deselected, so the test does not recurse.
      - `::test_no_labcore_import`: an AST scan of every `.py` under `src/ctx` finds no `labcore` import, including inside functions. Files scanned equals the pathlib count, and the count is above 0.
    - `tests/test_packaging.py::test_wheel_runs_outside_repo`:
      - build the wheel (the command is recorded at first run) and unzip it to tmp;
      - a child process asserts the repo's `src` is on `sys.path`, removes it, then imports `ctx`;
      - `ctx.__file__` is under tmp, and `ctx --version` equals the version read from pyproject.toml with tomllib.
    - `tests/test_repo.py`:
      - `::test_layout`:
        - exactly one root entry lowercases to `todo.md`, and it is `TODO.md`;
        - every `@path` in `CLAUDE.md` resolves;
        - a case-sensitive check finds no bare `lessons.md` and no lowercase `todo.md` in `CONVENTIONS.md`, `AGENTS.md` or `CLAUDE.md`;
        - the precedence sentence is in both `CLAUDE.md` and `AGENTS.md`.
      - `::test_gitignore_rules` (`git check-ignore --no-index`):
        - ignored: `data/x.csv`, `resources/manifest.json`, `tests/fixtures/chr22/.DS_Store`, `tests/fixtures/x.key`;
        - not ignored: `tests/fixtures/x.fa`, `x.fq`, `x.gz`, `src/ctx/data/x.py` (the last proves the anchoring).
      - `::test_no_cohort_identifiers`: working-tree files under `src/`, `tests/` and `scripts/`, plus pyproject.toml and pixi.toml, match none of `E-MTAB-\d+|GSE\d+|TCGA-[A-Z]+|PXD\d+|PDC\d+`. The file count is above 0.
    - `pixi run preflight` exits 0.
  - Evidence, 2026-09-26, commit 782817e:
    - The lock resolved pydeseq2 0.5.4 (bioconda), typer 0.27.2, jinja2 3.1.6, seaborn 0.13.2, mypy 2.3.1, interrogate 1.7.0, pypdf 6.19.0, hatchling 1.32.4, threadpoolctl 3.7.0 and typing_extensions 4.16.0. The typing_extensions `--locked` trap did not fire. Pinning `==` left every resolved package unchanged.
    - `pixi install` took 34.09 s wall at 788,873,216 bytes max RSS (`/usr/bin/time -l`). A fresh clone installs `--locked -e default` with exit 0 and passes 16/16.
    - `pixi run preflight` exits 0. Seven hooks pass, including `mypy-src` and `interrogate-src`.
    - Guard check: the same `-k zz_no_such_zz` run exits 5 under `--noconftest` and 1 with the guard.
    - `prod` still cannot install on osx-arm64, exactly as at c5: the template's `tools` feature is linux-64 only.
    - Deviations:
      - Two features instead of one: `ctx`, the runtime in `default` and `laptop`, and `ctx-dev`, the gates in `default` only, because the template defines `laptop` as analysis-only.
      - `test_hooks_run_the_gates_on_src` was added, so a hook whose `files` filter matches nothing fails a test.
      - The wheel is built with `python -m hatchling build -t wheel -d <dir>`.

- [x] P1.3 Synthetic count matrix with planted truth (after P1.2)
  - `tests/synth.py::counts(seed=42)` lives in the tests, not the wheel; P1.5 loads it by path. Dimensions are fixed:
    - 2 cell types x 2 conditions;
    - 6 donors per condition, nested in condition, each donor giving both cell types: 24 samples;
    - 1,000 genes: 50 cell-type markers, 50 shared condition effects, 50 condition-in-one-cell-type effects, 100 low-count, 750 null;
    - negative-binomial counts with per-sample size factors;
    - a `batch` column assigned with no planted effect.
  - Returns Int32 counts, a sample sheet (`sample, donor, cell_type, condition, batch`) and a truth table. It is regenerated each session and nothing is committed.
  - Proven by `tests/test_synth.py`:
    - `::test_design_nested`: each donor has one condition and both cell types; there are 24 unique samples; counts are non-negative and fit in Int32.
    - `::test_seed_determinism`: the same seed gives identical frames and a different seed gives different ones. Stability across numpy versions is not asserted; see P1.5.
    - `::test_planted_truth`:
      - within each cell type, the sign of the normalised log2 ratio matches the truth for every gene with planted |LFC| >= 2;
      - no null gene reaches that threshold;
      - PC1 of log-CPM separates the cell types.
  - Evidence, 2026-09-26: 3/3 pass.
    - Mutation check: planting the condition effect in the controls instead of the cases fails `test_planted_truth` and passes the other two.
    - Deviation: markers are 10 % of genes at |LFC| 4, and nulls 70 %. With 5 % at a smaller effect, the condition axis would carry about as much variance as cell type, and PC1 could split on condition.
    - Added beyond the plan: a per-donor effect (log2 SD 0.25), so the nesting matters, and a `batch` column balanced against condition with no planted effect.
    - `counts()` also takes `n_genes` and `donors_per_condition`, for the P1.5 size grid.

- [x] P1.4 `ctx.io`: naming lint and README companion (after P1.2)
  - `src/ctx/io/naming.py`:
    - `MODULES`: ordered (token, NN) pairs, the single list. Hyphens become underscores, as in `lite_build`.
    - `stem(module, description, version, *, day)`: the only name builder.
    - `figure_files(stem)`: the D16 file set, the only place it is defined.
    - `check(name)`.
    - `lint(directory) -> (scanned, violations)`.
  - `src/ctx/io/readme_companion.py::write_table(frame, directory, stem, *, columns={name: (unit, meaning)}, data_releases, formats=("parquet", "csv"))`
    - Order: validate, then write, then hash, and write the README last.
    - The README records: UTC timestamp; the command line; ctx version and source; `versions()`; data releases; SHA-256 per file; rows; `estimated_size()`; and a dictionary of name, type, unit and meaning.
    - The source is the git SHA only when `ctx.__file__` is `<toplevel>/src/ctx/__init__.py`; otherwise it says "installed wheel".
  - `src/ctx/runtime.py::versions()`: the one definition. The reproduction block reuses it.
  - Provenance lives in a sidecar README, never as a stamp inside the table. CLAUDE.md overrides the global stamp rule here, and every CSV stays safe for readers.
  - Proven by `tests/test_io.py`:
    - `::test_stem_token_order`:
      - `stem("de", "volcano", 1, day=date(2026, 9, 26)) == "2026-09-26_de_volcano_v1"`;
      - for every MODULES token crossed with a set of descriptions and versions, every name in `figure_files()` passes `check()`.
    - `::test_rejects`: `_white_v1`, `_light`, `_dark`, uppercase, a space, a hyphen in the description, `2026-02-30`, a missing or non-integer version, and an unknown module are all rejected.
    - `::test_lint_counts`: a tmp tree with k good and m bad files gives `scanned == k + m` and flags exactly the m bad ones. An empty tree gives `scanned == 0`.
    - `::test_companion_validates_before_writing`: a missing unit, a missing meaning or a bad stem raises and leaves the directory empty.
    - `::test_companion_describes_file_as_written`:
      - planted Int32 and Float32 columns are reported as such for Parquet;
      - for CSV, the reported types are what a CSV reader returns;
      - every README SHA-256 equals hashlib of the file on disk;
      - Parquet and CSV read back equal.
    - `::test_provenance_not_from_cwd`: run from inside an unrelated git repo, and with the unzipped wheel placed inside another git repo, the README carries neither repo's SHA.
  - Evidence, 2026-09-26: `tests/test_io.py` passes 23/23, and preflight passes with 42 tests.
    - Mutation check: running git in the working directory instead of ctx's checkout fails `test_provenance_not_from_cwd`.
    - Added:
      - `test_readme_written_last`, an order test: a CSV write that fails leaves only the Parquet file and no README.
      - `test_stem_raises_on_illegal_parts`.
    - Deviations:
      - The dirty flag is `git status --porcelain -- src`, which counts untracked source files. The first version ignored untracked files and reported a clean SHA for modules that were not committed yet.
      - `versions()` lists the distributions behind every imported module, not a fixed list.
      - `report` is 09 and `lite_build` is 10 in MODULES. todo.md numbers only M0 to M8.
      - `estimated_size()` is recorded in each README rather than printed.

- [ ] P1.5 Pyodide feasibility spike and the JupyterLite DE strategy (after P1.3, P1.4)
  - `scripts/python/spike_pyodide.py` has an `exempt:` codebase-meta block. It holds `decide(record)`, which is committed before the run.
  - The rule, fixed before the run. Choose live pydeseq2 in JupyterLite only if all of these hold:
    - micropip installs it (explicit pins are allowed and recorded);
    - `~ cell_type + condition` fits;
    - the result agrees with a CPython fit of the same pydeseq2 version: every gene whose padj lies outside [0.04, 0.06] falls on the same side of 0.05, and max |dLFC| <= 1e-3;
    - install plus fit takes <= 120 s in the browser, with wheels plus data <= 50 MB at the declared size.
  - Otherwise the fallback is the D4 fallback: precomputed grids plus live scipy tests. If plotly fails, the lite notebooks use matplotlib.
  - Stage A runs in Node with npm `pyodide@314.0.6`, the version jupyterlite-pyodide-kernel 0.8.6 uses (github.com/jupyterlite/pyodide-kernel/releases, retrieved 2026-09-26). `npm install` runs in the scratchpad.
    - Install with `micropip.install(..., keep_going=True)`.
    - Record the npm integrity hash; the numpy, scipy and scikit-learn versions on both sides; wall times; each wheel's filename and bytes; and the first blocking error verbatim.
    - Known risk: anndata 0.13.4 requires typing-extensions >= 4.16, and Pyodide 314.0.6 bundles 4.15.0 (pypi.org/pypi/anndata/0.13.4/json and pyodide.org/en/314.0.6/usage/packages-in-pyodide.html, both retrieved 2026-09-26).
    - Timed on a size grid taken from no cohort: 1,000 x 24, 20,000 x 24 and 20,000 x 200 (genes x samples). The verdicts hold only at those sizes, and Phase 5 reads its bundle's cost off the resulting curve.
  - Stage B runs only if stage A passes: a minimal JupyterLite site in the scratchpad, served on 127.0.0.1. Timing comes from `performance.now()` deltas read in the browser-pane console.
  - Outputs:
    - `deliverables/YYYY-MM-DD_lite_build_pyodide_spike_v1.csv` plus its `_README.md`, written with `write_table(formats=("csv",))` as a long table: stage, component, field, value.
    - The D-005 entry.
    - A D17 row in the decisions table.
    - The D4 bullet under Interactive HTML design, rewritten with the outcome and anchored on the whole bullet.
  - Proven by `tests/test_spike_record.py`, which runs in CI with no network:
    - `::test_decide_boundaries`, on hand-built records:
      - every criterion passing gives live;
      - each single failing criterion gives the fallback;
      - missing stage-B metrics give a provisional verdict.
    - `::test_decision_follows_rule`: `decide()` applied to the recorded metrics returns the recorded decision.
    - `::test_rule_predates_record`: `git merge-base --is-ancestor <commit adding decide()> <commit adding the record>`.
    - `::test_record_complete`:
      - the Pyodide version equals the pinned constant;
      - scipy, plotly and pydeseq2 each have a status plus a version or a verbatim error;
      - the wheel bytes sum to the recorded total;
      - the recorded `counts()` call, with every argument explicit, regenerates the recorded input SHA-256.
      - A numpy bump that changes `counts()` fails this test by design; the spike is then rerun.

- [ ] P1.6 Runtime profiles `laptop` / `hpc` and resource measurement (after P1.4)
  - `src/ctx/runtime.py`:
    - Profile: `--profile`, then `CTX_PROFILE`, then `laptop`.
    - Root: `--root`, then `CTX_ROOT`, then the working directory, on laptop only. `hpc` raises unless a root is given, and the message names both options. The root is printed and recorded in every README.
    - `laptop`: threads = min(4, cpu count); budgets of 6 GB for Track A and 4 GB for Track B, read from the budgets line in CONVENTIONS.md.
    - `hpc`: threads and memory come from the Slurm variables, with no defaults. It raises outside a job.
    - The profile is applied by one decorator at the start of each command, before any lazy heavy import. It is not applied in the group callback, because click runs the callback before a subcommand's `--help`, so `--help` would fail on a login node.
    - It sets `POLARS_MAX_THREADS`, `OMP_NUM_THREADS` and `OPENBLAS_NUM_THREADS`.
    - `measure()`: wall time from `perf_counter`; peak RSS from `getrusage` SELF and CHILDREN (bytes on macOS, KiB on Linux). It is valid once per fresh process. The reproduction block states that CHILDREN reports the largest reaped child.
  - Targeted edits, anchored on whole lines:
    - the budget-measurement sentence in CONVENTIONS.md, and the peak-RAM line of the Verification list;
    - both become: the ctx reproduction block, or for other commands `/usr/bin/time -l` on macOS and `-v` on Linux.
  - Proven by `tests/test_runtime.py`:
    - `::test_precedence`: for both profile and root, the flag beats the environment variable, which beats the default (an order test). `hpc` with no root raises.
    - `::test_hpc_requires_slurm`: with the Slurm variables unset, it raises and names the variable. `SLURM_CPUS_PER_TASK=7` gives 7 threads.
    - `::test_laptop_threads`: with cpu count patched to 2 it gives 2; with 16 it gives 4.
    - `::test_budgets_match_conventions`: the budgets in `runtime.py` equal the numbers parsed from CONVENTIONS.md.
    - `::test_import_is_light`: `import ctx.cli` in a subprocess leaves polars, numpy, matplotlib, plotly and seaborn out of `sys.modules`.
    - `::test_child_peak_attributed`, run in a fresh `python -c` process:
      - a child writes one byte per page of 256 MiB; the children peak is >= 256 MiB and < 1 GiB, and the self peak is < 256 MiB;
      - a control child that allocates nothing leaves the children peak < 256 MiB;
      - it runs on macOS locally and on Linux in CI, which proves the unit handling on both.

- [ ] P1.7 `ctx.viz` style, annotation and the box figure, static PDF and SVG (after P1.3, P1.4)
  - `src/ctx/viz/style.py` is the only style config:
    - Data colours in the CLAUDE.md:48 order, first seven only. `#000000` is the light-theme ink and never a category, so the cap is 7.
    - A role map, such as condition level to colour.
    - Per theme: transparent canvas, grid and ink.
    - Font chain: Helvetica, Arial, Liberation Sans, DejaVu Sans.
    - 14 pt base with role sizes.
    - Line widths: 1 pt for axes and lines, 1.5 pt for box outlines and error bars.
    - Widths of 85 or 170 mm.
    - rc applied in this order: `sns.set_context("paper")`, then the sizes and widths.
    - `pdf.fonttype` 42 and `svg.fonttype` none.
    - A plotly light template on `plotly_white` and a custom dark template.
    - CSS variables.
  - `#F0E442` is 1.32:1 against white, the lowest of the seven (WCAG formula, computed 2026-09-26). style.py marks it not safe for projection, and the setup report prints every data colour's contrast on both backgrounds.
  - `src/ctx/viz/annotate.py::compare_groups` runs, in order:
    - drop NaN, then count n;
    - Mann-Whitney U;
    - rank-biserial r, with one documented sign;
    - BH over the plotted genes (`scipy.stats.false_discovery_control`);
    - one Caption string, which both renderers print.
  - `src/ctx/viz/box.py`:
    - `BoxData` contract: a long frame, units, hover columns and the Caption.
    - `static(data, theme)` builds the matplotlib figure.
    - `interactive(data)` builds the plotly figure. Hover shows gene, group, sample, value with unit, n, adjusted p and r.
  - `src/ctx/viz/save.py::write_figure` writes `figure_files(stem)`: PDF and SVG, `_white` and `_black`, all from one Figure, plus the CSV and README through `write_table`. P1.8 adds the HTML.
  - Proven by `tests/test_viz.py`:
    - `::test_palette_contract`: the data colours equal the CLAUDE.md:48 list without `#000000`, in order; more than 7 groups raises; the CSS variables equal the palette in both themes.
    - `::test_no_style_literals`: no file in `src/ctx` other than `style.py` has a hex colour, `rgb(`, or a numeric `fontsize=`, `linewidth=` or `figsize=`. Files scanned equals the pathlib count, and the count is above 0.
    - `::test_context_then_sizes`: rcParams sizes and widths equal style.py's values, not seaborn's paper-scaled ones (an order test).
    - `::test_font_fallback`: with the chain (`NoSuchFontXYZ`, DejaVu Sans), the PDF embeds DejaVu Sans and the SVG names it.
    - `::test_bh_input_order`: an unsorted p vector gives the hand-worked BH result, in input order and capped at 1.
    - `::test_bh_after_subset`: BH runs over the plotted subset only.
    - `::test_nan_dropped_before_n`.
    - `::test_rank_biserial_truth`: complete separation gives +1 or -1, and identical groups give 0.
    - `::test_planted_recovered`:
      - first assert that the smallest attainable p for 6 vs 6 (2/924) times the number of plotted genes is below 0.05;
      - then assert that the planted genes reach BH < 0.05 with the planted sign, and none of 5 null genes does.
    - `::test_static_files_measured`, using pypdf and xml.etree.
      - Both PDFs:
        - no Type 3 font; every BaseFont is in the chain;
        - Tf sizes are a subset of the style sizes, and `w` operands a subset of the style widths;
        - no image XObject, and no filled rectangle covering the MediaBox;
        - MediaBox width is 85 or 170 mm, within 0.5 pt;
        - the ink colour equals the theme's ink;
        - the caption tokens (test, n, p, r) appear in the text after whitespace normalisation;
        - the y-axis label names its unit.
      - Both SVGs: every font-family is in the chain, and there is no `<image>`.
    - Both PDFs are opened once by eye, and that is noted here.

- [ ] P1.8 `ctx.report.builder`, the figure HTML and the round trip (after P1.6, P1.7)
  - `src/ctx/report/templates/base.html.j2`, with autoescape on:
    - a theme script in `<head>`: the stored `ctx-theme` if present, else `prefers-color-scheme`;
    - a toggle that calls `Plotly.relayout` with the embedded light and dark templates;
    - a relayout to the resolved theme once the plots are drawn;
    - plotly.js inlined once;
    - CSS variables from style.py.
  - Two page kinds are built on it:
    - `figure.html.j2`: the figure, its caption and a data link.
    - `module.html.j2`: the eight CLAUDE.md:56 sections, with macros for the parameter card, think-about-it, `<details>` self-check, common mistakes and the reproduction block.
  - `src/ctx/report/builder.py`: `render(page) -> str`, and `write(page, path)`, which raises before writing when the page is over 15 MB. Templates load through `importlib.resources`.
  - `write_figure` adds `<stem>.html` through the figure page, so there is one toggle implementation. The data table is embedded as its own JSON script, taken from the frame that wrote the CSV.
  - Proven by `tests/test_report.py`:
    - `::test_round_trip`:
      - one `BoxData` gives exactly `figure_files(stem)`;
      - `lint(outdir)` gives scanned equal to the file count and no violations;
      - the caption tokens are identical in the PDF text, the SVG text, the plotly JSON and the README;
      - the embedded table rows equal the CSV rows;
      - the hover covers every hover column in the contract.
    - `::test_html_offline`:
      - the inlined block equals `plotly.offline.get_plotlyjs()` and occurs exactly once;
      - outside that block there is no `<script src`, no `<link href` and no CSS `url(` pointing to http or `//`;
      - the file is <= 15 MB;
      - the embedded templates equal style.py's.
    - `::test_size_guard`: a page over 15 MB raises and writes nothing.
    - `::test_skeleton_order`: the module page's section ids appear in the CLAUDE.md:56 order, and a missing section raises (an order test).
    - `::test_escaping`: `</script><b>x` in the figure data and in the quiz text adds no element.
    - `::test_quiz_hidden`: at least one `<details>` answer exists, and every answer sits inside `<details><summary>`.
    - `test_packaging.py` is extended to render a page from the unzipped wheel.
    - Manual, once: serve the page on 127.0.0.1 in the browser pane under light and dark colour-scheme emulation. The request log must show one origin, and screenshots before and after the toggle must show the background and ink changing. The results are pasted here.

- [ ] P1.9 `ctx setup` environment report as the end-to-end consumer (after P1.6, P1.8)
  - This is the M0 subset from `TODO.md:102`; it downloads nothing.
  - Every CLI parameter gets a parameter card (meaning and effect), including the global `--profile` and `--root`, which are rendered once as a shared block.
  - Outputs under `<root>/outputs/00_setup/`:
    - `YYYY-MM-DD_setup_environment_v1.html`, a module page showing:
      - versions and hardware;
      - which font in the chain resolved;
      - colour contrasts;
      - a reproduction block with the exact command, wall time, peak RSS for self and children, the profile and its budget, and versions.
    - `YYYY-MM-DD_setup_versions_v1` as Parquet and CSV, with its README.
  - Proven by `tests/test_cli.py`:
    - `::test_parameter_cards_match_cli`:
      - the card set is non-empty and equals every click parameter except `--help` and `--version`;
      - the defaults read back from the written HTML equal the click defaults;
      - a throwaway typer app with an option that has no card makes the checker fail.
    - `::test_subcommands_in_modules`: every registered subcommand maps to a MODULES token.
    - `::test_threads_set_before_import`: a subprocess runs with `CTX_PROFILE=hpc`, the Slurm variables set, and a thread count v chosen so that v != `os.cpu_count()`. The polars pool is v, and threadpoolctl reports v for BLAS (threadpoolctl is UNVERIFIED until the lock shows it).
    - `::test_help_without_slurm`: `CTX_PROFILE=hpc ctx setup --help` exits 0 with no Slurm variables set.
    - `::test_code_is_not_data`: with `CTX_ROOT=tmp`, run from another directory, every output lands under tmp and nothing is written inside the package.
    - `::test_repro_block_shows_measured`: with `measure()` stubbed to 12.345 s and 678 MiB, the block shows exactly those values.
    - `::test_setup_end_to_end`:
      - exit 0;
      - lint scans >= 3 files with no violations;
      - every data file has a README;
      - the P1.8 offline checks pass;
      - wall time and peak RSS are above 0 and within the laptop budget.
    - `test_packaging.py` is extended: `ctx setup` runs from the unzipped wheel.

- [ ] P1.10 Chr22 slice and a simulated FASTQ pair with known truth (after P1.2; can run in parallel with P1.3–P1.9)
  - `scripts/python/make_chr22_fixture.py` has an `exempt:` block. It runs instructor-side and needs the network.
    - It pins a GENCODE release, looked up and cited when the script is written. Checksums are the publisher's if GENCODE provides a checksum file; otherwise they are labelled as first-download values.
    - It takes a small chr22 window holding multi-exon protein-coding genes on both strands.
    - It writes `slice.fa`, `slice.gtf` (rebased onto the contig `chr22_slice`) and `transcripts.fa`.
    - `reads(transcripts, n_pairs, seed=42)` lives in this script:
      - it uses only the standard-library `random`, whose `random()` Python documents as reproducible for a seed across versions;
      - it simulates paired-end, reverse-stranded reads with no errors;
      - it encodes each read's origin in its name;
      - it writes `truth.tsv`;
      - it refuses a directory that is not empty.
  - `tests/fixtures/chr22/` holds those files plus `reads_R1.fq`, `reads_R2.fq`, `MANIFEST.sha256` and `README.md`. The README states the release, URLs, hashes, window, offset and seed, and what the fixture can and cannot test. Files are uncompressed, <= 2 MB in total.
  - Run `pixi run docs` again for the new script.
  - Proven by `tests/test_fixtures.py`:
    - `::test_manifest`: the files on disk equal the manifest (count above 0); the hashes match; the total is <= 2 MB; the README hashes equal the manifest.
    - `::test_tracked_not_ignored`: every fixture is in `git ls-files` and not ignored (`--no-index`).
    - `::test_gtf_fasta_consistent`:
      - every exon lies inside the slice;
      - `transcripts.fa` equals the strand-aware splice of `slice.fa` at the rebased 1-based closed intervals;
      - an order test: rebase first, then extract.
    - `::test_reads_known_truth`:
      - pairs parse, names match, and lengths agree;
      - each read equals its transcript substring at the encoded offset, in the orientation the protocol demands;
      - counts equal `truth.tsv` and sum to n pairs.
    - `::test_regenerate_byte_identical`: `reads()`, loaded by path and run from a tmp copy of `transcripts.fa` into an empty tmp directory, gives hashes equal to the manifest.
    - `::test_refuses_nonempty_dir`.
    - `::test_slice_rebuilds` (marker `network`): rebuilding from the pinned release reproduces the manifest. It is run once and the output pasted here.

- [ ] P1.11 Close Phase 1
  - Divergence inventory: in a scratchpad clone, run `copier recopy --trust --defaults --skip-tasks`. `git diff --name-only` must list exactly `pixi.toml`, `.pre-commit-config.yaml`, `AGENTS.md`, `CLAUDE.md` and `knowledge/decisions.md`. Any other file is an unplanned divergence.
  - Real update, in a scratchpad clone of the template:
    - Tag a synthetic release that bumps the labcore and ruff pins, edits the last lines of `CLAUDE.md.jinja` and `AGENTS.md.jinja`, and touches the pre-commit local block.
    - Point a scratch project clone's `_src_path` at that clone and run `copier update --trust --skip-tasks`.
    - It must exit 0 with no `*.rej` files, no conflict markers, and every project line intact.
    - The real template repo is never touched.
  - In a fresh clone, `pixi run preflight` exits 0, and `/usr/bin/time -l pixi run test` takes < 5 min.
  - Review section:
    - what was built, and what deviated and why;
    - a table of measured wall time and peak RSS for copier copy, pixi install, pytest, the spike stages and `ctx setup`, with machine and OS;
    - the repository tree in this file updated by targeted edits;
    - `knowledge/lessons.md` updated;
    - the D-entries listed.
  - The CI run URL is recorded only after Rytis creates the remote and pushes.

### Conflicts with codebase_template

| Topic | Template | This project | Resolution |
|---|---|---|---|
| `CLAUDE.md` | Owned by Copier, not skip-if-exists (copier.yml:120-127) | Holds every rule; "merge, do not overwrite" | D14: moved to `CONVENTIONS.md`, imported, overrides AGENTS.md |
| Plan files | `TODO.md` and `knowledge/lessons.md` are skip-if-exists; `copier update` recreates them if missing | Root `todo.md` and `lessons.md` on a case-insensitive disk | D14: moved |
| Code layout | `scripts/<lang>/`; "no shared helper module" (AGENTS.md.jinja:46-52); labdocs walks only scripts, bin and containers (labcore walk.py:30) | Typer package in `src/ctx` | D-001. Conformance stays at 1; level 5 ("labcore used for viz/IO") is unreachable under D15 |
| Helper channel | Reusable helpers go in labcore | ctx is self-contained | D15; AST test; labcore kept for labdocs |
| `knowledge/api_index.md` | CLAUDE.md.jinja says to read it before writing a helper | It never indexes `src/` | Read `src/ctx` directly (stated in the CONVENTIONS pointer) |
| Output names | `out_<script>_<desc>` in `outputs/`; dated names in `deliverables/`; trailing `_white`/`_black` (AGENTS.md.jinja:62-70) | Dated names everywhere; "light/dark" | Dated names with the theme token after `_vN`; run outputs in `outputs/NN_<module>/`, which labdocs does not check |
| Interactive HTML | Opt-in (figures.toml.jinja:23-26); labcore writes one CDN-linked file per theme (save.py:147) | One HTML with a toggle, plotly embedded | ctx builder, plotly inlined once |
| Tables | Interchange format is TSV (formats.toml.jinja:7); labcore refuses CSV unless formats.toml lists it or `override=True` | Parquet + CSV | ctx writer; `settings/*.toml` untouched and not read by ctx |
| Style | `settings/figures.toml`, which labcore does not read; paper scale 8 pt; labcore's Okabe-Ito swaps positions 4 and 7 | Helvetica 14 pt, CLAUDE.md:48 order | `style.py` is the one config and ships in the wheel |
| Gates | `lint` covers `scripts` only; exit 5 passes; no mypy or interrogate; `PL` rejects lazy imports; the length hook matches only `^(scripts\|bin)/` | mypy and interrogate on `src`; lazy imports needed | prek local hooks, conftest guard, `# noqa: PLC0415` with a reason; `src/` sits outside the length ratchet |
| Pins | Ranges in pixi.toml, exact versions in pixi.lock | "Pin every version in pixi.toml" | Project dependencies `==` in `[feature.ctx]`; the template's ranges pinned by the lock (recorded exception) |
| `.gitignore` | Ignores `*.fa *.fq *.gz *.parquet`; no `data/` or `resources/` rule | Fixtures committed; `data/` and `resources/` never | Anchored rules and a fixture re-include, with secrets and `.DS_Store` re-ignored. Committed deliverables are CSV only |
| Container | Dockerfile builds `prod` from pixi.toml and pixi.lock alone | ctx installed editable | `ctx` feature kept out of `prod`; Phase 6 installs the wheel |
| Peak RAM | POSIX tooling | GNU `time -v`, which is absent on macOS | In-process `getrusage` (P1.6); `-l` on macOS |
| Tree | `deliverables/`, `containers/`, `outputs/NN_<step>/` | `reports/`, `figures/`, root `Dockerfile` | `deliverables/reports`, `deliverables/figures`, `containers/`; tree updated in P1.11 |
| Module layout | none | `viz/palette.py static.py interactive.py` | One module per figure type with `static()` and `interactive()`, plus `style.py` and `save.py`; tree updated in P1.11 |
| Root | No guessed default (paths.toml.jinja:3-7) | none | `hpc` requires a root; `laptop` defaults to the working directory (D-004) |
| Platforms | linux-64 (glibc 2.35) and osx-arm64 | Linux, WSL2, Intel Mac | Unchanged in Phase 1; osx-64 attempted in Phase 6 |
| "laptop" | A pixi environment | A runtime profile | Both kept; ctx reads only `--profile` and `CTX_PROFILE` |
| Licence | One answer, no LICENSE file | Code MIT, materials CC BY 4.0 | MIT and `LICENSE` now; the CC BY file with the first materials (Phase 5) |
| Report sections | none | TODO.md:98 lists 7, CLAUDE.md:56 lists 8 | CLAUDE.md:56's eight are the contract |
| CI comment | "while the repo is private" (ci.yml.jinja:33) | none | labcore is public; the comment is harmless and reported upstream, not edited |

Where the project overrides the global rules:
- Helvetica, not Arial.
- A 14 pt base, not 12/10/9/8.
- The project's palette order, which keeps yellow in position 4.
- Sidecar README provenance instead of a stamp in the table.

Where the project is silent, the global rules apply:
- Black is ink, with at most 7 categories.
- Widths of 85 or 170 mm.
- Line widths of 1 and 1.5 pt.
- PDF and SVG from one render (D16).
- Every axis names its unit, and every figure ships its numbers.

If 14 pt proves unworkable at 85 mm, that gets raised, not changed silently.

Not in Phase 1:
- anything tied to a dataset;
- resource fetching and `resources/manifest.json`;
- `polars_utils`;
- figure types other than the box plot, and panel letters;
- JupyterLite build, `ctx lite-build`, `ctx.function.lite`, and headless-browser tests;
- Slurm submission and containers;
- osx-64, WSL2 and Intel-Mac smoke tests;
- aligner or Salmon recovery of the chr22 truth, and bad-QC reads;
- `CITATION.cff`;
- upstream fixes to labcore or the template, which are filed as issues;
- raising conformance above 1;
- the FINAL package.
