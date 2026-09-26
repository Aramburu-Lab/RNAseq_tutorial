# Decisions — Cancer Transcriptomics Teaching Suite

Analysis decisions that a reader could reasonably have made differently, with the
evidence that settled them. Not a changelog: record the ones where the outcome
would change if the decision changed.

Newest first. One entry per decision, and never delete one — a superseded
decision is marked superseded, with a pointer to the entry that replaced it.

Template for a new entry:

## D-001 — <one-line statement of the decision>

**Date:** YYYY-MM-DD · **Status:** decided | superseded by D-00N

**Context.** What made this a question.

**Options.** What was actually considered, including the one not taken.

**Decision.** What was chosen.

**Evidence.** The number, plot or run that settled it. A decision with no
evidence line is a preference — say so rather than dressing it up.

**Consequence.** What now depends on this, and what would have to be redone if
it were reversed.

---

## D-006 — Raw reads enter only through nf-core/fetchngs 1.12.0; `ctx gather` is the only ctx reader of its layout (TODO.md D18)

**Date:** 2026-09-26 · **Status:** decided

**Context.** `knowledge/recommendation1` (Rytis, 2026-09-26; moved out of the repo on 2026-09-27 to `/Users/rytis/Scripts/Pipeline/OLD/RNAseq_tutorial_recommendation1_2026-09-26/`) replaced custom ENA/`fastq-dl` download logic with nf-core/fetchngs 1.12.0 and split M1 into find, fetch and make usable. Its layout description followed the fetchngs docs; parts of it do not match the 1.12.0 code.

**Options.** Custom download code in ctx; fetchngs 1.12.0, which the lab runs on Kebnekaise; fetchngs 1.13.0 (released 2026-09-15, needs Nextflow >= 25.10.4, renames the versions file, changes the samplesheet row order, untested in the lab).

**Decision.** fetchngs 1.12.0 with `--nf_core_pipeline rnaseq`, one fresh outdir per run under `data/raw/<slug>/`, read-only to ctx; `ctx gather` (package `ctx/gather/`) is the only ctx module that reads it; no launcher and no ENA/SRA access code in ctx. The reader stays in ctx under D-003.

**Evidence.** From the 1.12.0 code (raw.githubusercontent.com/nf-core/fetchngs/1.12.0, retrieved 2026-09-26):
- File names `<EXP>_<RUN>[_1|_2].fastq.gz` (`bin/sra_runinfo_to_ftp.py:143-145`); `.md5` files only from FTP and Aspera, written from ENA's md5 and checked with `md5sum -c` (`modules/local/sra_fastq_ftp/main.nf:29,44,52`); the sra-tools fallback is automatic per run and writes no md5 (`workflows/sra/main.nf:79-81`).
- `*.runinfo.tsv` is not published (`modules/local/sra_ids_to_runinfo/nextflow.config:3-6`); `metadata/<INPUT_ID>.runinfo_ftp.tsv` is, one per input line.
- `--nf_core_pipeline rnaseq` only inserts `strandedness` = `auto` (`modules/local/sra_to_samplesheet/main.nf:33-49`); sheet paths are `<outdir as given>/fastq/<name>` (same file, `:24-53`).
- The input regex rejects `E-MTAB-*` (`subworkflows/local/utils_nfcore_fetchngs_pipeline/main.nf:144-165`); an SRR pulls every run of its experiment.
- Minimum Nextflow `!>=23.04.0` (`nextflow.config:223`); container images are amd64-only, so Apple silicon needs `-profile docker,arm` (`nextflow.config:110-112`).
- On Kebnekaise the lab runs 1.12.0 with `-profile singularity` on Nextflow 26.04.4, which needs `NXF_SYNTAX_PARSER=v1` (`General_tools/cluster_utils/SHARED_run_fetchngs.sbatch:56-77`; `cluster_utils/README.md`).

Observed while building the Phase 1b fixtures on this Mac (2026-09-26, `-profile docker,arm`, Nextflow 26.04.6):
- `metadata/*.runinfo_ftp.tsv` rows come in no fixed order: a rebuild swapped two runs of GSM4907283. `samplesheet.csv` came back byte-identical. ctx reads metadata as a set of runs; the rebuild test compares those files by rows.
- The sra-tools fasterq-dump step asks for 36 GB by default and fails on a 24 GB laptop; the fixture script caps the run at `--max_cpus 4 --max_memory 8.GB`.
- In one rebuild, `prefetch` inside the amd64-emulated sra-tools container failed to open the run from NCBI's S3 bucket and then sat idle, while the host reached the same URL in 0.4 s. The FTP runs were unaffected. Class runs should prefer runs with an ENA FASTQ link.

**Consequence.** An upgrade to 1.13.0 is one deliberate change: `FETCHNGS_VERSION`, a new fixture, and the renamed versions file. `ctx gather` refuses any other release until then.

## D-005 — JupyterLite DE: precomputed grids plus live scipy tests, not live pydeseq2 (TODO.md D17)

**Date:** 2026-09-26 · **Status:** decided

**Context.** D4 left open whether the M5 notebook re-runs pydeseq2 in the browser. P1.5 fixed the
rule as `decide()` in `scripts/python/spike_pyodide.py`, committed in 247d3f5 before anything ran.
The record is `deliverables/2026-09-26_lite_build_pyodide_spike_v1.csv`.

**Options.** Live pydeseq2 0.5.4 in JupyterLite; the D4 fallback, precomputed grids plus live
scipy tests; live pydeseq2 under an in-process joblib backend, which the rule does not cover.

**Decision.** The fallback at all three grid sizes, 1,000 x 24, 20,000 x 24 and 20,000 x 200, with
plotly for figures. Stage B, the browser timing, did not run, because no size passed stage A.

The rule's "the `~ cell_type + condition` fit" was read as pydeseq2's default call. That reading
was fixed in the harness commit 4a54ff0, after a scratch dry run had shown the failure, not in the
rule commit 247d3f5. Under pydeseq2's public backend argument (the shim below) the 1,000 x 24 fit
meets every stage-A criterion, so the verdict at that size rests on the reading; how it was
settled is TODO.md D17. (The "D-005a" in AGENTS.md is the template's numbering, not this entry.)

**Evidence.** Stage A ran npm `pyodide@314.0.6`, Python 3.14.2, under Node v24.15.0 on macOS
15.7.9 arm64. npm integrity:
`sha512-BKDTJyIqFxC4BExLqeRS3f5xvXZIjOt8C3zGLN/Cc7tFxSwvKVhVkchQQ2AGLOtR4YrVOIFxbV8poyDOOmWwxQ==`.
- micropip installed pydeseq2 0.5.4 as 47 files of 41,800,855 bytes. anndata 0.13.4 arrived with
  typing-extensions 4.16.0 from PyPI, because the bundled 4.15.0 was never loaded, so the known
  risk did not fire in Node.
- The fit failed at every size. First blocking error, verbatim:
  `AssertionError: ThreadingBackend does not accept setting the inner_max_num_threads argument.`
  pydeseq2 0.5.4, the newest release on PyPI (https://pypi.org/pypi/pydeseq2/json, uploaded
  2026-01-23, retrieved 2026-09-26), wraps every inference call in
  `parallel_backend(self._backend, inner_max_num_threads=1)` (`default_inference.py:96`), where
  `backend` is a public `DefaultInference` argument that defaults to `"loky"`. Pyodide has no
  loky, joblib 1.5.3 falls back to threading, and threading rejects that argument. No version pin
  reaches it; the backend argument does.
- Loaded versions, Pyodide then CPython: numpy 2.4.6 and 2.5.3, scipy 1.18.0 and 1.18.1,
  scikit-learn 1.8.0 and 1.9.1, pandas 3.0.2 and 3.0.6, pydeseq2 0.5.4 in both.
- At 20,000 x 200, wheels plus data are 56,385,918 bytes, over the 50,000,000 limit, so that size
  fails even with a working fit.
- Fallback arm: in a Pyodide where `import pydeseq2` raises `ModuleNotFoundError`, scipy 1.18.0,
  pandas 3.0.2 and plotly 7.1.0 installed as 10 files of 32,226,431 bytes. Mann-Whitney U with BH
  matched CPython at every size: max |dp| 2.2e-16, max |dpadj| 5.6e-16, and no gene crossed 0.05.
- Recorded beside the rule, which never reads it: with a sequential joblib backend registered as
  `in_process`, pydeseq2 fits in Pyodide in 6.6, 38.5 and 44.8 s of Node time. It agrees at
  1,000 x 24 (0 disagreements, max |dLFC| 2.2e-5) and at 20,000 x 200 (0, 2.8e-4) but not at
  20,000 x 24: 389 disagreements. 388 are low-count genes, baseMean 0.27 to 1.34, whose padj is
  NaN in Pyodide only (CPython has none NaN); one of them, g05499, has CPython padj 0.0307. The
  389th is well expressed: g15471, baseMean 947.8, p 0.01027 against 0.003241, padj 0.0876
  against 0.0318. The largest p difference is 0.0381 (g01437, baseMean 1208.6). So the fits
  differ beyond independent filtering. Recomputed on 2026-09-26 from the run's raw results files.
- The record's `wheels_downloaded_bytes` holds the bytes of unique package files by name, not
  the bytes downloaded: formulaic, interface_meta and python_dotenv were each fetched twice
  (41,957,904 bytes of wheels fetched against 41,800,855 recorded), and 785,537 bytes of index
  pages and `.metadata` files are excluded. The Pyodide runtime is excluded as well:
  `pyodide.asm.wasm` 9,598,218 and `python_stdlib.zip` 2,545,564 bytes (npm package 314.0.6).
- The fallback's own budget at 20,000 x 200: 32,226,431 bytes of wheels plus 14,585,063 bytes of
  input leave 3,188,506 bytes of the 50,000,000 for the precomputed grids, before the runtime.
- Known gap in the pre-registered rule, left unedited: when both pydeseq2 version rows are
  missing, the version criterion compares None with None and passes. The record carries both
  versions, so this verdict is unaffected; a rerun must fix the rule first, in its own commit.
- Sources, retrieved 2026-09-26: https://registry.npmjs.org/pyodide (314.0.6 integrity; npm
  `latest` is 314.0.7, published 2026-09-14);
  https://registry.npmjs.org/@jupyterlite/pyodide-kernel/-/pyodide-kernel-0.8.6.tgz (devDependency
  `pyodide` 314.0.6); https://pypi.org/pypi/anndata/0.13.4/json (`typing-extensions>=4.16`).

**Consequence.** Phase 5 bundles precomputed DE grids and runs the scipy tests live. Reversing this
needs a new decision: the backend shim changes the notebook's code path, 20,000 x 24 still
disagrees under it, its wheels alone take 41.8 MB of the 50 MB budget, and no browser timing
exists. A numpy change that alters `counts()` fails `tests/test_spike_record.py`; the spike is
then rerun. The proof that the rule predates the record needs full history, so it runs locally
and skips in CI's depth-1 checkout.

## D-004 — Output root: the laptop profile defaults to the working directory; hpc needs an explicit root

**Date:** 2026-09-26 · **Status:** decided

**Context.** `settings/paths.toml` refuses guessed roots, because on the cluster a fallback
writes where nobody looks. Students run `ctx` from a course folder on a laptop.

**Options.** Always require `--root` or `CTX_ROOT`; default to the working directory everywhere;
default to it on the laptop profile only.

**Decision.** `--root`, then `CTX_ROOT`, then the working directory on `laptop`. `hpc` raises
without an explicit root. The resolved root is printed and written into every `_README.md`.

**Evidence.** None measured; this is a preference. The cluster risk is the one the template's
`settings/paths.toml` comment records.

**Consequence.** P1.6 tests the precedence. Reversing it changes where every student output lands.

## D-003 — ctx does not import labcore at runtime (TODO.md D15)

**Date:** 2026-09-26 · **Status:** decided

**Context.** `AGENTS.md` routes every reusable helper through the pinned `labcore` package.

**Options.** Depend on labcore after upstream changes (inline plotly, a theme toggle, a unit
column), or keep ctx self-contained with labcore in the environment for `labdocs` only.

**Decision.** Self-contained. Chosen by Rytis at the 2026-09-26 check-in.

**Evidence.** In labcore v0.6.6: `viz/save.py:147` writes plotly HTML with
`include_plotlyjs="cdn"`, one file per theme; the `write_table` README dictionary is
`Column | Dtype | Description`, with no unit (`io.py:239-273`); a CSV write is refused unless
`settings/formats.toml` lists it or `override=True` is passed; the Okabe-Ito order at
`viz/palette.py:60-61` puts `#CC79A7` fourth and `#F0E442` seventh, the reverse of
`CONVENTIONS.md`.

**Consequence.** Conformance level 5 ("labcore used for viz/IO") is unreachable. P1.2 adds an
AST test that fails on any `labcore` import under `src/ctx`.

## D-002 — Plan files moved to the template's locations; project rules live in CONVENTIONS.md (TODO.md D14)

**Date:** 2026-09-26 · **Status:** decided

**Context.** Copier owns `CLAUDE.md` and does not skip it when it exists. `TODO.md` and
`knowledge/lessons.md` are skip-if-exists, and on this case-insensitive disk `TODO.md` already
resolved to the root `todo.md`.

**Options.** Move and import (chosen); merge the project rules into `AGENTS.md`; keep the files at
the root and hand-merge `CLAUDE.md`.

**Decision.** `todo.md` became `TODO.md`, `lessons.md` became `knowledge/lessons.md`, and
`CLAUDE.md` became `CONVENTIONS.md`, which the generated `CLAUDE.md` imports. Where
`CONVENTIONS.md` conflicts with `AGENTS.md`, `CONVENTIONS.md` wins. `/data/`, `/resources/` and a
`tests/fixtures/**` re-include were added to `.gitignore`.

**Evidence.** Commit 45632ed is exactly three `R100` renames. The `copier copy` log (v0.9.10) has
0 overwrites and `conflict`/`skip` for `TODO.md` and `knowledge/lessons.md` only; all three files
kept their SHA-256.

**Consequence.** `copier update` merges against two short project additions in `CLAUDE.md` and
`AGENTS.md` instead of the whole rule set. They sit at the top of each file, after the import and
after the intro paragraph, not at the end: a real update to a synthetic template release that
appended to both files left conflict markers in each while the lines were appended, and merged
cleanly once they were moved (P1.11, 2026-09-26). A template change to those top lines would still
conflict; that is the smallest surface available.

## D-001 — `src/ctx` is an installable package, an exception to "no shared helper module"

**Date:** 2026-09-26 · **Status:** decided

**Context.** `AGENTS.md` puts every runnable step in `scripts/<lang>/` and forbids a shared
helper module in the repo; `labdocs` walks only `scripts/`, `bin/` and `containers/`.

**Options.** A typer CLI package in `src/ctx`, or scripts only.

**Decision.** The teaching product is the `ctx` package in `src/` (`CONVENTIONS.md`, Stack).
One-off steps such as fixture builders and spikes stay in `scripts/python/` with codebase-meta
blocks.

**Evidence.** None measured; this is a requirement: students install and run one CLI.

**Consequence.** `labdocs` does not document `src/ctx`, `knowledge/api_index.md` does not index
it, and the file-length ratchet does not cover it. Conformance stays at 1 until raised on purpose.
