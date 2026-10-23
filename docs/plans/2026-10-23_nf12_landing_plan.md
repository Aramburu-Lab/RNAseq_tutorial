# NF.12: landing nf-core/rnaseq's results in the tutorial pages

Prepared 2026-10-23 while NF.2's six runs go, so each run's results become pages with a copy and
one command. The copy list is NF.1's `ctx.gather.rnaseq.COPIED_FILES`; the step that needs it is
marked, and is filled once NF.1 merges.

## Ready now

- The core root `~/ctx_core`: every resource of the six datasets, set up under the core edition
  on 2026-10-23 with the network blocked (proxies at 127.0.0.1:9), copied by SHA-256 from
  `~/ctx_tier1/resources`, `~/ctx_breast/resources` and `~/ctx_lung/resources`; nothing
  downloaded. GTEx V11's files are a link, `~/ctx_core/data/raw/gtex_v11` to
  `~/ctx_tier1/data/raw/gtex_v11`, as the review build did. Every `ctx setup` exited 0.
- `ctx run --check` on Non-cancerous blood's local inputs, under the new names and the core
  edition, found only resources not yet set up (before this root existed): the renaming leaves
  `ctx run` working.
- M9 is not in `ctx run` by design (`pipeline.py:31`: one command fetches nothing); it is its own
  `ctx genecard` step below.

## When a run ends (Kebnekaise)

1. Its head job's log ends "Pipeline completed successfully": paste the last lines.
2. The first finished run saves its indexes once (`save_index.sh`, NF.10) into
   `$Project_d/Rytis/Data/Reference/Human/Gencode/v50_primary_assembly/nf-core-rnaseq-3.27.0`.
3. When all six have ended, the run folders take the new names (NF.5): `mv` each, with a
   symlink at the old path, since nf-core's and fetchngs' outputs hold absolute paths.

## The copy to the Mac (Rytis runs it on the Mac; one block per dataset)

Two folders per dataset under `~/ctx_core/data/raw/`:
- `<Name>_meta`: the fetchngs outdir without FASTQs (`samplesheet/`, `metadata/`,
  `pipeline_info/`), what `ctx gather` reads.
- `<Name>_rnaseq`: NF.1's copy list (TO FILL at NF.1's merge), never a BAM; its byte count is
  measured on the first dataset.

## The pages (on the Mac, per dataset)

    export CTX_EDITION=core
    ctx --root ~/ctx_core run --dataset <Name> --dataset-file conf/datasets/<Name>.toml \
        --fetchngs-dir ~/ctx_core/data/raw/<Name>_meta --rnaseq-dir ~/ctx_core/data/raw/<Name>_rnaseq --check
    ctx --root ~/ctx_core run ...the same, without --check
    ctx --root ~/ctx_core genecard --dataset <Name> --dataset-file conf/datasets/<Name>.toml

M9 calls HPA and Ensembl for genes the cache does not hold; it is asked for before it runs.
Then `ctx report` per dataset, the lite review set from `~/ctx_core`, and the kit
(`build_kit.py --edition core`). Every step's wall time and peak RSS is recorded
(`/usr/bin/time -l`).
