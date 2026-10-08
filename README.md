# Cancer Transcriptomics Teaching Suite

<!-- Generated once from the lab template, then yours. Template updates never
     touch this file again, so say what this project is for in the first
     paragraph and keep it true. -->

One paragraph: what question this project answers, on what data, and what comes
out the other end.

## Run it

```bash
pixi install
pixi run check     # metadata + naming lint, and the report staleness check
pixi run docs      # regenerate explain_codebase.html and knowledge/
```

Start at `explain_codebase.html` — it is generated from the scripts themselves,
so unlike a README it cannot drift from what the code does. `knowledge/overview.md`
is the ten-line version for a new reader.

## Your own data

```bash
pixi run -e laptop ctx --root ~/ctx_four_runs guide
```

writes the guide page, `outputs/09_report/<date>_report_guide_v1.html`: every command from a
study's accessions, through nf-core/fetchngs and nf-core/rnaseq on a laptop or on Kebnekaise,
to its own course pages and kit. `--dataset-file` writes it for your own dataset file, started
from `conf/dataset_template.toml`. Every kit's `START.html` links the page.

## Layout

| Path | Holds |
|---|---|
| `scripts/<lang>/` | every runnable step |
| `settings/` | paths, figure overrides, canonical output formats |
| `outputs/NN_<step>/` | intermediate results, one directory per step |
| `deliverables/` | dated, versioned things that leave the project |
| `knowledge/` | orientation for humans and agents |

## Provenance

- Scaffold: the lab Copier template. Update with `copier update --skip-tasks`.
- Shared logic: `labcore`, pinned by tag in `pixi.toml`. Never copied in.
- Author: Rytis Stakaitis <rytis.stakaitis@umu.se>
- License: MIT
- Repository: https://github.com/aramburu-lab/cancer_transcriptomics_teaching
