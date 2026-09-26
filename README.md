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
