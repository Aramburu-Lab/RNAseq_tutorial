# Agent instructions — Cancer Transcriptomics Teaching Suite

Canonical rules for any coding agent here. `CLAUDE.md` imports this file, so there is one
rule set and not two. Read it before touching a file.

`CONVENTIONS.md` holds this project's own rules. Where `CONVENTIONS.md` conflicts with `AGENTS.md`, `CONVENTIONS.md` wins; each divergence is a D-entry in `knowledge/decisions.md`.

## Commands you cannot guess

```bash
pixi install                # build the environment; CI adds --locked
pixi run test               # pytest, the one command that runs this project's tests
pixi run lint               # ruff over scripts/ (preflight runs both)
pixi run docs               # labdocs render && labdocs api
pixi run check              # labdocs lint && labdocs render --check
pixi run hooks              # every prek hook, incl. the file-length ratchet, as CI runs them
pixi run preflight          # everything CI runs; before every push, not after
labdocs render              # regenerate the codebase report: HTML + manifest + map
labdocs lint                # metadata and naming violations
copier update --trust --skip-tasks  # take template updates
```

`--skip-tasks` is mandatory: without it the update re-runs `git init` and re-scaffolds
over your customised tree. Environments: `default` (analysis + dev), `laptop` (analysis
only), `prod` (adds linux-64 cluster tools) — select with `pixi run -e <name>`.

## Never hand-edit these — they are generated

An edit is lost on the next run, and until then the file disagrees with its source.

- `.copier-answers.yml` → change answers by re-running `copier update --trust --skip-tasks`
- `explain_codebase.html` → `labdocs render`
- `codebase_manifest.json` → `labdocs render`
- `knowledge/api_index.md` → `labdocs api`
- `knowledge/codebase_map.md` → `labdocs render`
- `settings/site.env`, if present → generated one-way from `settings/paths.toml`

## Three sync channels, disjoint ownership

No file is ever owned by two channels. That invariant is what keeps updates mergeable.

| Channel | Owns | Update by |
|---|---|---|
| Copier | `AGENTS.md`, `CLAUDE.md`, `knowledge/`, `settings/`, `bin/`, lint config, `.labtemplate.yml` | `copier update --trust --skip-tasks` |
| Package pin | every reusable helper — `labcore.*` | bump the tag in `pixi.toml`, then `pixi update` |
| nf-core | nothing — this is an analysis project | `nf-core pipelines sync` |

So: **do not add a shared helper module to this repo.** A function a second project would
want belongs in `labcore`, where a fix ships as a version bump instead of a scaffold
update nobody asked for. Project-specific helpers stay in `scripts/`, marked `exempt:`.

## Where things live

`scripts/<lang>/` every runnable step, one language per subdirectory ·
`bin/` sourceable bash helpers, `lib_slurm.sh` is the bash channel ·
`settings/` paths, figure overrides, canonical formats, report CSS ·
`outputs/NN_<step>/` intermediates, one directory per step ·
`deliverables/` anything leaving the project ·
`knowledge/` orientation for humans and agents, half generated ·
`logs/` Slurm writes here; must exist before submission

## Naming — three regimes, enforced by `labdocs lint`

| Regime | Applies to | Pattern | Example |
|---|---|---|---|
| Script | `scripts/**` | `<single_string>.<ext>` | `alignment.sh` |
| Step output | `outputs/**` | `out_<script>_<desc>.<ext>` | `out_alignment_unmapped.sam` |
| Deliverable | `deliverables/**`, anything leaving the project | `YYYY-MM-DD_description_vN.ext` | `2026-08-11_de_volcano_v2.pdf` |
| nf-core | `results/**` | nf-core native — **exempt**, never rename | `results/star_salmon/…` |

Figures append a reserved theme token: `out_<script>_<desc>_white.<ext>` or `_black`,
the only two legal (D-005a). A CSV sidecar carries none.
Underscores only, no spaces, no date or `v2` in a script name. A script's `name:` must
equal its filename stem, so rename both together. Every script under `scripts/` carries one
`codebase-meta` block; a helper that is not a pipeline step carries `exempt: <why>` instead,
and every input needs `from: <script>` or `external: true` — silence fails lint.
## Style

`ruff.toml` owns the mechanical half: ≤50 statements and ≤12 branches per function,
complexity ≤10, ≤6 arguments, Google docstrings on public functions, no bare `except`, no
commented-out code. The file-length ratchet lives in `.labtemplate.yml:file_length` and
nowhere else — read the two numbers there rather than trusting a figure quoted here, which
goes stale the first time a project ratchets. Over the limit is the signal to extract a
helper. Comments say **why**, never what; no banner comments; no bare
`TODO`, write `TODO(name): reason` or `ponytail:` for a deliberate deferral.

What no linter checks, and so the first review question: **is the name intuitive?**
`unmapped_reads` not `df2`; `drop_low_quality_cells()` not `process()`.
