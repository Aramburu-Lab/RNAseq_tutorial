# check_layout: move or exempt, one line each, for Rytis to approve

Template v0.9.16's `bin/check_layout.py`, run on 2026-11-05 against this working tree. The run
used `python3 <worktree>/bin/check_layout.py /Users/rytis/Scripts/Pipeline/RNAseq_tutorial`, and
it found 5 problems. The template update waits on branch `template-v0.9.16` (ba5ea83) and is not
merged. Nothing has been moved, exempted, committed or pushed.

Mark each row approved, or give it another verdict.

| # | Path | Tracked | Verdict | What reads or writes it today | References a MOVE would change |
|---|---|---|---|---|---|
| 1 | `src/` | yes | **EXEMPT**: the `ctx` package. CONVENTIONS.md, Stack: "CLI package `ctx` (typer)"; Code rules: "`interrogate` at 100 % for `src/`", "`mypy` clean on `src/`" | `pyproject.toml:21` `packages = ["src/ctx"]` and `:25` `files = ["src"]`; the pixi editable install of `.`; `.pre-commit-config.yaml` mypy-src and interrogate-src, `files: '^src/.*\.py$'`; `containers/analysis.Dockerfile:58-72` and its `.dockerignore` `!src/`; `scripts/python/build_kit.py:103` `CODE`; `make_fetchngs_fixture.py:63`, `make_quant_fixture.py:31`, `make_rnaseq_fixture.py:44`; every test through `import ctx` | n/a |
| 2 | `docs/` | yes: 26 plans in `docs/plans/`, logs under `docs/plans/logs/` | **EXEMPT**: the plans and their run logs, the evidence D-entries and TODO.md cite. CONVENTIONS.md, Data and scope boundaries, Licences: licence terms "documented in `docs/`". Only the pipeline flavour sanctions `docs/` | `tests/test_proteo_tumour_normal.py:375` reads `docs/plans/2026-10-14_cptac_ccrcc_m7_rules.md` through `git log --diff-filter=A`; after a move the first add would be the move commit, so the registration check would compare against post-move text. 122 text references to `docs/plans` sit in 22 files: TODO.md, knowledge/decisions.md, knowledge/lessons.md, 11 dataset tomls, qc_metrics.toml, ci.yml, containers/README.md and 4 tests | n/a. The only sanctioned alternative is `knowledge/`, which copier owns (AGENTS.md, Three sync channels) |
| 3 | `docs/plans/logs/*.log`, 48 preflight logs from 2026-10-22 to 2026-11-04 | no | **EXEMPT** under row 2 | written by the guarded preflight command; nothing reads them | none |
| 4 | `TEMP/`: your drafts v3.10 to v3.12 (`CLAUDE9.md`, `CLAUDE10.md`, `todo9.md`, `todo10.md`, `lessons9.md`, `glossary.yaml`) and `cleanup_older_runs.py`, 312 K | no | **MOVE** to `~/Documents/RNAseq_tutorial_drafts/`, outside the repo. The repo is public and these are private drafts, so they should not be committed | nothing reads it; `src/ctx/report/glossary.yaml` was copied from it once (L.1) | 11 historical mentions (TODO.md lines 868, 869, 989, 990, 3945; knowledge/decisions.md lines 1139, 1245, 2138, 2182, 2210, 2231). Proposal: leave them as records of where the drafts sat then, and add one line to TODO.md naming the new folder |
| 5 | `TODO.html` and `TODO_files/`: Quarto 1.9.38's render of TODO.md, 348 K and 888 K | no | **MOVE** to `~/Documents/RNAseq_tutorial_drafts/`; it can be re-rendered from TODO.md at any time | nothing | none; no file in git names them |

Two items made copier refuse the tree as dirty but are not layout problems:

| Path | Why it is dirty | Verdict |
|---|---|---|
| `.claude/`: `launch.json`, the `ctx-core` preview server I wrote 2026-10-25, and an empty `worktrees/` | untracked | none needed: check_layout skips dot entries. v0.9.16 also adds `.claude/settings.json`, a Stop hook on Claude sessions, which awaits your yes or no |
| `pixi.lock` | pixi 0.81.0 rewrites it on every run (lessons.md, 84b5ec9) | none: discard it, or keep running with `PIXI_FROZEN=true` |

If rows 1 to 3 are approved, `.labtemplate.yml` gets:

```yaml
layout_exempt:
  - docs
  - src
```
