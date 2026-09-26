@AGENTS.md

# Claude-specific notes

Every rule is in `AGENTS.md`, imported by the line above (an import, not a symlink —
symlinks need admin or developer mode on Windows). Only what is true for Claude Code
specifically belongs here; write anything else in `AGENTS.md`, where every other agent
and every human reading the repo cold will also see it.

- Verbose domain knowledge — assay protocols, organism quirks, one-off recipes — goes in
  `.claude/skills/<topic>/SKILL.md`, which loads on demand. This file and `AGENTS.md` are
  loaded on every turn, so anything added here is paid for on every turn.
- Before writing a helper, read `knowledge/api_index.md`. It indexes this project's
  functions *and* the installed `labcore`, so "does this exist already?" is one file read
  rather than a glob over the tree.
- `knowledge/codebase_map.md` is the run order. Read it before proposing where a new step
  goes; the `order:` field decides sequence, never the filename.
- Run `pixi run docs` after changing any `codebase-meta` block. A stale report is rejected
  by the prek hook.
- Record any correction you had to make twice as one line in `knowledge/lessons.md`. That
  file is yours to append to and template updates skip it.

@CONVENTIONS.md

Where `CONVENTIONS.md` conflicts with `AGENTS.md`, `CONVENTIONS.md` wins; each divergence is a D-entry in `knowledge/decisions.md`. `knowledge/api_index.md` does not index the `src/ctx` package; read that package directly.
