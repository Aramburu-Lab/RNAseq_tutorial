# Containers

The analysis image for Cancer Transcriptomics Teaching Suite: the pixi `prod` environment plus `labcore`, built
multi-stage from `analysis.Dockerfile` and run on Kebnekaise under Apptainer.

Build it with `containers/build.sh`. Never with a bare `docker build` — see
[the architecture trap](#the-architecture-trap) for what that produces on this laptop.

```bash
containers/build.sh dev --push          # fast loop, Docker Hub
containers/build.sh release --version 0.1.0   # CI / login node only
containers/build.sh --help
```

## Two paths

They do different jobs and must never be confused (D-004).

| | **dev** | **release** |
|---|---|---|
| Built where | this Mac, `--platform linux/amd64` under QEMU | GitHub Actions, native linux-64 |
| Registry | `docker.io/stakaitis/cancer_transcriptomics_teaching-analysis` | `ghcr.io/aramburu-lab/cancer_transcriptomics_teaching-analysis` **and** `oras://ghcr.io/aramburu-lab/cancer_transcriptomics_teaching-analysis-sif` |
| Pulled with | `apptainer pull docker://…` | `apptainer pull oras://…` |
| Reproducible | **No**, and that is fine | **Yes** — the ADR-7 guarantee |
| In `digests.lock` | never | always |
| Cite in a paper | never | yes, digest-pinned |
| Tag | short git SHA (`dirty-<sha>` if the tree is dirty) | semver, from `--version` |

Override the dev namespace with `DEV_REGISTRY=docker.io/<you> containers/build.sh dev`.

The dev path exists because there is no `apptainer` on macOS (D-003), so the only way to test an
image change locally would otherwise be a git push and a CI wait.

## The architecture trap

This Mac is `arm64`. Kebnekaise is `x86_64`. A default `docker build` here produces a `linux/arm64`
image, `apptainer pull` accepts it without a word, and it then fails to `exec` on a compute node
with an error that names neither Docker nor the architecture. It is the single most likely
first-run failure.

`build.sh` therefore hardcodes `--platform linux/amd64` and re-reads the built image's
`.Os/.Architecture` before any push, refusing to publish anything that is not `linux/amd64`.

Emulation costs are real. QEMU runs every x86 binary the *build itself* executes — conda post-link
scripts, aligner self-tests, anything that compiles and then runs a check. Expect roughly 3–10×
slower, and expect occasional segfaults that a native runner never produces. **An image whose build
runs compiled x86 binaries is CI-only.** Do not spend an hour chasing a QEMU crash; push and let the
release path build it natively.

## Why `oras://` and not `docker://` for releases

`apptainer pull docker://…` converts an OCI image to a `.sif` **on the machine doing the pull**, and
that conversion is not byte-deterministic: creation dates and other metadata vary per run, so two
pulls of the identical manifest give two `.sif` files with different `sha256sum`s.
`--reproducible` / `APPTAINER_REPRODUCIBLE` narrows the gap but does not close it.

So the release path converts **once**, in CI, from the pushed manifest digest, and publishes the
resulting `.sif` as an ORAS artifact. `apptainer pull oras://…` then transfers that finished file
verbatim: every pull is byte-identical, and the `sif` digest in `digests.lock` is the number to
check it against.

```bash
apptainer pull a.sif oras://ghcr.io/aramburu-lab/cancer_transcriptomics_teaching-analysis-sif:0.1.0
apptainer pull b.sif oras://ghcr.io/aramburu-lab/cancer_transcriptomics_teaching-analysis-sif:0.1.0
sha256sum a.sif b.sif    # must match each other and the sif row of digests.lock
```

Building on the cluster is not an option: HPC2N does not allow `--fakeroot`.

## Prefer an image you did not build

ADR-8: do not rebuild `samtools`. Standard tools come from BioContainers or Seqera Containers
(`community.wave.seqera.io`). Build your own only for (a) images that carry `labcore` and this
project's analysis environment, and (b) tool combinations that do not exist upstream.

**Caveat:** images produced by the *free hosted Wave service* are **deleted after one week** unless
the build sets `wave.freeze = true` or you mirror the image into your own registry. Anything a
manuscript depends on must be mirrored into the `aramburu-lab` GHCR namespace and recorded in
`digests.lock` first. Wave's on-demand build is also unusable from Kebnekaise compute nodes, which
have no outbound egress — pre-pull on the login node.

## Private `labcore`: SSH, not a password prompt

`pixi.toml` pins `labcore` as a git dependency over `https://`. If that repo is private, the
resolver clones it during the build, and an unattended `docker build` would hang forever on a
credential prompt it has no terminal to answer.

The build uses a **BuildKit SSH mount** (`RUN --mount=type=ssh`), forwarded by `build.sh --ssh`:

```bash
ssh-add -l                       # the agent must actually hold a key with access
containers/build.sh dev --ssh
```

The mount is optional (`type=ssh` defaults to `required=false`), and the Dockerfile rewrites
`https://github.com/` to `git@github.com:` **only when the socket is present** — so a public
`labcore` still builds with no agent, no key and no flag. In CI, forward the deploy key with
`webfactory/ssh-agent` and pass `ssh: default` to `docker/build-push-action`.

The alternative — a PAT in a build ARG — is rejected: build args are visible in `docker history`
of the published image.

## Image shape

- Base `ubuntu:24.04`. **Not alpine**: musl breaks every glibc-linked bioconda binary. **Not
  distroless**: Nextflow's container contract requires `bash` ≥3.0 and `ps` inside the image, and
  distroless ships neither. Ubuntu 26.04 LTS exists; 24.04 stays the pin because it is the smaller
  glibc jump from Kebnekaise's own Ubuntu 22.04.
- `procps` supplies `ps`, which Nextflow shells out to for process metrics.
- `fonts-liberation` is not cosmetic. No container ships Arial or Helvetica; without a resolvable
  fallback every figure emits dozens of matplotlib `findfont` warnings into the Slurm log.
- The pixi prefix is `/app/.pixi/envs/prod` in **both** stages. conda/pixi environments are not
  relocatable — absolute paths are baked into binaries at link time and Apptainer preserves the
  layout verbatim — so moving the prefix produces a failure at exec time on the cluster rather than
  at build time here.
- `pixi install --locked`, never `--frozen`. `--locked` aborts when `pixi.lock` is stale relative to
  `pixi.toml`; `--frozen` installs the stale lock silently, which is how an image quietly stops
  matching the environment a paper was written against. If the build fails for a stale lock, run
  `pixi install` and commit the new `pixi.lock` — do not switch the flag.
- Build context is the **repo root**, because the Dockerfile copies `pixi.toml` and `pixi.lock` from
  there: `docker buildx build -f containers/analysis.Dockerfile -t <img> .`. `build.sh` handles it.
- Size. The 150–400 MB figure in the build plan is for a *single-tool* image. A full analysis
  environment does not fit in it: the stock `pixi.toml` (polars + pandas + duckdb + scipy +
  matplotlib + plotly) measured **599 MB compressed / 2.8 GB unpacked** on 2026-08-11, and the
  single largest avoidable slice is `pyside6` + `qt6-main` (~45 MB of site-packages plus their Qt
  libraries), pulled in because the dependency is `matplotlib` rather than `matplotlib-base`.
  Nothing here renders to a screen. Drop a tool from `[feature.tools]` into its own image before
  reaching for anything clever — several small images beat one large one both for pull time and for
  Nextflow, which pulls per process.

## On Kebnekaise

`apptainer` is on the login node with no `module load`; `singularity` is a symlink to it. Compute
nodes have no outbound network, so pull on the login node first and point the cache at project
storage, never at the 25 GB home quota.

```bash
export APPTAINER_CACHEDIR=$SNIC_NOBACKUP/cancer_transcriptomics_teaching/apptainer_cache
apptainer pull oras://ghcr.io/aramburu-lab/cancer_transcriptomics_teaching-analysis-sif:0.1.0
apptainer exec cancer_transcriptomics_teaching-analysis-sif_0.1.0.sif bash -lc 'python -c "import labcore"'
```

## Publication freeze

1. `git tag paper-submitted-v1`.
2. Commit `pixi.lock` and `containers/digests.lock`. Every image reference archived with the paper
   is pinned **`@sha256:`**, never a tag — a tag is a mutable pointer and says nothing about what
   ran. `docker://ghcr.io/…@sha256:…`, not `docker://ghcr.io/…:0.1.0`.
3. Archive the actual `.sif` files alongside the code on Zenodo. A registry is not an archive.
4. Set `frozen: true` in `.labtemplate.yml`. `copier update` is never run on the repo again, and
   bugfixes arrive only through a deliberate, documented re-pin.
