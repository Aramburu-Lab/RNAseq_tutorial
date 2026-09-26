# Analysis image for Cancer Transcriptomics Teaching Suite — pixi `prod` environment, multi-stage.
#
# Build from the REPO ROOT, never from containers/, because the COPY below reads pixi.toml and
# pixi.lock from the root:
#     docker buildx build --platform linux/amd64 -f containers/analysis.Dockerfile -t <img> .
# containers/build.sh does this for you and is the supported entry point.

FROM ghcr.io/prefix-dev/pixi:0.76.2 AS build

WORKDIR /app
COPY pixi.toml pixi.lock ./

# labcore is a git dependency (ADR-4, Channel 3). If the repo is private, the resolver clones it
# over SSH and an unattended build would otherwise block on a credential prompt forever. The agent
# socket is only forwarded when the caller passes `--ssh default`; `required=false` (the default for
# type=ssh) keeps the public-repo case building with no agent at all. The insteadOf rewrite is what
# makes an https:// URL in pixi.toml resolve over that socket, and it is applied only when the
# socket is actually there — otherwise a public build would break trying to use SSH it has no key
# for. accept-new avoids baking a known_hosts file that expires when GitHub rotates its host key.
RUN --mount=type=ssh set -eu; \
    if [ -n "${SSH_AUTH_SOCK:-}" ] && [ -S "${SSH_AUTH_SOCK}" ]; then \
        git config --global core.sshCommand "ssh -o StrictHostKeyChecking=accept-new"; \
        git config --global url."git@github.com:".insteadOf "https://github.com/"; \
    fi; \
    pixi install --locked -e prod

# --locked, never --frozen. --locked aborts when pixi.lock is stale relative to pixi.toml;
# --frozen installs the stale lock and says nothing, which is how a container silently stops
# matching the environment the paper was written against.

# The shell-hook is the activation script pixi would run interactively: PATH, LD_LIBRARY_PATH and
# every activation script the conda packages ship. Baking it into an entrypoint gives the runtime
# stage a correct environment without shipping pixi itself.
RUN pixi shell-hook -e prod -s bash > /shell-hook
RUN printf '#!/bin/bash\n' > /app/entrypoint.sh \
 && cat /shell-hook >> /app/entrypoint.sh \
 && echo 'exec "$@"' >> /app/entrypoint.sh


FROM ubuntu:24.04 AS production

# ubuntu:24.04, deliberately. Not alpine — musl breaks every glibc-linked bioconda binary. Not
# distroless — Nextflow's container contract requires bash >= 3.0 and `ps` inside the image, and
# distroless ships neither. Ubuntu 26.04 LTS exists as of 2026-04-23; 24.04 stays the pin because
# it is the smaller glibc jump from Kebnekaise's own Ubuntu 22.04 (knowledge/versions.md).
#
# procps is `ps`, which Nextflow shells out to for process metrics.
# fonts-liberation is not cosmetic: no container ships Arial or Helvetica, and without a resolvable
# fallback every figure emits dozens of matplotlib findfont warnings into the Slurm log
# (knowledge/prior_art.md §2).
RUN apt-get update \
 && apt-get install -y --no-install-recommends bash procps fonts-liberation \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# THE PREFIX PATH MUST BE IDENTICAL IN BOTH STAGES. pixi/conda environments are not relocatable:
# absolute paths are baked into binaries and .pc/.cmake files at link time, and Apptainer preserves
# the image layout verbatim, so a moved prefix fails at exec time on the cluster rather than here.
COPY --from=build /app/.pixi/envs/prod /app/.pixi/envs/prod
COPY --from=build --chmod=0755 /app/entrypoint.sh /app/entrypoint.sh

# Written by containers/build.sh so `docker inspect` can answer "which commit built this".
ARG VCS_REF=unknown
ARG IMAGE_VERSION=dev
LABEL org.opencontainers.image.title="cancer_transcriptomics_teaching-analysis" \
      org.opencontainers.image.description="Cancer Transcriptomics Teaching Suite analysis environment" \
      org.opencontainers.image.source="https://github.com/aramburu-lab/cancer_transcriptomics_teaching" \
      org.opencontainers.image.revision="${VCS_REF}" \
      org.opencontainers.image.version="${IMAGE_VERSION}" \
      org.opencontainers.image.licenses="MIT"

ENTRYPOINT ["/app/entrypoint.sh"]
CMD ["bash"]
