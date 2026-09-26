# Versions — Cancer Transcriptomics Teaching Suite

Where a version actually lives, and the site facts that are not versions of
anything installable.

**Nothing here is authoritative for a package.** `pixi.toml` declares the
constraints and `pixi.lock` records the exact resolution; both are committed.
This file exists for the pins that no lockfile can hold.

## Environment

| Thing | Where the pin lives |
|---|---|
| Python and every conda/PyPI package | `pixi.toml` + `pixi.lock` |
| `labcore` | `pixi.toml`, `[feature.analysis.pypi-dependencies]`, by git tag |
| ruff, prek | `pixi.toml` dev feature, and `rev:` in `.pre-commit-config.yaml` |
| Container images | `containers/digests.lock`, by `@sha256:` digest |

Bump `labcore` by editing the tag and running `pixi update labcore`. Never vendor
a copy of it into this repository to work around a pin.

## Reference data

| Dataset | Version / release | Where it came from |
|---|---|---|
| | | |

A genome build or annotation release belongs here, with its download date. It is
the pin most often forgotten and the one that most often invalidates a result.
