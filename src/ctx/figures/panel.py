"""The Reproduce panel, generated from a ``FigureSpec``, never typed.

Five steps, each commented from ``content/figure.toml``: the environment, the inputs copied
and checked, the params file as a heredoc, ``ctx figure``, and the check of the redrawn data
table against the page's. Three tabs differ in the environment step alone: Docker and
Apptainer run the figures image (``containers/analysis.Dockerfile``, target ``figures``) behind
a shell function named ``ctx``; pixi builds the environment from a clone. The pasted block has
no ``set -e``, which would close the reader's shell on the first error; the downloaded script
starts with ``set -euo pipefail``.
"""

import base64
import shlex
from pathlib import Path
from typing import Any

from ctx.figures.params import params_file
from ctx.figures.record import checkout
from ctx.figures.spec import FigureSpec, Input
from ctx.report import text

PAGE = "figure"
# The course repository the environment step clones (git remote of this checkout; its
# releases page is https://github.com/Aramburu-Lab/RNAseq_tutorial/releases, TODO.md K.11).
REPOSITORY = "https://github.com/Aramburu-Lab/RNAseq_tutorial.git"
CLONE = "ctx_course"
PARAMS = "params.toml"
SCRIPT_HEAD = "#!/usr/bin/env bash\nset -euo pipefail\n"
# The figures image as .github/workflows/containers.yml publishes it, tagged with ctx's
# version; its digest comes from the lock that workflow commits.
IMAGE = "ghcr.io/aramburu-lab/cancer_transcriptomics_teaching-figures"
DIGESTS = Path("containers") / "digests.lock"
# pixi.lock solves linux-64 alone, so the image has no arm64 variant.
PLATFORM = "linux/amd64"
TABS = ("docker", "apptainer", "pixi")


def _comment(key: str, **values: Any) -> str:
    """One comment line of the panel, from the content file."""
    return f"# {text.fill(PAGE, key, **values)}"


def _quoted(path: str) -> str:
    """A shown path quoted for bash, ``~/`` kept expandable as ``$HOME``."""
    if path == "~" or path.startswith("~/"):
        rest = path[1:].replace("\\", "\\\\").replace('"', '\\"').replace("$", "\\$")
        return f'"$HOME{rest}"'
    return shlex.quote(path)


def _copy(item: Input) -> str:
    """The line that copies one input into the clean folder."""
    if item.base == "root":
        return f'cp "$ROOT/{item.source}" .'
    if item.base == "repo":
        return f'cp "$CTX_REPO/{item.source}" .'
    return f"cp {_quoted(item.source)} ."


def image(version: str, lock: Path | None) -> tuple[str, bool]:
    """The figures image of one ctx version, and whether it is pinned by digest.

    Args:
        version: ctx's version, which is the image's tag.
        lock: ``containers/digests.lock``; None or a missing file names the tag alone.

    Returns:
        ``<repository>:<version>@sha256:<digest>`` and True when the lock has an ``oci`` row of
        the image at that version; else ``<repository>:<version>`` and False.
    """
    tagged = f"{IMAGE}:{version}"
    rows = lock.read_text().splitlines() if lock is not None and lock.is_file() else []
    for row in rows:
        kind, row_version, reference, digest = [*row.split("\t"), "", "", ""][:4]
        if (kind, row_version, reference) == ("oci", version, IMAGE) and digest[:7] == "sha256:":
            return f"{tagged}@{digest}", True
    return tagged, False


def _lock() -> Path | None:
    """The digests lock of the checkout ctx runs from; None for an installed ctx."""
    repo = checkout()
    return None if repo is None else repo / DIGESTS


def _clone(spec: FigureSpec) -> list[str]:
    """A clone at the drawing commit, for a container tab with an input in the repository."""
    if not spec.commit.startswith("git ") or all(i.base != "repo" for i in spec.inputs):
        return []
    commit = spec.commit.split()[1]
    return [
        _comment("steps.repository", commit=commit),
        f"git clone {REPOSITORY} {CLONE}",
        f"git -C {CLONE} checkout {commit}",
        f'CTX_REPO="$PWD/{CLONE}"',
    ]


def _container(spec: FigureSpec, tab: str) -> list[str]:
    """Step 1 of the Docker or the Apptainer tab: the image, and ``ctx`` run inside it."""
    ref, pinned = image(spec.version, _lock())
    lines = [
        _comment(f"steps.environment_{tab}", step=1, version=spec.version),
        _comment("steps.image_digest" if pinned else "steps.image_tag"),
        _comment("steps.image_platform"),
        f'IMG="${{CTX_FIGURES_IMAGE:-{ref}}}"',
    ]
    if tab == "docker":
        # -u: the redrawn files belong to the reader, not to root.
        lines += [
            f'docker image inspect "$IMG" >/dev/null 2>&1 || docker pull --platform {PLATFORM} '
            '"$IMG"',
            f'ctx() {{ docker run --rm --platform {PLATFORM} -u "$(id -u):$(id -g)" '
            '-v "$PWD":/work -w /work "$IMG" ctx "$@"; }',
        ]
    else:
        # run, not exec: the image's entrypoint is what puts ctx on the PATH.
        lines += [
            _comment("steps.image_login"),
            f'SIF="$PWD/ctx_figures_{spec.version}.sif"',
            '[ -f "$SIF" ] || apptainer pull "$SIF" "docker://$IMG"',
            'ctx() { apptainer run --cleanenv "$SIF" ctx "$@"; }',
        ]
    return [*lines, *_clone(spec)]


def _environment(spec: FigureSpec, tab: str) -> list[str]:
    """Step 1: the image, ctx at the drawing commit through pixi, or an installed ctx."""
    folder = f"repro_{spec.figure_id}_{spec.dataset}"
    clean = [f"mkdir -p {folder} && cd {folder}"]
    if tab != "pixi":
        return [*_container(spec, tab), *clean]
    if not spec.commit.startswith("git "):
        return [_comment("steps.environment_wheel", step=1, version=spec.version), *clean]
    commit = spec.commit.split()[1]
    key = "steps.environment_dirty" if "uncommitted" in spec.commit else "steps.environment"
    return [
        _comment(key, step=1, version=spec.version, commit=commit),
        f"git clone {REPOSITORY} {CLONE}",
        f"git -C {CLONE} checkout {commit}",
        f"pixi install --manifest-path {CLONE}/pixi.toml --frozen",
        f'CTX_REPO="$PWD/{CLONE}"',
        'ctx() { pixi run --manifest-path "$CTX_REPO/pixi.toml" ctx "$@"; }',
        *clean,
    ]


def body(spec: FigureSpec, tab: str = "pixi") -> str:
    """The five steps of one tab, as the pasted block holds them."""
    csv = f"????-??-??_{spec.stem[11:]}.csv"
    lines = [
        *_environment(spec, tab), "",
        _comment("steps.inputs", step=2), _comment("steps.root"), f"ROOT={_quoted(spec.root)}",
        *(_copy(i) for i in spec.inputs),
        _comment("steps.check_tool"),
        'sums() { if command -v sha256sum >/dev/null; then sha256sum -c "$1"; '
        'else shasum -a 256 -c "$1"; fi; }',
        "cat > inputs.sha256 <<'EOF'", *(f"{i.sha256}  {i.file}" for i in spec.inputs), "EOF",
        "sums inputs.sha256", "",
        _comment("steps.params", step=3), f"cat > {PARAMS} <<'EOF'",
        params_file(spec).rstrip("\n"), "EOF", "",
        _comment("steps.generate", step=4, plotting=spec.plotting, file=spec.plotting_file),
        f"ctx figure {spec.figure_id} --params {PARAMS}", "",
        _comment("steps.check", step=5), "ls",
        f"printf '%s  %s\\n' {spec.data_sha256} {csv} > data.sha256", "sums data.sha256",
    ]  # fmt: skip
    return "\n".join(lines) + "\n"


def script(spec: FigureSpec, tab: str = "pixi") -> str:
    """The downloadable script of one tab: the same steps under ``set -euo pipefail``."""
    return f"{SCRIPT_HEAD}\n{body(spec, tab)}"


def _tab(spec: FigureSpec, tab: str) -> dict[str, str]:
    """One tab of the panel: its label, its pasted block and its downloadable script."""
    name = f"reproduce_{spec.figure_id}_{spec.dataset}_{tab}.sh"
    return {
        "key": tab,
        "label": text.word(PAGE, f"tabs.{tab}"),
        "pasted": body(spec, tab),
        "script": base64.b64encode(script(spec, tab).encode()).decode(),
        "name": name,
        "download": text.fill(PAGE, "panel.download", name=name),
    }


def block(spec: FigureSpec) -> dict[str, Any]:
    """What the page's Reproduce panel prints for one figure, its tabs in ``TABS`` order."""
    return {
        "summary": text.word(PAGE, "panel.summary"),
        "lead": text.fill(PAGE, "panel.lead", figure_id=spec.figure_id),
        "note": text.word(PAGE, "panel.note"),
        "tabs": [_tab(spec, tab) for tab in TABS],
        "upstream_title": text.word(PAGE, "panel.upstream"),
        "upstream": spec.upstream,
    }
