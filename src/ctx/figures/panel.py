"""The Reproduce panel, generated from a ``FigureSpec``, never typed.

Five steps, each commented from ``content/figure.toml``: the environment, the inputs copied
and checked, the params file as a heredoc, ``ctx figure``, and the check of the redrawn data
table against the page's. The pasted block has no ``set -e``, which would close the reader's
shell on the first error; the downloaded script starts with ``set -euo pipefail``.
"""

import base64
import shlex
from typing import Any

from ctx.figures.params import params_file
from ctx.figures.spec import FigureSpec, Input
from ctx.report import text

PAGE = "figure"
# The course repository the environment step clones (git remote of this checkout; its
# releases page is https://github.com/Aramburu-Lab/RNAseq_tutorial/releases, TODO.md K.11).
REPOSITORY = "https://github.com/Aramburu-Lab/RNAseq_tutorial.git"
CLONE = "ctx_course"
PARAMS = "params.toml"
SCRIPT_HEAD = "#!/usr/bin/env bash\nset -euo pipefail\n"


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


def _environment(spec: FigureSpec) -> list[str]:
    """Step 1: ctx at the drawing commit through pixi, or an installed ctx; a clean folder."""
    folder = f"repro_{spec.figure_id}_{spec.dataset}"
    clean = [f"mkdir -p {folder} && cd {folder}"]
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


def body(spec: FigureSpec) -> str:
    """The five steps, as the pasted block holds them."""
    csv = f"????-??-??_{spec.stem[11:]}.csv"
    lines = [
        *_environment(spec), "",
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


def script(spec: FigureSpec) -> str:
    """The downloadable script: the same steps under ``set -euo pipefail``."""
    return f"{SCRIPT_HEAD}\n{body(spec)}"


def block(spec: FigureSpec) -> dict[str, str]:
    """What the page's Reproduce panel prints for one figure."""
    name = f"reproduce_{spec.figure_id}_{spec.dataset}.sh"
    return {
        "summary": text.word(PAGE, "panel.summary"),
        "lead": text.fill(PAGE, "panel.lead", figure_id=spec.figure_id),
        "note": text.word(PAGE, "panel.no_container"),
        "pasted": body(spec),
        "script": base64.b64encode(script(spec).encode()).decode(),
        "name": name,
        "download": text.fill(PAGE, "panel.download", name=name),
        "upstream_title": text.word(PAGE, "panel.upstream"),
        "upstream": spec.upstream,
    }
