"""Phase 2D: the Reproduce panel, generated from each figure's FigureSpec (ctx.figures.panel).

The downloaded script of one panel is run as a reader would run it, in a clean folder, from the
session pipeline's root; it checks its inputs and the redrawn data table itself, and the test
compares that table with the module's byte for byte. Each test names its mutation.
"""

import base64
import html
import os
import re
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

from conftest import module_page
from test_figures import specs

from ctx.figures import panel
from ctx.figures.params import params_file
from ctx.figures.spec import parse

PANEL = re.compile(
    r'<code id="ctx-repro-(\d+)">(.*?)</code>.*?href="data:text/x-shellscript;base64,([^"]+)"',
    re.S,
)


def panels(page: str) -> list[tuple[str, str]]:
    """(pasted block, downloaded script) of every Reproduce panel on a page."""
    return [(html.unescape(p), base64.b64decode(s).decode()) for _, p, s in PANEL.findall(page)]


def test_every_spec_has_its_panel_on_the_module_page(pipeline):
    """M5's and M6's pages carry one panel per figure with a spec, each the panel its spec
    generates. Catches the module page's block, built in a second call, without the panel the
    figure page has, and a panel typed or kept from another figure."""
    root = pipeline[0]
    found = specs(root)
    shown = [*panels(module_page(root, "05_de", "desynth_overview")),
             *panels(module_page(root, "06_function", "desynth_overview"))]  # fmt: skip
    assert len(shown) == len(found) > 0
    for spec in found:
        # The pipeline drew from the wheel, this process from the checkout: the environment
        # step and the dataset file's copy line differ; the figure's own lines may not.
        csv = f"????-??-??_{spec.stem[11:]}.csv"
        own = [params_file(spec), f"ctx figure {spec.figure_id} --params params.toml\n",
               f"printf '%s  %s\\n' {spec.data_sha256} {csv}"]  # fmt: skip
        alone = panels((next(root.glob(f"outputs/*/{spec.stem}.html"))).read_text())
        assert len(alone) == 1 and alone[0] in shown, spec.stem
        assert all(line in alone[0][0] for line in own), spec.stem


def test_the_downloaded_script_redraws_the_data_byte_for_byte(pipeline, tmp_path):
    """The volcano's downloaded script, run in a clean folder, copies and checks its inputs,
    redraws the figure and checks the data table's SHA-256 itself; the table equals M5's byte
    for byte. Catches a copy line under the wrong base, a hash of another file than the data
    table, and an input the panel leaves out."""
    root = pipeline[0]
    spec = next(s for s in specs(root) if parse(s.figure_id)[2] == "volcano")
    page = module_page(root, "05_de", "desynth_overview")
    (script,) = [s for p, s in panels(page) if f"ctx figure {spec.figure_id} " in p]
    clean = tmp_path / "reader"
    clean.mkdir()
    (clean / "reproduce.sh").write_text(script)
    # ctx as the reader's environment step leaves it: the command on the PATH.
    env = {**os.environ, "PATH": f"{Path(sys.executable).parent}{os.pathsep}{os.environ['PATH']}"}
    env = {k: v for k, v in env.items() if not k.startswith("CTX_")}
    done = subprocess.run(["bash", "reproduce.sh"], cwd=clean, env=env, capture_output=True,
                          text=True, check=False)  # fmt: skip
    assert done.returncode == 0, (done.stdout[-2000:], done.stderr[-3000:])
    assert "data.sha256" not in done.stderr and ": OK" in done.stdout, done.stdout
    (csv,) = clean.glob(f"repro_*/????-??-??_{spec.stem[11:]}.csv")
    assert csv.read_bytes() == next(root.glob(f"outputs/05_de/{spec.stem}.csv")).read_bytes()


def test_pasted_block_never_stops_the_shell_and_the_script_does(pipeline):
    """No pasted block holds set -e, which would close the reader's shell on the first error;
    every script starts with the shebang and set -euo pipefail. Catches either line moved
    into the other block."""
    root = pipeline[0]
    for pasted, script in panels(module_page(root, "05_de", "desynth_overview")):
        assert not re.search(r"^\s*set -[a-z]*e", pasted, re.M), pasted[:200]
        assert script.startswith("#!/usr/bin/env bash\nset -euo pipefail\n"), script[:80]
        assert script.endswith(pasted)


def test_environment_follows_where_ctx_came_from(pipeline):
    """A git checkout's panel clones and checks out its commit; one with uncommitted changes
    says a clone does not hold them; an installed ctx clones nothing. Catches the commit of
    another checkout typed in, and a dirty tree presented as reproducible."""
    spec = specs(pipeline[0])[0]
    commit = "0123456789abcdef0123456789abcdef01234567"
    clean = panel.body(replace(spec, commit=f"git {commit}"))
    assert f"git -C {panel.CLONE} checkout {commit}\n" in clean
    assert "--frozen" in clean and "uncommitted" not in clean
    dirty = panel.body(replace(spec, commit=f"git {commit} (uncommitted changes)"))
    assert "uncommitted changes, which a clone does not hold" in dirty
    wheel = panel.body(replace(spec, commit="installed wheel"))
    assert "git clone" not in wheel and "an installed ctx" in wheel
    assert str(Path.home()) not in clean + dirty + wheel
