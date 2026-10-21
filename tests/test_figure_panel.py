"""Phase 2D: the Reproduce panel, generated from each figure's FigureSpec (ctx.figures.panel).

The downloaded script of one panel is run as a reader would run it, in a clean folder, from the
session pipeline's root; it checks its inputs and the redrawn data table itself, and the test
compares that table with the module's byte for byte. The Docker tab's run inside the figures
image is containers/prove_figures.py's, which CI runs (.github/workflows/containers.yml). Each
test names its mutation.
"""

import html
import importlib.util
import os
import re
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

from conftest import ROOT, module_page
from test_figures import specs

from ctx.figures import panel
from ctx.figures.params import params_file
from ctx.figures.spec import parse
from ctx.report.builder import unpacked

PANEL = re.compile(
    r'<details class="repro" data-script-head="([^"]*)">(.*?)<details class="repro-upstream"', re.S
)
CODE = re.compile(r'<code id="ctx-repro-\d+-([a-z]+)"[^>]*>(.*?)</code>', re.S)
TAB = re.compile(r'role="tab" id="ctx-repro-tab-(\d+)-([a-z]+)"[^>]*>([^<]*)</button>')
DIGEST = "sha256:" + "ab" * 32


def panels(page: str, tab: str = "pixi") -> list[tuple[str, str]]:
    """(pasted block, downloaded script) of one tab of every Reproduce panel on a page.

    Joined as the page's script joins them, which ``test_reproduce_browser`` holds to the
    text Chromium copies and downloads: the tab's step 1, a blank line, the shared steps.
    """
    found = []
    for head, inside in PANEL.findall(page):
        codes = {key: html.unescape(code) for key, code in CODE.findall(inside)}
        pasted = f"{codes[tab]}\n{codes['shared']}"
        found.append((pasted, html.unescape(head) + pasted))
    return found


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
        alone = panels(unpacked((next(root.glob(f"outputs/*/{spec.stem}.html"))).read_text()))
        assert len(alone) == 1 and alone[0] in shown, spec.stem
        assert all(line in alone[0][0] for line in own), spec.stem


def test_tabs_are_docker_apptainer_pixi_and_differ_in_the_environment_alone(pipeline):
    """Every panel shows three tabs, Docker first and selected, each with its own block and
    script, Docker's and Apptainer's step 1 each running ctx through its own tool alone; from
    step 2 on the three blocks are one text. Catches a tab dropped or reordered, a tab's script
    holding another tab's block (every tab printing Docker's step 1 survived without the tool
    check, F.2d), and a step after the environment drifting between tabs."""
    page = module_page(pipeline[0], "05_de", "desynth_overview")
    tabs = TAB.findall(page)
    by_figure: dict[str, list[tuple[str, str]]] = {}
    for index, key, label in tabs:
        by_figure.setdefault(index, []).append((key, label))
    assert by_figure and all(v == [("docker", "Docker"), ("apptainer", "Apptainer"),
                                   ("pixi", "pixi")] for v in by_figure.values())  # fmt: skip
    selected = re.findall(r'id="ctx-repro-tab-\d+-([a-z]+)"[^>]*aria-selected="true"', page)
    assert selected == ["docker"] * len(by_figure)
    blocks = {key: panels(page, key) for key in panel.TABS}
    assert len({len(b) for b in blocks.values()}) == 1
    runs = {"docker": "ctx() { docker run ", "apptainer": "ctx() { apptainer run "}
    for pasted_scripts in zip(*blocks.values(), strict=True):
        rests = set()
        for tab, (pasted, script) in zip(panel.TABS, pasted_scripts, strict=True):
            assert script == f"{panel.SCRIPT_HEAD}\n{pasted}"
            step_1, rest = pasted.split("\n# 2. ", 1)
            rests.add(rest)
            assert all((run in step_1) == (key == tab) for key, run in runs.items()), tab
        assert len(rests) == 1


def test_container_tabs_run_ctx_in_the_image_as_the_reader(pipeline):
    """The Docker tab runs ctx in the image as the reader's own user, on linux/amd64, in the
    clean folder; Apptainer runs the image's entrypoint from a SIF pulled once; both name the
    image by tag with no digest in the lock. Catches a root-owned redraw (-u dropped), exec in
    place of run (ctx off the PATH), and a digest typed in before one is recorded."""
    spec = specs(pipeline[0])[0]
    docker, apptainer = panel.body(spec, "docker"), panel.body(spec, "apptainer")
    tagged = f"{panel.IMAGE}:{spec.version}"
    for block in (docker, apptainer):
        assert f'IMG="${{CTX_FIGURES_IMAGE:-{tagged}}}"\n' in block
        assert "@sha256:" not in block and "named by its tag" in block
        assert "pixi" not in block.split("\n# 2. ")[0]
    assert ('ctx() { docker run --rm --platform linux/amd64 -u "$(id -u):$(id -g)" '
            '-v "$PWD":/work -w /work "$IMG" ctx "$@"; }\n') in docker  # fmt: skip
    assert 'ctx() { apptainer run --cleanenv "$SIF" ctx "$@"; }\n' in apptainer
    assert 'apptainer pull "$SIF" "docker://$IMG"' in apptainer


def test_container_tabs_clone_only_for_an_input_in_the_repository(pipeline):
    """A figure drawn from a checkout with its dataset file in the repository clones that
    commit in the container tabs too, without pixi; one drawn from a wheel clones nothing.
    Catches the clone dropped, which leaves $CTX_REPO empty, and pixi run in a container tab."""
    spec = specs(pipeline[0])[0]
    commit = "0123456789abcdef0123456789abcdef01234567"
    assert any(i.base == "repo" for i in spec.inputs)
    checked_out = panel.body(replace(spec, commit=f"git {commit}"), "docker")
    assert f"git -C {panel.CLONE} checkout {commit}\n" in checked_out
    assert "pixi install" not in checked_out and 'CTX_REPO="$PWD/' in checked_out
    wheel = panel.body(replace(spec, commit="installed wheel"), "apptainer")
    assert "git clone" not in wheel


def test_image_is_pinned_by_digest_once_the_lock_records_it(pipeline, tmp_path, monkeypatch):
    """The lock's oci row of the figures image at ctx's version pins it; the analysis image's
    row, another version's and a row with no digest do not. Catches the analysis image's digest
    on the figures panel, a digest from another release, and the template's commented example."""
    spec = specs(pipeline[0])[0]
    other = "ghcr.io/aramburu-lab/cancer_transcriptomics_teaching-analysis"
    rows = [f"oci\t{spec.version}\t{other}\tsha256:{'cd' * 32}\tT\tr",
            f"oci\t9.9.9\t{panel.IMAGE}\tsha256:{'ef' * 32}\tT\tr",
            f"#   oci\t{spec.version}\t{panel.IMAGE}\t{DIGEST}\tT\tr",
            f"oci\t{spec.version}\t{panel.IMAGE}\tskipped\tT\tr"]  # fmt: skip
    lock = tmp_path / "digests.lock"
    lock.write_text("\n".join(rows) + "\n")
    assert panel.image(spec.version, lock) == (f"{panel.IMAGE}:{spec.version}", False)
    assert panel.image(spec.version, None) == (f"{panel.IMAGE}:{spec.version}", False)
    lock.write_text("\n".join([*rows, f"oci\t{spec.version}\t{panel.IMAGE}\t{DIGEST}\tT\tr"]))
    pinned = f"{panel.IMAGE}:{spec.version}@{DIGEST}"
    assert panel.image(spec.version, lock) == (pinned, True)
    monkeypatch.setattr(panel, "_lock", lambda: lock)
    docker = panel.body(spec, "docker")
    assert f'IMG="${{CTX_FIGURES_IMAGE:-{pinned}}}"' in docker and "pinned by its digest" in docker


def test_the_image_proof_reads_the_docker_tab_the_page_prints(pipeline):
    """containers/prove_figures.py, which CI runs inside the image, finds the volcano's Docker
    tab on M5's page, its script byte for byte, and the stem its data check names; the workflow
    publishes and locks the repository the panel names. Catches a template change the CI proof
    no longer parses, and the image renamed in one of the two places."""
    root = pipeline[0]
    spec = next(s for s in specs(root) if parse(s.figure_id)[2] == "volcano")
    loader = importlib.util.spec_from_file_location("prove", ROOT / "containers/prove_figures.py")
    prove = importlib.util.module_from_spec(loader)
    loader.loader.exec_module(prove)
    script, stem = prove.volcano_tab(root)
    page = module_page(root, "05_de", "desynth_overview")
    assert (script, stem) == (
        next(s for p, s in panels(page, "docker") if f"ctx figure {spec.figure_id} " in p),
        spec.stem[11:],
    )
    workflow = (ROOT / ".github/workflows/containers.yml").read_text()
    assert f"  FIGURES_REPO: {panel.IMAGE}\n" in workflow
    assert '"$FIGURES_REPO" "${{ needs.figures.outputs.digest }}"' in workflow


def test_the_downloaded_script_redraws_the_data_byte_for_byte(pipeline, tmp_path):
    """The volcano's downloaded pixi script, run in a clean folder, copies and checks its
    inputs, redraws the figure and checks the data table's SHA-256 itself; the table equals
    M5's byte for byte. Catches a copy line under the wrong base, a hash of another file than
    the data table, and an input the panel leaves out."""
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
    """No pasted block of any tab holds set -e, which would close the reader's shell on the
    first error; every script starts with the shebang and set -euo pipefail. Catches either
    line moved into the other block."""
    page = module_page(pipeline[0], "05_de", "desynth_overview")
    for tab in panel.TABS:
        for pasted, script in panels(page, tab):
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
    every = [panel.body(replace(spec, commit=c), t) for c in (f"git {commit}", "installed wheel")
             for t in panel.TABS]  # fmt: skip
    assert str(Path.home()) not in clean + dirty + wheel + "".join(every)
