"""F.2d: plant F2's named mutations and the lock writer's, run the tests named to catch each.

Usage: <env wrapper> python docs/plans/logs/2026-10-21_f2d_mutations.py <scratch dir> [name ...]

Each pytest session rebuilds the session pipeline (about 10 min), so the mutations run in
parallel, one copy of this worktree per worker under <scratch dir>, each with its own src on
PYTHONPATH. A copy keeps the worktree's .git file, so runtime.source() reads the same commit;
GIT_OPTIONAL_LOCKS=0 keeps their git status from writing the shared index. Each mutation
replaces one exact text, found exactly once, and is caught when its tests fail. One extra
worker runs every named test on an unmutated copy: the baseline.
"""

import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
LOG = Path(__file__).with_suffix(".log")
BIN = "/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin"
WORKERS = 6
P, SAVE, RECORD = "src/ctx/figures/panel.py", "src/ctx/viz/save.py", "src/ctx/figures/record.py"
REPRO, LOCK = "src/ctx/report/templates/_reproduce.html.j2", "containers/lock.sh"
F = "tests/test_figure_panel.py::"
L = "tests/test_containers_lock.py::"
EVERY = F + "test_every_spec_has_its_panel_on_the_module_page"
TABS = F + "test_tabs_are_docker_apptainer_pixi_and_differ_in_the_environment_alone"
RUN = F + "test_container_tabs_run_ctx_in_the_image_as_the_reader"
CLONE = F + "test_container_tabs_clone_only_for_an_input_in_the_repository"
PIN = F + "test_image_is_pinned_by_digest_once_the_lock_records_it"
PROOF = F + "test_the_image_proof_reads_the_docker_tab_the_page_prints"
SCRIPT = F + "test_the_downloaded_script_redraws_the_data_byte_for_byte"
SET_E = F + "test_pasted_block_never_stops_the_shell_and_the_script_does"
ENV = F + "test_environment_follows_where_ctx_came_from"
KEEP = L + "test_release_replaces_its_rows_and_keeps_every_other_line"
FIRST = L + "test_release_into_a_lock_with_no_rows_writes_them_under_its_header"
KEY = "new[$1 FS $2 FS $3] = 1; next }\n                !(($1 FS $2 FS $3) in new)"

# (name, file, old, new, tests)
MUTATIONS: list[tuple[str, str, str, str, list[str]]] = [
    ("the module page's block without the panel the figure page has", SAVE,
     "reproduce=_panels().get(figure_stem),", "reproduce=None,", [EVERY]),
    ("a panel kept from another figure", SAVE,
     "_panels()[figure_stem] = panel.block(spec)",
     '_panels()[figure_stem] = _panels().setdefault("", panel.block(spec))', [EVERY]),
    ("a tab dropped", P, 'TABS = ("docker", "apptainer", "pixi")',
     'TABS = ("docker", "pixi")', [TABS]),
    ("the tabs reordered", P, 'TABS = ("docker", "apptainer", "pixi")',
     'TABS = ("pixi", "docker", "apptainer")', [TABS]),
    ("every tab's block holding the Docker tab's", REPRO,
     "{{ t.head }}</code>", "{{ r.tabs[0].head }}</code>", [TABS]),
    ("a step after the environment drifting in one tab", P,
     '"head": head(spec, tab),',
     '"head": head(spec, tab) + ("\\n# 2. drift\\n" if tab == "pixi" else ""),', [TABS]),
    ("a root-owned redraw (-u dropped)", P,
     '--platform {PLATFORM} -u "$(id -u):$(id -g)" \'', "--platform {PLATFORM} '", [RUN]),
    ("exec in place of run", P, "ctx() { apptainer run --cleanenv",
     "ctx() { apptainer exec --cleanenv", [RUN]),
    ("a digest typed in before one is recorded", P,
     "f'IMG=\"${{CTX_FIGURES_IMAGE:-{ref}}}\"',",
     "f'IMG=\"${{CTX_FIGURES_IMAGE:-{ref}@sha256:{\"0\" * 64}}}\"',", [RUN]),
    ("the clone dropped", P, "return [*lines, *_clone(spec)]", "return lines", [CLONE]),
    ("pixi run in a container tab", P, 'if tab != "pixi":', 'if tab == "apptainer":', [CLONE]),
    ("the analysis image's digest on the figures panel", P,
     '(kind, row_version, reference) == ("oci", version, IMAGE)',
     '(kind, row_version) == ("oci", version)', [PIN]),
    ("a digest from another release", P,
     '(kind, row_version, reference) == ("oci", version, IMAGE)',
     '(kind, reference) == ("oci", IMAGE)', [PIN]),
    ("the template's commented example read as a row", P, '[*row.split("\\t"), "", "", ""]',
     '[*row.lstrip("# ").split("\\t"), "", "", ""]', [PIN]),
    ("a row with no digest pinning the image", P, 'and digest[:7] == "sha256:"', "", [PIN]),
    ("a template change the CI proof no longer parses", REPRO,
     '<code id="ctx-repro-{{ f.index }}-shared">', '<code id="ctx-repro-{{ f.index }}-steps">',
     [PROOF]),
    ("the image renamed in one of the two places", P,
     'IMAGE = "ghcr.io/aramburu-lab/cancer_transcriptomics_teaching-figures"',
     'IMAGE = "ghcr.io/aramburu-lab/cancer_transcriptomics_teaching-figure"', [PROOF]),
    ("a copy line under the wrong base", P, "return f'cp \"$ROOT/{item.source}\" .'",
     "return f'cp \"$CTX_REPO/{item.source}\" .'", [SCRIPT]),
    ("a hash of another file than the data table", RECORD, "data_sha256=sha256(data_csv),",
     "data_sha256=sha256(next(iter(paths.values()))),", [SCRIPT]),
    ("an input the panel leaves out", P, "*(_copy(i) for i in spec.inputs),",
     "*(_copy(i) for i in spec.inputs[1:]),", [SCRIPT]),
    ("set -e moved into the pasted block", P, '_comment("steps.inputs", step=2), ',
     '"set -euo pipefail", _comment("steps.inputs", step=2), ', [SET_E]),
    ("set -euo pipefail left out of the script", P,
     'SCRIPT_HEAD = "#!/usr/bin/env bash\\nset -euo pipefail\\n"',
     'SCRIPT_HEAD = "#!/usr/bin/env bash\\n"', [SET_E]),
    ("the commit of another checkout typed in", P,
     "    commit = spec.commit.split()[1]\n    key =",
     '    commit = "f" * 40\n    key =', [ENV]),
    ("a dirty tree presented as reproducible", P,
     'key = "steps.environment_dirty" if "uncommitted" in spec.commit else "steps.environment"',
     'key = "steps.environment"', [ENV]),
    ("the lock rewritten without its rows (F.2d)", LOCK, "!(($1 FS $2 FS $3) in new)", "0",
     [KEEP]),
    ("rows of every version replaced", LOCK, KEY,
     "new[$1 FS $3] = 1; next }\n                !(($1 FS $3) in new)", [KEEP]),
    ("a row appended twice", LOCK, "!(($1 FS $2 FS $3) in new)", "1", [KEEP]),
    ("the header lost, rows kept", LOCK, "!(($1 FS $2 FS $3) in new)",
     "!/^#/ && !(($1 FS $2 FS $3) in new)", [FIRST]),
    ("nothing written on a first release", LOCK, "        printf '%s\\n' \"$@\"\n    } >",
     "    } >", [FIRST]),
]  # fmt: skip


def _copy(scratch: Path, index: int) -> Path:
    """A fresh copy of this worktree, bytecode and outputs left out."""
    copy = scratch / f"copy_{index}"
    shutil.rmtree(copy, ignore_errors=True)
    shutil.copytree(ROOT, copy, ignore=shutil.ignore_patterns("__pycache__", "copy_*"))
    return copy


def _run(copy: Path, tests: list[str]) -> int:
    """The tests' exit code, run serially in one copy with that copy's src."""
    env = {**os.environ, "PATH": f"{BIN}:{os.environ['PATH']}", "PYTHONPATH": str(copy / "src"),
           "PYTHONDONTWRITEBYTECODE": "1", "GIT_OPTIONAL_LOCKS": "0"}  # fmt: skip
    argv = [f"{BIN}/python", "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider", *tests]
    return subprocess.run(argv, cwd=copy, env=env, capture_output=True, check=False).returncode


def _plant(copy: Path, mutation: tuple[str, str, str, str, list[str]]) -> str:
    """Plant one mutation in a copy, run its tests, restore; the log line."""
    name, file, old, new, tests = mutation
    path = copy / file
    source = path.read_text()
    assert source.count(old) == 1, (name, source.count(old))
    path.write_text(source.replace(old, new))
    try:
        code = _run(copy, tests)
    finally:
        path.write_text(source)
    verdict = "caught" if code != 0 else "SURVIVED"
    return f"{verdict} (exit {code}): {name} [{' '.join(t.split('::')[-1] for t in tests)}]"


def _worker(scratch: Path, index: int, chosen: list) -> list[str]:
    """One copy: its share of the mutations, or the baseline when it has none."""
    copy = _copy(scratch, index)
    if not chosen:
        tests = sorted({t for m in MUTATIONS for t in m[4]})
        return [f"baseline, every named test on an unmutated copy: exit {_run(copy, tests)}"]
    return [_plant(copy, m) for m in chosen]


def main(scratch: Path, names: list[str]) -> int:
    """Run the chosen mutations (all when none is named); return how many survived."""
    chosen = [m for m in MUTATIONS if not names or m[0] in names]
    shares = [share for i in range(WORKERS) if (share := chosen[i::WORKERS])] + [[]]
    with ThreadPoolExecutor(len(shares)) as pool:
        lines = [line for done in pool.map(lambda a: _worker(scratch, *a), enumerate(shares))
                 for line in done]  # fmt: skip
    with LOG.open("a") as log:
        for line in lines:
            print(line)
            log.write(line + "\n")
    return sum(line.startswith("SURVIVED") for line in lines)


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1]), sys.argv[2:]))
