"""TODO.md F.2d: ``build.sh release`` keeps every row of digests.lock it does not rebuild.

The script runs whole, in a committed copy of the two container scripts, with docker and
apptainer stubbed on the PATH: nothing is built or pushed. The lock it writes is read back by
``ctx.figures.panel.image``, the reader of the row the Reproduce panels pin by.
"""

import hashlib
import os
import shutil
import subprocess
from pathlib import Path

from conftest import ROOT

from ctx.figures import panel

VERSION = "0.1.0"
ANALYSIS = "ghcr.io/aramburu-lab/cancer_transcriptomics_teaching-analysis"
SIF = f"oras://{ANALYSIS}-sif"
PUSHED = "sha256:" + "12" * 32
SIF_BYTES = b"a stub sif\n"
DOCKER = f"""#!/bin/sh
case "$1 $2" in
  "version "*) echo amd64 ;;
  "image inspect") case "$*" in *Architecture*) echo linux/amd64 ;; *) echo 1048576 ;; esac ;;
  "image ls") echo 1GB ;;
  "buildx imagetools") echo {PUSHED} ;;
esac
"""
# build --reproducible --force <sif> <uri> writes the .sif; push does nothing.
APPTAINER = f"""#!/bin/sh
[ "$1" = build ] && printf '{SIF_BYTES.decode().strip()}\\n' > "$4"
exit 0
"""


def _row(kind: str, version: str, reference: str, digest: str) -> str:
    """One lock row, as both writers print it."""
    return f"{kind}\t{version}\t{reference}\t{digest}\t2026-10-01T00:00:00Z\told1234"


def _repo(tmp_path: Path, rows: list[str]) -> tuple[Path, Path]:
    """A committed copy of the container scripts and a lock holding ``rows``; (repo, stub bin)."""
    repo, stubs = tmp_path / "repo", tmp_path / "bin"
    (repo / "containers").mkdir(parents=True)
    stubs.mkdir()
    for name in ("build.sh", "lock.sh"):
        shutil.copy(ROOT / "containers" / name, repo / "containers" / name)
    template = (ROOT / "containers/digests.lock").read_text()
    (repo / "containers/digests.lock").write_text(template + "\n".join(rows) + "\n")
    for name in ("pixi.toml", "pixi.lock"):
        (repo / name).write_text("")
    for name, text in (("docker", DOCKER), ("apptainer", APPTAINER)):
        (stubs / name).write_text(text)
        (stubs / name).chmod(0o755)
    git = ["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t"]
    subprocess.run([*git, "init", "-q"], check=True)
    subprocess.run([*git, "add", "-A"], check=True)
    subprocess.run([*git, "commit", "-qm", "fixture"], check=True)
    return repo, stubs


def _release(repo: Path, stubs: Path) -> None:
    """``build.sh release --version VERSION`` in the copy, the stubs first on the PATH."""
    env = {**os.environ, "PATH": f"{stubs}{os.pathsep}{os.environ['PATH']}"}
    argv = ["bash", "containers/build.sh", "release", "--version", VERSION]
    done = subprocess.run(argv, cwd=repo, env=env, capture_output=True, text=True, check=False)
    assert done.returncode == 0, done.stderr[-2000:]


def test_release_replaces_its_rows_and_keeps_every_other_line(tmp_path):
    """A release of 0.1.0 replaces the analysis image's oci and sif rows at 0.1.0 and keeps the
    figures image's row, an earlier release's row and every comment byte for byte; the panel
    then pins the figures image by the kept digest. Catches the lock rewritten wholesale (the
    figures row dropped, F.2d), rows of every version replaced, and a row appended twice."""
    figures = _row("oci", VERSION, panel.IMAGE, "sha256:" + "fe" * 32)
    earlier = _row("oci", "0.0.9", ANALYSIS, "sha256:" + "09" * 32)
    rows = [earlier, _row("oci", VERSION, ANALYSIS, "sha256:" + "aa" * 32), figures,
            _row("sif", VERSION, SIF, "sha256:" + "bb" * 32)]  # fmt: skip
    repo, stubs = _repo(tmp_path, rows)
    lock = repo / "containers/digests.lock"
    before = lock.read_text().splitlines()
    _release(repo, stubs)
    after = lock.read_text().splitlines()
    assert figures in after and earlier in after, after
    assert [line for line in before if line.startswith("#")] == [
        line for line in after if line.startswith("#")
    ]
    data = [line.split("\t") for line in after if line and not line.startswith("#")]
    keys = [tuple(r[:3]) for r in data]
    assert len(keys) == len(set(keys)) == 4, after
    by_key = {tuple(r[:3]): r[3] for r in data}
    assert by_key["oci", VERSION, ANALYSIS] == PUSHED
    assert by_key["sif", VERSION, SIF] == f"sha256:{hashlib.sha256(SIF_BYTES).hexdigest()}"
    assert panel.image(VERSION, lock) == (f"{panel.IMAGE}:{VERSION}@{'sha256:' + 'fe' * 32}", True)


def test_release_into_a_lock_with_no_rows_writes_them_under_its_header(tmp_path):
    """The lock as committed today, no rows, gains the two rows below the column line. Catches
    a writer that keeps rows but loses the header, and one that writes nothing on a first
    release."""
    repo, stubs = _repo(tmp_path, [])
    lock = repo / "containers/digests.lock"
    before = lock.read_text()
    _release(repo, stubs)
    after = lock.read_text()
    assert after.startswith(before)
    added = after[len(before) :].splitlines()
    assert [r.split("\t")[:3] for r in added] == [["oci", VERSION, ANALYSIS], ["sif", VERSION, SIF]]
