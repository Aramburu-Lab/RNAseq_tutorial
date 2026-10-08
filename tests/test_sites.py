"""The site files: what fetchngs.sbatch and rnaseq.sbatch change per computer (U.4).

The launchers run here behind the stubs of test_cluster_scripts, with uname and docker stubbed
too, so a laptop's arm64 and x86_64 are both tested on any machine. Each test names the
mutation it catches.
"""

import hashlib
import re
import subprocess

import pytest
from test_cluster_scripts import (
    BASH,
    FETCHNGS,
    _after,
    _ids,
    _recorded,
    _rnaseq,
    _run,
    _spool,
    _stub,
)
from test_cluster_scripts import site as _cluster_site

from ctx.gather.fetchngs import FETCHNGS_VERSION
from ctx.gather.rnaseq import RNASEQ_VERSION


@pytest.fixture
def site(tmp_path):
    """test_cluster_scripts' site: the clone, stubs and job environment."""
    return _cluster_site.__wrapped__(tmp_path)


SITES = ("site_kebnekaise.sh", "site_laptop.sh")
LAUNCHERS = ("fetchngs.sbatch", "rnaseq.sbatch")
UNAME = """#!/bin/sh
if [ "$1" = -m ]; then echo "${STUB_ARCH:-x86_64}"; else exec /usr/bin/uname "$@"; fi
"""
# docker info --format '{{.NCPU}} {{.MemTotal}}', as Docker Desktop answers it, or refused.
DOCKER = """#!/bin/sh
echo "$@" >> "$STUB_LOG/docker.calls"
if [ -n "${STUB_DOCKER_DOWN:-}" ]; then echo "Cannot connect to the Docker daemon" >&2; exit 1; fi
echo "6 12884901888"
"""
NAME = re.compile(r"\b(SITE_[A-Z_]+)=|^(site_\w+)\(\)", re.M)
USED = re.compile(r"\$\{?(SITE_[A-Z_]+)|\b(site_[a-z_]+)\b")
ARCHES = [("arm64", "docker,arm", "docker,arm64"), ("aarch64", "docker,arm", "docker,arm64"),
          ("x86_64", "docker", "docker")]  # fmt: skip


def _defined(name: str) -> set[str]:
    return {a or b for a, b in NAME.findall((BASH / name).read_text())}


def test_sites_define_what_the_launchers_use():
    """Every site sets the same names, and the launchers use no other. Catches a name dropped
    from one site, which fails under set -u on that site alone, and a launcher reading a name no
    site sets."""
    kebnekaise, laptop = (_defined(n) for n in SITES)
    assert kebnekaise == laptop and len(kebnekaise) == 7, (kebnekaise, laptop)
    code = "\n".join(s for n in LAUNCHERS for s in (BASH / n).read_text().splitlines()
                     if not s.lstrip().startswith("#"))  # fmt: skip
    used = {a or b for a, b in USED.findall(code)}
    assert used and used <= kebnekaise | {"SITE_FILE"}, used - kebnekaise


def _laptop(site, arch, **given):
    """The site fixture as a laptop: uname and docker stubbed, no store, a home without spaces."""
    tmp_path, env = site
    _stub(tmp_path / "stubs" / "uname", UNAME)
    _stub(tmp_path / "stubs" / "docker", DOCKER)
    home = tmp_path / "home"
    home.mkdir()
    env = {k: v for k, v in env.items() if k != "CTX_STORE"}
    return tmp_path, {**env, "HOME": str(home), "CTX_SITE": "laptop", "STUB_ARCH": arch,
                      **given}, home  # fmt: skip


@pytest.mark.parametrize(("arch", "fetchngs", "rnaseq"), ARCHES)
def test_laptop_fetchngs(site, arch, fetchngs, rnaseq):
    """A download on a laptop, its run folder in home. Catches the Kebnekaise environment on a
    laptop (it needs a store and sets the Slurm executor), the two arm profiles swapped (fetchngs
    1.12.0 names it arm, rnaseq 3.27.0 arm64), home refused on a laptop, and the fasterq-dump cap
    dropped."""
    del rnaseq
    tmp_path, env, home = _laptop(site, arch)
    run = home / "ctx_runs" / "four"
    run.mkdir(parents=True)
    (run / "ids.txt").write_text("".join(f"{i}\n" for i in _ids("download")))
    done = _run([_spool(tmp_path, "fetchngs.sbatch"), str(run), "ids.txt"], run, env,
                STUB_FIXTURE=str(FETCHNGS / "download"))  # fmt: skip
    assert done.returncode == 0, done.stderr
    argv, nf_env, _ = _recorded(tmp_path)
    assert argv[:4] == ["run", "nf-core/fetchngs", "-r", FETCHNGS_VERSION]
    assert _after(argv, "-profile") == fetchngs
    assert _after(argv, "--max_cpus") == "4" and _after(argv, "--max_memory") == "8.GB"
    assert nf_env["NXF_SYNTAX_PARSER"] == "v1" and "NXF_EXECUTOR" not in nf_env
    assert not (tmp_path / "log" / "pixi.args").exists()
    cores = subprocess.run(["getconf", "_NPROCESSORS_ONLN"], capture_output=True, text=True,
                           check=True).stdout.strip()  # fmt: skip
    config = (run / "nfcore_fetchngs" / "nextflow.config").read_text()
    assert re.search(rf"^\s*queueSize = {cores}$", config, re.M), config
    assert (run / "nfcore_fetchngs" / "results.complete").is_file()


@pytest.mark.parametrize(("arch", "fetchngs", "rnaseq"), ARCHES)
def test_laptop_rnaseq(site, arch, fetchngs, rnaseq):
    """Counting on a laptop. Catches the profiles swapped, the Slurm executor or the Kebnekaise
    queue on a laptop, and tasks let ask for more than Docker has (Nextflow refuses them)."""
    del fetchngs
    tmp_path, env, _ = _laptop(site, arch)
    done, _, launch = _rnaseq((tmp_path, env))
    assert done.returncode == 0, done.stderr
    argv, nf_env, _ = _recorded(tmp_path)
    assert argv[:4] == ["run", "nf-core/rnaseq", "-r", RNASEQ_VERSION]
    assert _after(argv, "-profile") == rnaseq and "NXF_EXECUTOR" not in nf_env
    config = (launch / "site.config").read_text()
    assert "resourceLimits = [cpus: 6, memory: 12288.MB]" in config, config
    assert "queueSize = 50" not in config


def test_laptop_rnaseq_without_docker(site):
    """Docker not running stops the launch before anything is written. Catches the limits read
    after the launch folder is made, or an empty answer taken as limits."""
    tmp_path, env, _ = _laptop(site, "arm64", STUB_DOCKER_DOWN="1")
    done, _, launch = _rnaseq((tmp_path, env))
    assert done.returncode == 1 and "start Docker" in done.stderr, done.stderr
    assert not launch.exists()


def test_unknown_site(site):
    """CTX_SITE names a site file that exists. Catches the default taken for a typing error."""
    tmp_path, env = site
    done, _, _ = _rnaseq((tmp_path, {**env, "CTX_SITE": "lapotp"}))
    assert done.returncode == 1 and "site_lapotp.sh" in done.stderr, done.stderr


def _sums(**names: bytes) -> bytes:
    return "".join(f"{hashlib.md5(data).hexdigest()}  {n}\n" for n, data in names.items()).encode()


REFS = {"genome.fa.gz": b"genome", "annotation.gtf.gz": b"annotation"}
CASES = {
    "absent": (None, 0, "no MD5SUMS beside it"),
    "equal": (_sums(**REFS), 0, "MD5 equal to"),
    "binary_marker": (_sums(**{f"*{k}": v for k, v in REFS.items()}), 0, "MD5 equal to"),
    "names_neither": (_sums(**{"other.fa.gz": b"x", "other.gtf.gz": b"y"}), 1, "no single line"),
    "empty": (b"", 1, "no single line"),
    "listed_twice": (_sums(**REFS) * 2, 1, "no single line"),
    "wrong": (_sums(**{**REFS, "annotation.gtf.gz": b"other"}), 1, "download it again"),
}


@pytest.mark.parametrize("case", CASES)
def test_rnaseq_md5sums(site, case):
    """The references against the publisher's MD5SUMS beside them, before the pipeline runs.
    Catches a list that names neither file passing (an empty md5sum -c passes on macOS), a hash
    that differs let through, a missing list passing unsaid, and a file listed twice."""
    sums, code, said = CASES[case]
    tmp_path, _ = site
    if sums is not None:
        (tmp_path / "submit" / "run" / "refs").mkdir(parents=True)
        (tmp_path / "submit" / "run" / "refs" / "MD5SUMS").write_bytes(sums)
    done, _, _ = _rnaseq(site)
    assert done.returncode == code and said in done.stdout + done.stderr, done
    assert (tmp_path / "log" / "nextflow.argv").exists() == (code == 0)
