"""The four Kebnekaise scripts in scripts/bash, run here behind stub executables.

No job is submitted and nothing reaches the network: pixi, nextflow, apptainer and wget are
stubs on a PATH of the system folders, and every script runs from a copy outside the
repository, as sbatch runs it from its spool folder. What only the cluster shows is G.5's.
Each test names the mutation it catches.
"""

import csv
import itertools
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from ctx.dataset import load
from ctx.gather.fetchngs import FETCHNGS_VERSION
from ctx.gather.rnaseq import RNASEQ_VERSION

ROOT = Path(__file__).resolve().parents[1]
BASH = ROOT / "scripts" / "bash"
SCRIPTS = ("kebnekaise_env.sh", "fetchngs.sbatch", "rnaseq.sbatch",
           "site_kebnekaise.sh", "site_laptop.sh", "fetch_gdc.sbatch", "job_name.sh")  # fmt: skip
JOBS = tuple(name for name in SCRIPTS if name.endswith(".sbatch"))
FETCHNGS = ROOT / "tests" / "fixtures" / "fetchngs"
SYSTEM_PATH = "/usr/bin:/bin:/usr/sbin:/sbin"
# A script students copy must not carry the lab's account, storage path, user or address.
LAB = re.compile(r"hpc2n20\d\d-\d+|/proj/nobackup|/pfs/|rstak|/Users/|[\w.+-]+@[\w-]+\.[\w.]+")
# The one exit the sourced file may hold: it runs only when the file is run as a program.
GUARD = (
    '[[ "${BASH_SOURCE[0]}" != "$0" ]] || '
    '{ echo "kebnekaise_env.sh: source this file, do not run it" >&2; exit 1; }'
)
HOME_CHECK = re.compile(r"^# Home is 25 GB.*?^esac$", re.M | re.S)
# What the end check prints and the Usage tells the student to read before the download.
SIZE_LINE = "bytes of FASTQ at ENA"
CUT = ("sample", "fastq_1", "fastq_2", "strandedness")

PIXI = """#!/bin/sh
# pixi 0.76.2's shell-hook as the environment file meets it: exit 0 and a full hook whether or
# not the environment was installed, PATH by the physical path (as pixi printed it on the
# cluster), a prompt unless --change-ps1 false, and a pixi function.
echo "$@" > "$STUB_LOG/pixi.args"
if [ -n "${STUB_PIXI_FAIL:-}" ]; then echo "stub pixi: refused" >&2; exit 1; fi
if [ -n "${STUB_PIXI_EMPTY:-}" ]; then exit 0; fi
manifest=""; ps1=true; prev=""
for a in "$@"; do
    case "$prev" in --manifest-path) manifest=$a ;; --change-ps1) ps1=$a ;; esac
    prev=$a
done
[ -f "$manifest" ] || exit 1
root=$(cd "$(dirname "$manifest")" && pwd -P)
echo "export PATH=\\"$root/.pixi/envs/fetch/bin:\\$PATH\\""
echo "export CONDA_PREFIX=\\"$root/.pixi/envs/fetch\\""
[ "$ps1" = false ] || echo 'PS1="(fetch) ${PS1-}"'
echo 'pixi() { :; }'
"""
NEXTFLOW = """#!/usr/bin/env bash
# nf-core/fetchngs or nf-core/rnaseq, as far as the scripts can tell: argv, environment, folder
# and any marker left from an earlier launch recorded, then an outdir published.
set -eu
printf '%s\\n' "$@" > "$STUB_LOG/nextflow.argv"
env > "$STUB_LOG/nextflow.env"
pwd -P > "$STUB_LOG/nextflow.cwd"
if [ -e results.complete ]; then echo yes; else echo no; fi > "$STUB_LOG/marker_at_start"
prev=""; out=""; params=""
for a in "$@"; do
    case "$prev" in --outdir) out=$a ;; -params-file) params=$a ;; esac
    prev=$a
done
if [ -n "$out" ]; then
    cp -R "$STUB_FIXTURE/." "$out/"
    if [ -n "${STUB_DROP:-}" ]; then rm "$out/$STUB_DROP"; fi
else
    out=$(sed -n 's/^outdir: //p' "$params")
    mkdir -p "$out/salmon"
    if [ -z "${STUB_NO_MATRIX:-}" ]; then
        printf 'gene_id\\tgene_name\\ts1\\ts2\\ng1\\tA\\t1\\t2\\ng2\\tB\\t3\\t4\\n' \\
            > "$out/salmon/salmon.merged.gene_counts_length_scaled.tsv"
        if grep -q '^aligner: star_salmon$' "$params"; then
            mkdir -p "$out/star_salmon"
            cp "$out/salmon/salmon.merged.gene_counts_length_scaled.tsv" "$out/star_salmon/"
        fi
    fi
    if [ -n "${STUB_LINK:-}" ]; then ln -s /nonexistent "$out/salmon/report.json"; fi
fi
"""
LOGIN_NEXTFLOW = "#!/bin/sh\necho 'the login Nextflow ran' >&2\nexit 99\n"
APPTAINER = """#!/usr/bin/env bash
# apptainer exec [options] IMAGE command: the command on this machine.
echo "$@" >> "$STUB_LOG/apptainer.calls"
shift
while [[ $# -gt 0 && "$1" != *.sif ]]; do shift; done
shift
exec "$@"
"""
# The environment's own Python, which the rnaseq sheet cut must run on (repo L5).
PYTHON3 = '#!/bin/sh\necho called >> "$STUB_LOG/python3.calls"\nexec "{python}" "$@"\n'
WGET = """#!/bin/sh
# wget -nv -O FILE URL: the file written, nothing fetched, or half of it and a failure.
echo "$@" >> "$STUB_LOG/wget.calls"
if [ -n "${STUB_WGET_FAIL:-}" ]; then echo half > "$3"; exit 4; fi
echo "stub image" > "$3"
"""


def _stub(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    path.chmod(0o755)


def _run(args, cwd, env, **given):
    return subprocess.run(
        ["bash", *map(str, args)], cwd=cwd, env={**env, **given}, capture_output=True,
        text=True, timeout=120, check=False,
    )  # fmt: skip


@pytest.fixture
def site(tmp_path):
    """A clone reached through a symlink, a store, a home, stubs, and the environment a job sees.

    The clone, its link, the store and the home have spaces in their names, so an unquoted path
    in a script splits. The submit folder has none: fetchngs refuses a spaced run folder.
    """
    repo = tmp_path / "my clone"
    (repo / "scripts" / "bash").mkdir(parents=True)
    for name in SCRIPTS:
        shutil.copy2(BASH / name, repo / "scripts" / "bash" / name)
    (repo / "pixi.toml").write_text("")
    fetch = repo / ".pixi" / "envs" / "fetch" / "bin"
    _stub(fetch / "nextflow", NEXTFLOW)
    log = tmp_path / "log"
    log.mkdir()
    _stub(fetch / "python3", PYTHON3.format(python=sys.executable))
    (tmp_path / "course repo").symlink_to(repo)
    stubs = tmp_path / "stubs"
    _stub(stubs / "pixi", PIXI)
    _stub(stubs / "nextflow", LOGIN_NEXTFLOW)
    _stub(stubs / "apptainer", APPTAINER)
    _stub(stubs / "wget", WGET)
    for folder in ("my home", "my store", "submit", "spool"):
        (tmp_path / folder).mkdir()
    # A home reached through a symlink, which a check on the strings misses.
    (tmp_path / "submit link").symlink_to(tmp_path / "submit")
    env = {
        "HOME": str(tmp_path / "my home"), "PATH": f"{stubs}:{SYSTEM_PATH}",
        "CTX_REPO": str(tmp_path / "course repo"), "CTX_STORE": str(tmp_path / "my store"),
        "PIXI": str(stubs / "pixi"), "STUB_LOG": str(log),
        "ENV_SH": str(repo / "scripts" / "bash" / "kebnekaise_env.sh"),
    }  # fmt: skip
    return tmp_path, env


def _spool(tmp_path: Path, name: str) -> Path:
    """The script as sbatch runs it: a copy with no sibling files."""
    copy = tmp_path / "spool" / "slurm_script"
    shutil.copy2(BASH / name, copy)
    return copy


def _recorded(tmp_path: Path) -> tuple[list[str], dict[str, str], str]:
    log = tmp_path / "log"
    argv = (log / "nextflow.argv").read_text().splitlines()
    env = dict(line.split("=", 1) for line in (log / "nextflow.env").read_text().splitlines()
               if "=" in line)  # fmt: skip
    return argv, env, (log / "nextflow.cwd").read_text().strip()


def _after(argv: list[str], flag: str) -> str:
    return argv[argv.index(flag) + 1]


def _ids(mode: str) -> list[str]:
    return (FETCHNGS / mode / "ids.csv").read_text().split()


def _usage(name: str) -> list[str]:
    """The Usage block's lines, the comment marker removed."""
    lines = (BASH / name).read_text().splitlines()
    after = lines[lines.index("# Usage:") + 1 :]
    block = itertools.takewhile(lambda s: s.startswith("#   "), after)
    return [s[4:] for s in block]


def _queue_size(config: str) -> bool:
    return bool(re.search(r"^\s*queueSize\s*=\s*50$", config, re.M))


@pytest.mark.parametrize("name", SCRIPTS)
def test_shell_rules(name):
    """Catches a dropped fi, a directive after set -e, an account line (spaced or attached), a
    mail line, a lab path or address, an em dash or arrow (any byte outside ASCII), a shebang on
    the sourced file, a second exit in it or its guard dropped, and an exempt reason naming a ctx
    option the CLI does not have."""
    path = BASH / name
    assert path.read_bytes().isascii(), f"{name}: not ASCII"
    text = path.read_text()
    assert subprocess.run(["bash", "-n", str(path)], check=False).returncode == 0
    lines = text.splitlines()
    first = next(i for i, s in enumerate(lines) if s.strip() and not s.lstrip().startswith("#"))
    directives = [i for i, s in enumerate(lines) if s.startswith("#SBATCH")]
    assert bool(directives) == (name in JOBS) == text.startswith("#!")
    assert all(i < first for i in directives), "Slurm stops reading #SBATCH at the first command"
    assert not [s for s in lines if re.match(r"#SBATCH\s*(-A|--account|--mail)", s)]
    assert not LAB.findall(text), LAB.findall(text)
    (exempt,) = re.findall(r"^# exempt: (.+)$", text, re.M)
    cli = (ROOT / "src" / "ctx" / "cli.py").read_text()
    assert all(f'"{flag}"' in cli for flag in re.findall(r"--[\w-]+", exempt)), exempt
    if name == "kebnekaise_env.sh":
        code = [s for s in lines if s.strip() and not s.lstrip().startswith("#")]
        assert code[0] == GUARD, "without it, `bash kebnekaise_env.sh` runs past every refusal"
        exits = [s for s in code if re.search(r"\bexit\b", s)]
        assert exits == [GUARD], "exit closes a login shell"


# What each head job asks Slurm for: one task, its cores, its memory (the head job's 4 GB heap
# under it) and a walltime its pipeline fits in. rnaseq's 12 GB is Rytis's P4.5 answer: its head
# peaked at 6.16, 6.61 and 7.49 GiB of 8 on three breast studies.
RESOURCES = {
    "fetchngs.sbatch": {"-n": "1", "-c": "1", "--mem": "8G", "-t": "1-00:00:00"},
    "rnaseq.sbatch": {"-n": "1", "-c": "1", "--mem": "12G", "-t": "1-00:00:00"},
    # R1: one core for curl and the end check; six hours, 3.6 times the slowest laptop fetch,
    # 5,965 s for 2.65 GB on 2026-10-18 (GDC served 665 s for 2.6 GB on 2026-10-01).
    "fetch_gdc.sbatch": {"-n": "1", "-c": "1", "--mem": "2G", "-t": "0-06:00:00"},
}


@pytest.mark.parametrize("name", JOBS)
def test_job_resources(name):
    """Catches a walltime cut short, which ends every student's run mid-pipeline, and a core,
    memory or walltime request inflated for a whole class."""
    text = (BASH / name).read_text()
    assert dict(re.findall(r"^#SBATCH (-[nct]|--mem)[ =](\S+)$", text, re.M)) == RESOURCES[name]


def test_home_check_single_source():
    """The run-folder home check is one block, the same text in all three job scripts. Catches
    one copy edited alone."""
    blocks = [HOME_CHECK.findall((BASH / name).read_text()) for name in JOBS]
    assert all(len(block) == 1 for block in blocks) and len({b[0] for b in blocks}) == 1, blocks


def test_env_usage_pins_pixi():
    """The install line pins the pixi CI tests with and leaves ~/.bashrc alone, and the install
    caches on the store. Catches another pixi, PIXI_NO_PATH_UPDATE dropped, and PIXI_CACHE_DIR
    dropped (G.0: pixi's home cache sits on a network filesystem)."""
    ci = (ROOT / ".github" / "workflows" / "ci.yml").read_text()
    (pinned,) = re.findall(r'^\s*PIXI_VERSION:\s*"(v[\d.]+)"', ci, re.M)
    usage = _usage("kebnekaise_env.sh")
    (install,) = [s for s in usage if "pixi.sh/install.sh" in s]
    assert install.endswith(f"| PIXI_VERSION={pinned} PIXI_NO_PATH_UPDATE=1 bash"), install
    (frozen,) = [s for s in usage if " install --manifest-path " in s]
    assert frozen.startswith('PIXI_CACHE_DIR="$CTX_STORE/'), frozen
    assert frozen.endswith(" -e fetch --frozen"), frozen


PROBE = """
PS1='student$ '
before=$PATH
source "$ENV_SH"
echo "status=$?"
echo "nextflow=$(command -v nextflow)"
echo "options=$-"
echo "pipefail=$(set -o | awk '$1 == "pipefail" {print $2}')"
echo "path_kept=$([[ "$PATH" == "$before" ]] && echo yes || echo no)"
echo "PS1=$PS1"
echo "helpers=$( { compgen -A variable _ctx_; compgen -A function _ctx_; } | tr '\\n' ' ')"
for v in NXF_VER NXF_ASSETS NXF_JAVA_HOME NXF_HOME NXF_EXECUTOR NXF_OPTS NXF_ANSI_LOG \\
    SBATCH_ACCOUNT NXF_SINGULARITY_CACHEDIR NXF_APPTAINER_CACHEDIR APPTAINER_CACHEDIR \\
    APPTAINER_TMPDIR; do
    echo "$v=${!v-unset}"
done
echo alive
"""


def _probe(env: dict[str, str], **given) -> dict[str, str]:
    """Source the file by its own path, never through CTX_REPO, which a case may empty."""
    done = subprocess.run(
        ["bash", "-c", PROBE], env={**env, **given}, capture_output=True, text=True,
        timeout=60, check=False,
    )  # fmt: skip
    assert done.returncode == 0, done.stderr
    return dict(line.split("=", 1) for line in done.stdout.splitlines() if "=" in line) | {
        "alive": str("alive" in done.stdout.splitlines()),
        "stderr": done.stderr,
    }


@pytest.mark.parametrize("pixi_from", ["PIXI", "PIXI_HOME", "HOME"])
def test_env_activation(site, pixi_from):
    """Catches the nextflow path compared as a string (a symlinked clone fails it), NXF_VER
    inherited, set -euo pipefail left on in the caller's shell, the prompt changed, PIXI_HOME
    ignored (G.0: the installer honours it), NXF_OPTS or NXF_ANSI_LOG dropped, and Apptainer's
    scratch on the store (Lustre fails with unpriv.lsetxattr)."""
    tmp_path, env = site
    env = dict(env)
    if pixi_from != "PIXI":
        del env["PIXI"]
        found = tmp_path / "pixi home" if pixi_from == "PIXI_HOME" else Path(env["HOME"]) / ".pixi"
        _stub(found / "bin" / "pixi", PIXI)
        if pixi_from == "PIXI_HOME":
            env["PIXI_HOME"] = str(found)
    node = tmp_path / "node tmp"
    got = _probe(env, NXF_VER="23.04.0", NXF_ASSETS="/elsewhere", NXF_JAVA_HOME="/jdk",
                 SLURM_JOB_ACCOUNT="compute_project", SNIC_TMP=str(node))  # fmt: skip
    assert got["status"] == "0", got["stderr"]
    store = Path(env["CTX_STORE"]).resolve()
    fetch = (tmp_path / "my clone" / ".pixi" / "envs" / "fetch" / "bin").resolve()
    assert got["nextflow"] == str(fetch / "nextflow")
    assert not {"e", "u"} & set(got["options"]) and got["pipefail"] == "off"
    assert got["PS1"] == "student$ " and got["helpers"].strip() == ""
    assert got["NXF_VER"] == got["NXF_ASSETS"] == got["NXF_JAVA_HOME"] == "unset"
    assert got["NXF_HOME"] == str(store / "nextflow") and got["NXF_EXECUTOR"] == "slurm"
    assert got["NXF_ANSI_LOG"] == "false" and got["SBATCH_ACCOUNT"] == "compute_project"
    # The head job's heap under the memory both head jobs ask Slurm for.
    heap = re.fullmatch(r"-Xms\d+M -Xmx(\d+)G", got["NXF_OPTS"])
    head = "".join((BASH / n).read_text() for n in ("fetchngs.sbatch", "rnaseq.sbatch"))
    mems = [int(m) for m in re.findall(r"^#SBATCH --mem=(\d+)G$", head, re.M)]
    assert heap and len(mems) == 2 and int(heap[1]) < min(mems), (got["NXF_OPTS"], mems)
    assert got["NXF_SINGULARITY_CACHEDIR"] == got["NXF_APPTAINER_CACHEDIR"] == str(store / "images")
    assert not got["APPTAINER_CACHEDIR"].startswith(got["NXF_SINGULARITY_CACHEDIR"])
    assert got["APPTAINER_TMPDIR"] == str(node), got["APPTAINER_TMPDIR"]
    assert all((store / d).is_dir() for d in ("nextflow", "images", "apptainer_cache"))
    args = (tmp_path / "log" / "pixi.args").read_text().split()
    assert _after(args, "-e") == "fetch" and "--as-is" in args


@pytest.mark.parametrize(
    "case", ["no_repo", "no_store", "store_in_home", "store_linked_into_home", "no_environment",
             "hook_fails", "hook_leaves_login_nextflow"],
)  # fmt: skip
def test_env_refuses(site, case):
    """Each refusal returns 1 with its own message, the calling shell goes on, and its PATH,
    prompt and variables are as they were. Catches `exit` for `return`, a check deleted (another
    one's message then shows), the hook evaluated for an environment never installed, the hook
    evaluated without the nextflow check (student L1), the home check made on the strings (a
    store linked into home passes it), and helpers left behind."""
    tmp_path, env = site
    home = Path(env["HOME"])
    (home / "store").mkdir()
    (tmp_path / "linked store").symlink_to(home / "store")
    store = {
        "no_store": f"{tmp_path}/missing", "store_in_home": f"{home}/store",
        "store_linked_into_home": f"{tmp_path}/linked store",
    }.get(case, env["CTX_STORE"])  # fmt: skip
    given, said = {
        "no_repo": ({"CTX_REPO": ""}, "export CTX_REPO and CTX_STORE first"),
        "no_store": ({"CTX_STORE": store}, f"CTX_STORE={store} is not a folder"),
        "store_in_home": ({"CTX_STORE": store}, f"CTX_STORE={store} is in your home"),
        "store_linked_into_home": ({"CTX_STORE": store}, f"CTX_STORE={store} is in your home"),
        "no_environment": ({}, f"no nextflow in {env['CTX_REPO']}/.pixi/envs/fetch/bin;"),
        "hook_fails": ({"STUB_PIXI_FAIL": "1"}, f"'{env['PIXI']} shell-hook' failed"),
        "hook_leaves_login_nextflow": ({"STUB_PIXI_EMPTY": "1"},
                                       f"nextflow is {tmp_path}/stubs/nextflow, not the one in"),
    }[case]  # fmt: skip
    if case == "no_environment":
        shutil.rmtree(tmp_path / "my clone" / ".pixi")
    got = _probe(env, **given)
    assert got["status"] == "1" and got["alive"] == "True", got
    mine = [s for s in got["stderr"].splitlines() if s.startswith("kebnekaise_env.sh: ")]
    assert len(mine) == 1 and mine[0].startswith(f"kebnekaise_env.sh: {said}"), got["stderr"]
    assert got["path_kept"] == "yes" and got["PS1"] == "student$ ", got
    assert got["helpers"].strip() == "" and got["NXF_HOME"] == "unset", got


def test_env_run_as_a_program(site):
    """`bash kebnekaise_env.sh` exits 1 before any check. Catches the guard dropped: every
    return is then a no-op, and a store refused as lying in home gets the cache folders."""
    tmp_path, env = site
    store = Path(env["HOME"]) / "store"
    store.mkdir()
    done = _run([env["ENV_SH"]], tmp_path, env, CTX_STORE=str(store))
    assert done.returncode == 1, done
    assert done.stderr == "kebnekaise_env.sh: source this file, do not run it\n"
    assert not list(store.iterdir())


def _write_list(path: Path, ids: list[str], bad: str | None = None) -> None:
    """A list as a spreadsheet saves it: a byte-order mark, CRLF, a blank, a comment, no final
    newline."""
    lines = [ids[0], "", "# chosen runs", *([bad] if bad else []), f"  {ids[1]} \t", *ids[2:]]
    path.write_bytes(b"\xef\xbb\xbf" + "\r\n".join(lines).encode())


def test_fetchngs_list(site):
    """The dry run keeps exactly the accessions. Catches a byte-order mark or CR kept."""
    tmp_path, env = site
    ids = _ids("download")
    _write_list(tmp_path / "submit" / "ids.txt", ids)
    done = _run([BASH / "fetchngs.sbatch", "run", "ids.txt"], tmp_path / "submit", env, DRY_RUN="1")
    assert done.returncode == 0, done.stderr
    written = tmp_path / "submit" / "run" / "nfcore_fetchngs" / "ids.csv"
    assert written.read_text() == "".join(f"{i}\n" for i in ids)
    assert f"{len(ids)} accessions" in done.stdout
    assert not (tmp_path / "submit" / "run" / "nfcore_fetchngs" / "results").exists()
    assert not (tmp_path / "log" / "nextflow.argv").exists()


@pytest.mark.parametrize(
    ("bad", "said"),
    [("SAME123", "not an SRA"), ("{0} {1}", "one column"), ("{0},{1}", "one column"),
     ("E-MTAB-{n}", "give its ENA study")],
)  # fmt: skip
def test_fetchngs_list_refuses(site, bad, said):
    """A line outside fetchngs's pattern exits 1 naming its line (line 4: after the first id,
    the blank and the comment). Catches the pattern check dropped and the E-MTAB hint lost."""
    tmp_path, env = site
    ids = _ids("download")
    # Built at run time: an ArrayExpress id written here would trip test_no_cohort_identifiers.
    _write_list(tmp_path / "submit" / "ids.txt", ids, bad.format(*ids, n=len(ids) * 1000))
    done = _run([BASH / "fetchngs.sbatch", "run", "ids.txt"], tmp_path / "submit", env, DRY_RUN="1")
    assert done.returncode == 1 and "line 4:" in done.stderr and said in done.stderr, done.stderr
    assert not (tmp_path / "submit" / "run" / "nfcore_fetchngs" / "ids.csv").exists()


@pytest.mark.parametrize(
    ("case", "said"),
    [("space_in_run", "has a space"), ("relative_run_in_spaced_folder", "has a space"),
     ("results_exist", "exists; move it aside"), ("only_comments", "no accession left"),
     ("run_in_home", "is in your home folder"), ("run_linked_into_home", "is in your home folder")],
)  # fmt: skip
def test_fetchngs_refuses_the_run(site, case, said):
    """Catches a relaunch into an existing results/ (a second params file), the space check made
    before the run folder is absolute (a relative run under a spaced folder then reaches
    fetchngs, whose --input pattern refuses it after the queue), and the home check dropped or
    made on the strings (a home reached through a symlink)."""
    tmp_path, env = site
    submit = tmp_path / ("my course" if case == "relative_run_in_spaced_folder" else "submit")
    submit.mkdir(exist_ok=True)
    _write_list(submit / "ids.txt", _ids("download"))
    run = "my run" if case == "space_in_run" else "run"
    given = {"DRY_RUN": "1"}
    if case == "results_exist":
        (submit / "run" / "nfcore_fetchngs" / "results").mkdir(parents=True)
    if case == "only_comments":
        (submit / "ids.txt").write_text("# none yet\n\n")
    if case in ("run_in_home", "run_linked_into_home"):
        given["HOME"] = str(tmp_path / ("submit" if case == "run_in_home" else "submit link"))
    done = _run([BASH / "fetchngs.sbatch", run, "ids.txt"], submit, env, **given)
    assert done.returncode == 1 and done.stderr.startswith("fetchngs.sbatch: "), done.stderr
    assert said in done.stderr, done.stderr
    if case not in ("results_exist", "only_comments"):
        assert not (submit / run).exists(), "refused after writing"


def test_fetchngs_relaunch_names_a_free_folder(site):
    """Three failed launches, each moved aside with the mv the dry run prints, leave three
    sibling folders. Catches a fixed results_try1, which nests the second and stops the third."""
    tmp_path, env = site
    submit = tmp_path / "submit"
    _write_list(submit / "ids.txt", _ids("download"))
    launch = submit / "run" / "nfcore_fetchngs"
    for attempt in ("1", "2", "3"):
        (launch / "results").mkdir(parents=True)
        (launch / "results" / "attempt").write_text(attempt)
        done = _run([BASH / "fetchngs.sbatch", "run", "ids.txt"], submit, env, DRY_RUN="1")
        (move,) = re.findall(r"\((mv \S+ \S+)\)", done.stderr)
        assert subprocess.run(move.split(), check=False).returncode == 0, move
    tries = sorted(launch.iterdir())
    assert [p.name for p in tries] == ["results_try1", "results_try2", "results_try3"]
    assert [(p / "attempt").read_text() for p in tries] == ["1", "2", "3"]


def test_fetchngs_usage_waits_for_the_size():
    """Between the metadata job and the download, the Usage reads the size the metadata job
    printed, from that job's own log. Catches the two sbatch lines back to back, which start the
    download unread and put two of one student's runs in the queue at once."""
    text = (BASH / "fetchngs.sbatch").read_text()
    usage = _usage("fetchngs.sbatch")
    meta = next(i for i, s in enumerate(usage) if s.startswith("SKIP_DOWNLOAD=1 sbatch "))
    download = next(i for i, s in enumerate(usage) if i > meta and s.startswith("sbatch "))
    (job,) = re.findall(r"^#SBATCH -J (\S+)$", text, re.M)
    (log,) = re.findall(r"^#SBATCH -o (\S+)$", text, re.M)
    log = log.replace("%x", job).replace("%j", "<job id>")
    between = [s for s in usage[meta + 1 : download] if not s.startswith("#")]
    assert between == [f"grep '{SIZE_LINE}' {log}"], between


@pytest.mark.parametrize(("mode", "forks"), [("metadata_only", None), ("download", "3")])
def test_fetchngs_stubbed(site, mode, forks):
    """From a spool copy, with a relative run folder. Catches parser v2, the environment found
    via $0, a relative outdir (student L2), another profile, FASTQs published as symlinks into
    work/ (D-006), MAX_FORKS or queueSize not taken, and an earlier launch's marker still in
    place while the relaunch runs (F5)."""
    tmp_path, env = site
    submit = tmp_path / "submit"
    shutil.copy2(FETCHNGS / mode / "ids.csv", submit / "ids.txt")
    meta = mode == "metadata_only"
    launch = submit.resolve() / "run" / ("nfcore_fetchngs_meta" if meta else "nfcore_fetchngs")
    launch.mkdir(parents=True)
    (launch / "results.complete").write_text("an earlier launch\n")
    done = _run([_spool(tmp_path, "fetchngs.sbatch"), "run", "ids.txt"], submit, env,
                SKIP_DOWNLOAD=str(int(meta)), STUB_FIXTURE=str(FETCHNGS / mode),
                SLURM_JOB_ACCOUNT="compute_project",
                **({"MAX_FORKS": forks} if forks else {}))  # fmt: skip
    assert done.returncode == 0, done.stderr
    argv, nf_env, cwd = _recorded(tmp_path)
    assert argv[:4] == ["run", "nf-core/fetchngs", "-r", FETCHNGS_VERSION]
    assert _after(argv, "-profile") == "singularity"
    assert _after(argv, "--publish_dir_mode") == "link"
    assert _after(argv, "--monochrome_logs") == "true"
    assert _after(argv, "--nf_core_pipeline") == "rnaseq" and "-resume" in argv
    assert ("--skip_fastq_download" in argv) == meta
    assert _after(argv, "--outdir") == str(launch / "results")
    assert _after(argv, "--input") == str(launch / "results" / "ids.csv")
    assert cwd == str(launch)
    assert nf_env["NXF_SYNTAX_PARSER"] == "v1" and nf_env["SBATCH_ACCOUNT"] == "compute_project"
    assert "NXF_VER" not in nf_env
    config = (launch / "nextflow.config").read_text()
    assert _queue_size(config) and re.search(rf"^\s*maxForks\s*= {forks or 5}$", config, re.M)
    assert (tmp_path / "log" / "marker_at_start").read_text() == "no\n"
    assert (launch / "results.complete").read_text() != "an earlier launch\n"
    assert SIZE_LINE in done.stdout


@pytest.mark.parametrize(
    ("mode", "drop"),
    [("metadata_only", "metadata"), ("download", "single_end"), ("download", "mate_2"),
     ("download", "samplesheet")],
)  # fmt: skip
def test_fetchngs_end_check(site, mode, drop):
    """An outdir missing one metadata file, a single-end FASTQ, a mate 2 or the samplesheet exits
    1 and leaves no marker. Catches a marker on a partial download (F9), mate 2 not checked, an
    empty sheet let through, and an earlier launch's marker kept."""
    tmp_path, env = site
    submit = tmp_path / "submit"
    shutil.copy2(FETCHNGS / mode / "ids.csv", submit / "ids.txt")
    meta = mode == "metadata_only"
    fastq = sorted(p.name for p in (FETCHNGS / "download" / "fastq").glob("*.fastq.gz"))
    single = [n for n in fastq if not re.search(r"_[12]\.fastq\.gz$", n)]
    path = {
        "metadata": f"metadata/{_ids(mode)[0]}.runinfo_ftp.tsv",
        "single_end": f"fastq/{single[0]}",
        "mate_2": "fastq/" + next(n for n in fastq if n.endswith("_2.fastq.gz")),
        "samplesheet": "samplesheet/samplesheet.csv",
    }[drop]
    said = "no rows in samplesheet/samplesheet.csv" if drop == "samplesheet" else f"no {path}"
    launch = submit / "run" / ("nfcore_fetchngs_meta" if meta else "nfcore_fetchngs")
    launch.mkdir(parents=True)
    (launch / "results.complete").write_text("an earlier launch\n")
    done = _run([_spool(tmp_path, "fetchngs.sbatch"), "run", "ids.txt"], submit, env,
                SKIP_DOWNLOAD=str(int(meta)), STUB_FIXTURE=str(FETCHNGS / mode),
                STUB_DROP=path)  # fmt: skip
    assert done.returncode == 1 and said in done.stderr, done.stderr
    assert not (launch / "results.complete").exists()


def _download(run: Path, reverse: bool = False) -> list[dict[str, str]]:
    """The download fixture as a finished fetchngs.sbatch run: its outdir and the marker."""
    results = run / "nfcore_fetchngs" / "results"
    shutil.copytree(FETCHNGS / "download", results)
    (run / "nfcore_fetchngs" / "results.complete").write_text("done\n")
    sheet = results / "samplesheet" / "samplesheet.csv"
    with sheet.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    if reverse:
        with sheet.open("w", newline="") as fh:
            out = csv.DictWriter(fh, fieldnames=list(rows[0])[::-1], quoting=csv.QUOTE_ALL)
            out.writeheader()
            out.writerows(rows)
    return rows


REFS = {"genome.fa.gz": b"genome", "annotation.gtf.gz": b"annotation"}


def _rnaseq(site, gencode="true", marker=True, refs=REFS, run="run", **given):
    """rnaseq.sbatch from a spool copy on the download fixture, the references under run/refs.

    A reference whose bytes are None is named but not written. Environment values in ``given``
    may name ``{tmp}``, the test's folder.
    """
    tmp_path, env = site
    submit = tmp_path / "submit"
    rows = _download(submit / run, reverse=True)
    if not marker:
        (submit / run / "nfcore_fetchngs" / "results.complete").unlink()
    for name, data in refs.items():
        path = submit / run / "refs" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        if data is not None:
            path.write_bytes(data)
    args = [run, *(f"{run}/refs/{name}" for name in refs), gencode]
    done = _run([_spool(tmp_path, "rnaseq.sbatch"), *args], submit, env,
                **{k: v.format(tmp=tmp_path) for k, v in given.items()})  # fmt: skip
    return done, rows, submit.resolve() / run / "nfcore_rnaseq"


@pytest.mark.parametrize("gencode", ["true", "false"])
def test_rnaseq_stubbed(site, gencode):
    """Catches a cut by position, the cut on the system Python (repo L5), parser v1, another
    profile, queueSize dropped, and any params key changed or dropped: publish_dir_mode, gencode
    not taken from the argument, Kebnekaise's STAR route (D-113) or its Salmon beside it lost,
    and a relative GTF."""
    done, rows, launch = _rnaseq(site, gencode=gencode)
    tmp_path = site[0]
    assert done.returncode == 0, done.stderr
    with (launch / "samplesheet.csv").open(newline="") as fh:
        cut = list(csv.DictReader(fh))
    assert list(cut[0]) == list(CUT)
    assert cut == [{k: r[k] for k in CUT} for r in rows]
    assert (tmp_path / "log" / "python3.calls").read_text().split() == ["called"]
    params = dict(line.split(": ", 1) for line in (launch / "params.yml").read_text().splitlines())
    refs = launch.parent / "refs"
    assert params == {
        "input": str(launch / "samplesheet.csv"), "outdir": str(launch / "results"),
        "fasta": str(refs / "genome.fa.gz"), "gtf": str(refs / "annotation.gtf.gz"),
        "gencode": gencode, "pseudo_aligner": "salmon", "aligner": "star_salmon",
        "skip_preseq": "false", "save_reference": "true", "publish_dir_mode": "copy",
    }  # fmt: skip
    argv, nf_env, cwd = _recorded(tmp_path)
    assert argv[:4] == ["run", "nf-core/rnaseq", "-r", RNASEQ_VERSION]
    assert _after(argv, "-profile") == "singularity"
    assert _after(argv, "-params-file") == str(launch / "params.yml") and "-resume" in argv
    config = Path(_after(argv, "-c")).read_text()
    assert _queue_size(config)
    assert "withName: '.*:FASTQ_FASTQC_UMITOOLS_TRIMGALORE:TRIMGALORE'" in config
    # Half of nf-core 3.27.0's task limits on Kebnekaise (Rytis 2026-10-23); not the head job's.
    assert "withLabel: process_high { time = { 8.h * task.attempt } }" in config
    assert "    time = { 2.h * task.attempt }" in config
    # Measured requests, memory and time 25 % more per retry, cores never (Rytis 2026-10-23).
    assert "withName: '.*:TRIMGALORE' { cpus = 3; memory = { 4148.MB * Math.pow(1.25" in config
    assert nf_env["NXF_SYNTAX_PARSER"] == "v2" and cwd == str(launch)
    assert "count matrix: 2 genes, 2 samples" in done.stdout
    assert "alignment count matrix: 2 genes, 2 samples" in done.stdout


@pytest.mark.parametrize(
    ("given", "said"),
    [({"marker": False}, "results.complete"),
     ({"refs": {**REFS, "genome.fa.gz": None}}, "no genome FASTA"),
     ({"refs": {**REFS, "annotation.gtf.gz": b""}}, "no GTF"),
     ({"gencode": "yes"}, "GENCODE is true or false"),
     ({"run": "my run"}, "has a space"),
     ({"refs": {"my refs/g.fa.gz": b"g", "a.gtf.gz": b"a"}}, "has a space"),
     ({"refs": {"g.fa.bgz": b"g", "a.gtf.gz": b"a"}}, "genome FASTA only if its name matches"),
     ({"refs": {"g.fa.gz": b"g", "a.gff3": b"a"}}, "GTF only if its name matches"),
     ({"HOME": "{tmp}/submit"}, "is in your home folder"),
     ({"HOME": "{tmp}/submit link"}, "is in your home folder"),
     ({"STUB_LINK": "1"}, "links into work/"),
     ({"STUB_NO_MATRIX": "1"}, "no count matrix")],
    ids=["no_marker", "no_fasta", "empty_gtf", "gencode_yes", "space_in_run", "space_in_refs",
         "fasta_name", "gtf_name", "run_in_home", "run_linked_into_home", "link_into_work",
         "no_matrix"],
)  # fmt: skip
def test_rnaseq_refuses(site, given, said):
    """Exits 1, and before the launch folder exists unless the pipeline has run. Catches rnaseq
    started on a partial download (student L8), an unchecked reference (student L2), a guessed
    GENCODE (student L11), a path nf-core/rnaseq's schema rejects after the queue, the home check
    dropped or made on the strings, a link into work/ let through, and no count matrix."""
    done, _, launch = _rnaseq(site, **given)
    assert done.returncode == 1 and said in done.stderr, done.stderr
    ran = "STUB_LINK" in given or "STUB_NO_MATRIX" in given
    assert launch.exists() == ran


def test_jobs_are_named_for_their_project_and_task(site):
    """A running head job renames itself <ANALYSIS>_<PROJECT>_nf-head_<VERSION>, and the config
    names each task the same way, both from job_name.sh. Catches the rename dropped, the run
    folder not used as the project, and the task line missing from the executor block."""
    tmp_path, _ = site
    _stub(tmp_path / "stubs" / "scontrol", '#!/bin/sh\necho "$@" >> "$STUB_LOG/scontrol.calls"\n')
    datasets = tmp_path / "my clone" / "conf" / "datasets"
    datasets.mkdir(parents=True)
    (datasets / "BRCA-immune-cells-macrophages.toml").write_text(
        'label = "BRCA immune cells: macrophages (TAM)"\n'
    )
    done, _, launch = _rnaseq(site, run="BRCA-immune-cells-macrophages", SLURM_JOB_ID="123")
    assert done.returncode == 0, done.stderr
    calls = (tmp_path / "log" / "scontrol.calls").read_text().split()
    want = "RNAseq-Tutorial_BRCA-immune-cells-macrophages-TAM_nf-head_NA"
    assert calls == ["update", "JobId=123", f"JobName={want}"]
    config = (launch / "site.config").read_text()
    assert 'jobName = { "RNAseq-Tutorial_BRCA-immune-cells-macrophages-TAM_nf-" + task.n' in config
    assert (".replaceAll(/[()]/, '').replaceAll(/[^A-Za-z0-9.-]/, '-').replaceAll(/-+/, '-')"
            in config)  # fmt: skip
    call = f'source "{BASH / "job_name.sh"}"; job_name "a_b, c: d" "nf-STAR_ALIGN (S1)"'
    named = subprocess.run(["bash", "-c", call], capture_output=True, text=True, check=True).stdout
    assert named == "RNAseq-Tutorial_a-b-c-d_nf-STAR-ALIGN-S1_NA" and len(named.split("_")) == 4


def _job(call: str, repo: Path = ROOT) -> str:
    script = f'source "{BASH / "job_name.sh"}"; {call}'
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True, check=True,
                          env={**os.environ, "CTX_REPO": str(repo)}).stdout  # fmt: skip


def test_a_dataset_name_is_its_label_by_the_job_name_rule():
    """Each dataset named from its label (D-114) carries the name job_field makes of that label,
    as its slug and its file's stem, and a run folder of that name is that project. A folder
    with no dataset file is its own name. Catches: a slug typed apart from its label, the rule
    keeping parentheses, and job_project losing the new-name file or the fallback."""
    named = 0
    for path in sorted((ROOT / "conf" / "datasets").glob("*.toml")):
        dataset = load(path)
        if dataset.slug.islower() and "-" not in dataset.slug:
            continue  # a slug from before D-114, not made from its label
        named += 1
        assert dataset.slug == path.stem == _job(f'job_field "{dataset.name}"'), path.name
        assert _job(f'job_field "$(job_project "{path.stem}")"') == path.stem
    assert named == 7  # noqa: PLR2004 - the seven datasets D-114 renamed
    assert _job('job_project "no-dataset-file_here"') == "no-dataset-file_here"
