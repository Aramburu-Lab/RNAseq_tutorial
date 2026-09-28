"""Runtime profiles, roots, thread limits and resource measurement."""

import json
import os
import re
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from ctx.runtime import (
    BUDGETS,
    HEAVY,
    LAPTOP_CORES,
    resolve_profile,
    resolve_root,
    resources_for,
)

ROOT = Path(__file__).resolve().parents[1]
MIB = 2**20


def _python(code: str, env: dict[str, str] | None = None) -> str:
    return subprocess.run(
        [sys.executable, "-c", textwrap.dedent(code)],
        capture_output=True,
        text=True,
        check=True,
        env={**os.environ, **(env or {})},
    ).stdout


@pytest.mark.parametrize(
    ("flag", "env", "expected"),
    [
        ("hpc", {"CTX_PROFILE": "laptop"}, "hpc"),
        (None, {"CTX_PROFILE": "hpc"}, "hpc"),
        (None, {}, "laptop"),
    ],
)
def test_profile_precedence(flag, env, expected):
    assert resolve_profile(flag, env) == expected


def test_unknown_profile_raises():
    with pytest.raises(ValueError, match="unknown profile"):
        resolve_profile("cloud", {})


def test_root_precedence(tmp_path):
    flag, env, cwd = tmp_path / "flag", tmp_path / "env", tmp_path / "cwd"
    ctx_root = {"CTX_ROOT": str(env)}
    assert resolve_root(flag, "laptop", ctx_root, cwd) == flag.resolve()
    assert resolve_root(None, "laptop", ctx_root, cwd) == env.resolve()
    assert resolve_root(None, "laptop", {}, cwd) == cwd.resolve()
    assert resolve_root(None, "hpc", ctx_root, cwd) == env.resolve()
    with pytest.raises(ValueError, match="--root or set CTX_ROOT"):
        resolve_root(None, "hpc", {}, cwd)


def test_hpc_requires_slurm():
    with pytest.raises(ValueError, match="SLURM_CPUS_PER_TASK"):
        resources_for("hpc", {})
    with pytest.raises(ValueError, match="SLURM_MEM_PER_CPU nor SLURM_MEM_PER_NODE"):
        resources_for("hpc", {"SLURM_CPUS_PER_TASK": "7"})
    per_cpu = {"SLURM_CPUS_PER_TASK": "7", "SLURM_MEM_PER_CPU": "1000"}
    assert resources_for("hpc", per_cpu) == (7, 7 * 1000 * MIB)
    per_node = {"SLURM_CPUS_PER_TASK": "7", "SLURM_MEM_PER_NODE": "5000"}
    assert resources_for("hpc", per_node) == (7, 5000 * MIB)


@pytest.mark.parametrize(("cores", "threads"), [(2, 2), (16, LAPTOP_CORES)])
def test_laptop_threads(monkeypatch, cores, threads):
    monkeypatch.setattr(os, "cpu_count", lambda: cores)
    assert resources_for("laptop", {}) == (threads, None)


def test_budgets_match_conventions():
    line = next(
        line for line in (ROOT / "CONVENTIONS.md").read_text().splitlines() if "Budgets:" in line
    )
    found = re.search(
        r"Track A < (\d+) min and < (\d+) GB peak RAM on (\d+) cores; "
        r"Track B < (\d+) min and < (\d+) GB",
        line,
    )
    assert found, line
    a_min, a_gb, cores, b_min, b_gb = map(int, found.groups())
    assert (BUDGETS["A"].minutes, BUDGETS["A"].ram_gb, LAPTOP_CORES) == (a_min, a_gb, cores)
    assert (BUDGETS["B"].minutes, BUDGETS["B"].ram_gb) == (b_min, b_gb)


def test_import_is_light():
    loaded = _python(
        """
        import sys
        import ctx.cli
        from ctx.runtime import HEAVY
        print(",".join(m for m in HEAVY if m in sys.modules))
        """
    )
    assert loaded.strip() == ""
    assert len(HEAVY) > 0


def test_start_refuses_after_heavy_import():
    with pytest.raises(subprocess.CalledProcessError) as failed:
        _python("import numpy\nfrom ctx.runtime import start\nstart(None, None)")
    assert "profile applied after importing ['numpy']" in failed.value.stderr


def test_threads_applied_before_import(tmp_path):
    threads = 2 if os.cpu_count() != 2 else 3
    slurm = {"CTX_PROFILE": "hpc", "SLURM_CPUS_PER_TASK": str(threads), "SLURM_MEM_PER_CPU": "100"}
    out = _python(
        f"""
        from pathlib import Path
        from ctx.runtime import start
        start(None, Path({str(tmp_path)!r}))
        # numpy loads OpenBLAS, and threadpoolctl reports only libraries already loaded.
        import numpy, polars, threadpoolctl
        pools = threadpoolctl.threadpool_info()
        blas = [p["num_threads"] for p in pools if p["user_api"] == "blas"]
        print(polars.thread_pool_size(), blas)
        """,
        env=slurm,
    )
    assert out.split(maxsplit=1) == [str(threads), f"[{threads}]\n"]


PARENT = """
import json, subprocess, sys
from time import perf_counter
from ctx.runtime import measure
started = perf_counter()
subprocess.run([sys.executable, "-c", {child!r}], check=True)
usage = measure(started)
print(json.dumps([usage.peak_rss_self, usage.peak_rss_children]))
"""
# One write per 4 KiB: a zero-filled bytearray does not raise RSS until its pages are touched.
TOUCH_256_MIB = "b = bytearray(256 * 2**20)\nb[::4096] = b'\\x01' * (len(b) // 4096)"


def test_child_peak_attributed():
    own, child = json.loads(_python(PARENT.format(child=TOUCH_256_MIB)))
    assert 256 * MIB <= child < 1024 * MIB
    assert own < 256 * MIB
    _, idle_child = json.loads(_python(PARENT.format(child="pass")))
    assert idle_child < 256 * MIB


LOKY_CHILD = textwrap.dedent(
    """
    import json, sys
    from time import perf_counter
    import numpy as np
    from joblib.externals.loky import get_reusable_executor
    from ctx.de.deseq import _shutdown_workers
    from ctx.runtime import measure

    def touch(mib):
        block = np.ones(mib * 2**20, dtype=np.uint8)
        return int(block.sum())

    started = perf_counter()
    get_reusable_executor(max_workers=1).submit(touch, 256).result()
    _shutdown_workers()
    usage = measure(started)
    print(json.dumps({"children": usage.peak_rss_children, "self": usage.peak_rss_self}))
    """
)


def test_loky_child_peak_attributed():
    # A fresh process, so no other child of the test run can set the children's peak.
    done = subprocess.run(
        [sys.executable, "-c", LOKY_CHILD], capture_output=True, text=True, check=True
    )
    usage = json.loads(done.stdout.strip().splitlines()[-1])
    assert usage["children"] >= 256 * 2**20, usage
    assert usage["children"] < 2**30, usage


def test_home_paths_shown_as_tilde():
    from ctx.runtime import retypable, shown_path  # noqa: PLC0415

    home = Path.home()
    assert shown_path(home / "ctx_tier1") == "~/ctx_tier1"
    assert shown_path(home) == "~"
    # A sibling whose name starts with the home path is not under it.
    assert shown_path(f"{home}2/x") == f"{home}2/x"
    command = retypable(["--root", str(home / "a b"), "de", "--x", "/opt/c d"])
    # The tilde stays outside the quotes, so a shell still expands it.
    assert command == "ctx --root ~/'a b' de --x '/opt/c d'"
    assert str(home) not in command
