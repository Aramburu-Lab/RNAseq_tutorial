"""What this run is: its profile and root, the resources it used, the software it imported.

Profiles (D2, D-004): ``laptop`` is the default; ``hpc`` takes threads and memory from the
Slurm allocation and never guesses. The command syntax is the same under both.
"""

import os
import platform
import resource
import subprocess
import sys
from collections.abc import Mapping
from dataclasses import dataclass
from importlib import metadata
from pathlib import Path
from time import perf_counter

import ctx

PROFILES = ("laptop", "hpc")
LAPTOP_CORES = 4
THREAD_VARIABLES = ("POLARS_MAX_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS")
# Libraries that size their thread pools on import: the profile must be applied before them.
HEAVY = ("polars", "numpy", "scipy", "sklearn", "matplotlib", "seaborn", "plotly", "pydeseq2")
TRACK_A_MODULES = frozenset({"qc", "process", "quant"})


@dataclass(frozen=True)
class Budget:
    """A track's limits (CONVENTIONS.md, Stack): wall minutes and peak RAM on LAPTOP_CORES."""

    minutes: int
    ram_gb: int


BUDGETS = {"A": Budget(minutes=15, ram_gb=6), "B": Budget(minutes=10, ram_gb=4)}


@dataclass(frozen=True)
class Settings:
    """The resolved runtime of one command."""

    profile: str
    root: Path
    threads: int
    memory_bytes: int | None
    started: float


@dataclass(frozen=True)
class Usage:
    """Measured wall time and peak resident set sizes, in seconds and bytes."""

    wall_seconds: float
    peak_rss_self: int
    peak_rss_children: int


def budget(module: str) -> Budget:
    """The budget of the track a module belongs to (A: read-level modules, B: the rest)."""
    return BUDGETS["A" if module in TRACK_A_MODULES else "B"]


def resolve_profile(flag: str | None, env: Mapping[str, str]) -> str:
    """Pick the profile: the flag, then ``CTX_PROFILE``, then ``laptop``.

    Args:
        flag: Value of ``--profile``, or None when not given.
        env: The environment to read ``CTX_PROFILE`` from.

    Returns:
        One of ``PROFILES``.

    Raises:
        ValueError: For an unknown profile name.
    """
    name = flag or env.get("CTX_PROFILE") or "laptop"
    if name not in PROFILES:
        raise ValueError(f"unknown profile {name!r}; choose one of {PROFILES}")
    return name


def resolve_root(flag: Path | None, profile: str, env: Mapping[str, str], cwd: Path) -> Path:
    """Pick the output root: the flag, then ``CTX_ROOT``, then the working directory on laptop.

    Args:
        flag: Value of ``--root``, or None when not given.
        profile: The resolved profile.
        env: The environment to read ``CTX_ROOT`` from.
        cwd: The working directory, the laptop default.

    Returns:
        The absolute root.

    Raises:
        ValueError: On ``hpc`` with no explicit root, because a guessed root on a cluster writes
            where nobody looks (settings/paths.toml).
    """
    chosen = flag or (Path(env["CTX_ROOT"]) if env.get("CTX_ROOT") else None)
    if chosen is None and profile == "hpc":
        raise ValueError("the hpc profile needs an explicit root: pass --root or set CTX_ROOT")
    return (chosen or cwd).resolve()


def resources_for(profile: str, env: Mapping[str, str]) -> tuple[int, int | None]:
    """Threads and memory for a profile.

    Args:
        profile: The resolved profile.
        env: The environment; on ``hpc`` the Slurm allocation variables (sbatch(1), OUTPUT
            ENVIRONMENT VARIABLES).

    Returns:
        Threads, and the memory allocation in bytes (None on laptop, which has a budget instead).

    Raises:
        ValueError: On ``hpc`` outside a Slurm job, naming the missing variable.
    """
    if profile == "laptop":
        return min(LAPTOP_CORES, os.cpu_count() or 1), None
    if not env.get("SLURM_CPUS_PER_TASK"):
        raise ValueError("hpc profile outside a Slurm job: SLURM_CPUS_PER_TASK is not set")
    threads = int(env["SLURM_CPUS_PER_TASK"])
    if env.get("SLURM_MEM_PER_CPU"):
        return threads, threads * int(env["SLURM_MEM_PER_CPU"]) * 2**20
    if env.get("SLURM_MEM_PER_NODE"):
        return threads, int(env["SLURM_MEM_PER_NODE"]) * 2**20
    raise ValueError("hpc profile: neither SLURM_MEM_PER_CPU nor SLURM_MEM_PER_NODE is set")


def start(profile_flag: str | None, root_flag: Path | None) -> Settings:
    """Resolve the profile and root and set thread limits; call before any heavy import.

    Args:
        profile_flag: Value of ``--profile``.
        root_flag: Value of ``--root``.

    Returns:
        The settings, with the start time for ``measure()``.

    Raises:
        RuntimeError: If a library in ``HEAVY`` was imported first; its pool is already sized.
    """
    started = perf_counter()
    if early := [name for name in HEAVY if name in sys.modules]:
        raise RuntimeError(f"profile applied after importing {early}; thread limits would not hold")
    profile = resolve_profile(profile_flag, os.environ)
    root = resolve_root(root_flag, profile, os.environ, Path.cwd())
    threads, memory = resources_for(profile, os.environ)
    for variable in THREAD_VARIABLES:
        os.environ[variable] = str(threads)
    return Settings(profile, root, threads, memory, started)


def measure(started: float) -> Usage:
    """Wall time since ``started`` and peak RSS of this process and of its largest reaped child.

    Valid once per fresh process: both peaks are lifetime maxima, so a second command in the
    same process would inherit the first one's.

    Args:
        started: A ``time.perf_counter()`` value.

    Returns:
        The measured usage.
    """
    # ru_maxrss is bytes on macOS and kibibytes on Linux (getrusage(2) on each).
    scale = 1 if sys.platform == "darwin" else 1024
    return Usage(
        wall_seconds=perf_counter() - started,
        peak_rss_self=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale,
        peak_rss_children=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * scale,
    )


def versions() -> dict[str, str]:
    """Versions of Python, ctx and every distribution behind a module imported so far.

    Called at write time, so a README lists what the run actually used, not what is installed.

    Returns:
        Distribution name to version, sorted by name.
    """
    found = {"python": platform.python_version(), "ctx": ctx.__version__}
    owners = metadata.packages_distributions()
    for top in {name.partition(".")[0] for name in list(sys.modules)}:
        for dist in owners.get(top, []):
            found.setdefault(dist, metadata.version(dist))
    return dict(sorted(found.items()))


def source() -> str:
    """Where the running ctx came from: a git checkout (with its commit) or an installed wheel.

    The commit is read only when ctx is imported from ``<checkout>/src/ctx``. Walking up from
    the working directory would report whatever repository the user happens to stand in.

    Returns:
        ``"git <sha>"``, ``"git <sha> (uncommitted changes)"``, ``"installed wheel"`` or
        ``"git unavailable"``.
    """
    init = Path(ctx.__file__).resolve()
    checkout = init.parents[2]
    if not (checkout / ".git").exists() or init != checkout / "src" / "ctx" / "__init__.py":
        return "installed wheel"
    try:
        sha = _git(checkout, "rev-parse", "HEAD")
        # Untracked files count only under src/: a new module that ran is a change to the code.
        dirty = _git(checkout, "status", "--porcelain", "--", "src")
    except (OSError, subprocess.CalledProcessError):
        return "git unavailable"
    return f"git {sha}" + (" (uncommitted changes)" if dirty else "")


def _git(checkout: Path, *args: str) -> str:
    """Run one git command in the checkout and return its stripped output."""
    return subprocess.run(
        ["git", "-C", str(checkout), *args], capture_output=True, text=True, check=True
    ).stdout.strip()
