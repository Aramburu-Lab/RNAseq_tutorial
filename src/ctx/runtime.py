"""What this run is: the software it imported and where ctx itself came from."""

import platform
import subprocess
import sys
from importlib import metadata
from pathlib import Path

import ctx


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
