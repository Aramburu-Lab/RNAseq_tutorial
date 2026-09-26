"""The built wheel works on its own, with the editable checkout out of reach."""

import subprocess
import sys
import textwrap
import tomllib
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def build_wheel(tmp_path: Path) -> Path:
    """Build the wheel and unzip it; return the directory that acts as site-packages."""
    dist = tmp_path / "dist"
    subprocess.run(
        [sys.executable, "-m", "hatchling", "build", "-t", "wheel", "-d", str(dist)],
        cwd=ROOT,
        check=True,
        capture_output=True,
    )
    (wheel,) = dist.glob("*.whl")
    site = tmp_path / "site"
    with zipfile.ZipFile(wheel) as zf:
        zf.extractall(site)
    return site


def run_from_wheel(site: Path, code: str, cwd: Path) -> str:
    """Run code in a child that first removes the editable src/ from sys.path."""
    prelude = textwrap.dedent(f"""
        import sys
        src = {str(ROOT / "src")!r}
        assert src in sys.path, "editable src/ not on sys.path: removing it would prove nothing"
        sys.path[:] = [{str(site)!r}] + [p for p in sys.path if p != src]
    """)
    return subprocess.run(
        [sys.executable, "-c", prelude + textwrap.dedent(code)],
        cwd=cwd,
        check=True,
        capture_output=True,
        text=True,
    ).stdout


def test_wheel_runs_outside_repo(tmp_path):
    site = build_wheel(tmp_path)
    (entry_points,) = site.glob("*.dist-info/entry_points.txt")
    assert "ctx = ctx.cli:app" in entry_points.read_text()

    out = run_from_wheel(
        site,
        """
        import ctx
        from typer.testing import CliRunner
        from ctx.cli import app
        print(ctx.__file__)
        print(CliRunner().invoke(app, ["--version"]).output.strip())
        """,
        cwd=tmp_path,
    ).split()
    assert Path(out[0]).is_relative_to(site)
    declared = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["version"]
    assert out[1] == declared
