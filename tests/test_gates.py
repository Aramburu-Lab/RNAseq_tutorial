"""The quality gates on src/ can fail, and ctx never imports labcore (D-003)."""

import ast
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src" / "ctx"
CONFIG = ROOT / "pyproject.toml"
GATES = {
    "mypy": ["mypy", "--config-file", str(CONFIG)],
    "interrogate": ["interrogate", "--config", str(CONFIG)],
}


def _gate(name: str, target: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", *GATES[name], str(target)],
        cwd=target.parent,
        capture_output=True,
        text=True,
        check=False,
    )


@pytest.mark.parametrize("name", GATES)
def test_gates_fail_on_violation(tmp_path, name):
    pkg = tmp_path / "ctx"
    shutil.copytree(SRC, pkg, ignore=shutil.ignore_patterns("__pycache__"))
    clean = _gate(name, pkg)
    assert clean.returncode == 0, clean.stdout + clean.stderr
    (pkg / "violation.py").write_text("def f(x):\n    return x\n")
    dirty = _gate(name, pkg)
    assert dirty.returncode != 0, f"{name} passed an untyped, undocumented function"


def test_hooks_run_the_gates_on_src():
    config = yaml.safe_load((ROOT / ".pre-commit-config.yaml").read_text())
    hooks = {h["id"]: h for repo in config["repos"] for h in repo["hooks"]}
    for hook_id, tool in (("mypy-src", "mypy"), ("interrogate-src", "interrogate")):
        hook = hooks[hook_id]
        assert hook["entry"].startswith(f"pixi run {tool}")
        # A files filter that matches nothing makes the hook pass on every commit.
        assert re.match(hook["files"], "src/ctx/cli.py")


def test_zero_collected_fails():
    inner = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", "-k", "zz_no_such_zz"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert inner.returncode == pytest.ExitCode.TESTS_FAILED, inner.stdout


def test_no_labcore_import():
    files = sorted(SRC.rglob("*.py"))
    assert files
    offenders = []
    for path in files:
        for node in ast.walk(ast.parse(path.read_text(), str(path))):
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [node.module or ""]
            else:
                continue
            offenders += [
                f"{path.name}:{node.lineno}" for n in names if n.split(".")[0] == "labcore"
            ]
    assert not offenders
