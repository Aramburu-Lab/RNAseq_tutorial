"""ctx.tools finds the pinned read-level tools, and refuses any other version (D-027)."""

import os
import stat
import tomllib
from pathlib import Path

import pytest

from ctx import tools

ROOT = Path(__file__).resolve().parents[1]


def test_pins_match_pixi_toml():
    """One pin per tool: the table ctx checks against equals what pixi installs."""
    pinned = tomllib.loads((ROOT / "pixi.toml").read_text())["feature"]["reads"]["dependencies"]
    assert {t.package: f"=={t.version}" for t in tools.TOOLS.values()} == {
        package: pin for package, pin in pinned.items() if package in
        {t.package for t in tools.TOOLS.values()}
    }  # fmt: skip
    assert {t.package for t in tools.TOOLS.values()} <= set(pinned)


def test_every_tool_is_the_pinned_one():
    if not tools.SIBLING.is_dir():
        pytest.fail(f"{tools.SIBLING} missing: run `pixi install -e reads` (D-027)")
    tools.require.cache_clear()
    for name in tools.TOOLS:
        assert tools.require(name).parent == tools.SIBLING, name


def _stub(folder: Path, name: str, answer: str) -> None:
    script = folder / name
    script.write_text(f"#!/bin/sh\necho '{answer}'\n")
    script.chmod(script.stat().st_mode | stat.S_IXUSR)


def test_another_version_is_refused(tmp_path, monkeypatch):
    _stub(tmp_path, "fastqc", "FastQC v0.11.9")
    monkeypatch.setenv(tools.TOOLS_ENV, str(tmp_path))
    tools.require.cache_clear()
    try:
        with pytest.raises(tools.ToolError, match="reports 0.11.9; this course pins fastqc 0.12.1"):
            tools.require("fastqc")
    finally:
        tools.require.cache_clear()


def test_ctx_tools_comes_first(tmp_path, monkeypatch):
    _stub(tmp_path, "fastqc", "FastQC v0.12.1")
    monkeypatch.setenv(tools.TOOLS_ENV, str(tmp_path))
    tools.require.cache_clear()
    try:
        assert tools.require("fastqc") == tmp_path / "fastqc"
    finally:
        tools.require.cache_clear()


def test_a_failed_run_says_why(tmp_path):
    if not tools.SIBLING.is_dir():
        pytest.fail(f"{tools.SIBLING} missing: run `pixi install -e reads` (D-027)")
    tools.require.cache_clear()
    with pytest.raises(tools.ToolError, match="seqkit exited"):
        tools.run("seqkit", ["stats", str(tmp_path / "absent.fq")], cwd=tmp_path,
                  log=tmp_path / "seqkit.log")  # fmt: skip
    assert os.path.getsize(tmp_path / "seqkit.log") > 0
