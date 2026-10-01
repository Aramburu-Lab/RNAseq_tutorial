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


def test_stdout_is_kept_apart_from_messages(tmp_path):
    """A tool's result on stdout goes to its own file, its messages to the log, and a failed
    run's error names its messages.

    seqkit sample writes reads to stdout and "[INFO] sample by proportion" to stderr.

    Mutations caught: stdout ignored (the names land in the log, and the file is never made);
    messages sent into the stdout file (the log stays empty, and the error loses its tail)."""
    if not tools.SIBLING.is_dir():
        pytest.fail(f"{tools.SIBLING} missing: run `pixi install -e reads` (D-027)")
    tools.require.cache_clear()
    reads = ROOT / "tests" / "fixtures" / "chr22_se" / "se_a.fq"
    names, log = tmp_path / "names.txt", tmp_path / "seqkit.log"
    tools.run("seqkit", ["seq", "-n", str(reads)], cwd=tmp_path, log=log, stdout=names)
    lines = names.read_text().splitlines()
    assert len(lines) == 1000
    assert lines[0] not in log.read_text()
    sampled = tmp_path / "sampled.fq"
    tools.run("seqkit", ["sample", "-p", "0.5", "-s", "11", str(reads)], cwd=tmp_path, log=log,
              stdout=sampled)  # fmt: skip
    records = sampled.read_text().splitlines()
    assert records and len(records) % 4 == 0
    assert all(line.startswith("@") for line in records[::4])
    assert "sample by proportion" in log.read_text()
    with pytest.raises(tools.ToolError, match="no such file"):
        tools.run("seqkit", ["seq", "-n", str(reads), str(tmp_path / "absent.fq")],
                  cwd=tmp_path, log=log, stdout=names)  # fmt: skip
    assert "no such file" in log.read_text() and "no such file" not in names.read_text()
