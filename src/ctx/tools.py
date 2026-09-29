"""The read-level tools ctx runs by subprocess: where each one is, and that it is the pinned one.

They live in the `reads` pixi environment, outside the analysis environment's solve (D-027),
so a run started from `default` or `laptop` finds them in the sibling environment. A tool is
looked up in ``$CTX_TOOLS`` (a bin directory) first, then in that sibling, then on PATH, and
refused unless it reports the version `pixi.toml` pins; a figure made with another FastQC is
not the figure this course describes.
"""

import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from functools import cache
from pathlib import Path

TOOLS_ENV = "CTX_TOOLS"
# pixi keeps each environment in a sibling folder: .pixi/envs/<name>.
SIBLING = Path(sys.prefix).parent / "reads" / "bin"


class ToolError(RuntimeError):
    """A tool is missing, reports another version, or failed; the message says which."""


@dataclass(frozen=True)
class Tool:
    """One pinned tool: its pixi package, how to ask its version, and what the answer holds."""

    package: str
    version: str
    probe: tuple[str, ...]
    pattern: str


# The versions are pixi.toml's [feature.reads.dependencies]; tests/test_repo.py holds them equal.
TOOLS = {
    "fastqc": Tool("fastqc", "0.12.1", ("--version",), r"FastQC v(\S+)"),
    "fastp": Tool("fastp", "1.1.0", ("--version",), r"fastp (\S+)"),
    "seqkit": Tool("seqkit", "2.13.0", ("version",), r"seqkit v(\S+)"),
    "multiqc": Tool("multiqc", "1.35", ("--version",), r"version (\S+)"),
    "STAR": Tool("star", "2.7.11b", ("--version",), r"^(\S+)$"),
    "samtools": Tool("samtools", "1.22.1", ("--version",), r"samtools (\S+)"),
    "featureCounts": Tool("subread", "2.1.1", ("-v",), r"featureCounts v(\S+)"),
    "salmon": Tool("salmon", "2.7.0", ("--version",), r"salmon (\S+)"),
    "trim_galore": Tool("trim-galore", "2.3.0", ("--version",), r"trim_galore (\S+)"),
    "hisat2": Tool("hisat2", "2.2.3", ("--version",), r"version (\S+)"),
}


def locate(name: str) -> Path:
    """The tool's executable: ``$CTX_TOOLS``, else the sibling `reads` environment, else PATH.

    Raises:
        ToolError: If no executable of that name is found anywhere.
    """
    folders = [Path(os.environ[TOOLS_ENV])] if os.environ.get(TOOLS_ENV) else []
    for folder in [*folders, SIBLING]:
        if (folder / name).is_file():
            return folder / name
    found = shutil.which(name)
    if found is None:
        raise ToolError(f"{name} not found; install the tools with `pixi install -e reads`")
    return Path(found)


@cache
def require(name: str) -> Path:
    """The tool's executable, after checking it reports the pinned version.

    Raises:
        ToolError: If it is missing or reports another version.
    """
    tool, path = TOOLS[name], locate(name)
    done = subprocess.run(
        [str(path), *tool.probe], capture_output=True, text=True, check=False,
        env=_env(path), timeout=120,
    )  # fmt: skip
    found = re.search(tool.pattern, done.stdout + done.stderr, re.MULTILINE)
    version = found.group(1) if found else "no version"
    if version != tool.version:
        raise ToolError(f"{path} reports {version}; this course pins {name} {tool.version}")
    return path


def _env(path: Path) -> dict[str, str]:
    """The tool's environment: its own bin folder first on PATH.

    Its helpers, the Java FastQC starts and the Python MultiQC runs on, then come from the
    tool's environment, not from the one ctx runs in.
    """
    return {**os.environ, "PATH": f"{path.parent}{os.pathsep}{os.environ.get('PATH', '')}"}


def run(name: str, args: list[str], *, cwd: Path, log: Path) -> None:
    """Run a pinned tool; its output goes to ``log``.

    Args:
        name: A key of ``TOOLS``.
        args: The arguments after the executable.
        cwd: The working directory.
        log: Where stdout and stderr go, kept beside the tool's own outputs.

    Raises:
        ToolError: If the tool is not the pinned one, or exits non-zero; the message ends
            with the last lines it wrote.
    """
    path = require(name)
    with log.open("w") as handle:
        done = subprocess.run([str(path), *args], cwd=cwd, stdout=handle,
                              stderr=subprocess.STDOUT, env=_env(path), check=False)  # fmt: skip
    if done.returncode != 0:
        tail = "\n".join(log.read_text().splitlines()[-15:])
        raise ToolError(f"{name} exited {done.returncode}; the end of {log}:\n{tail}")
