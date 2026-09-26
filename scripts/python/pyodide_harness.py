#!/usr/bin/env python3
# /// codebase-meta
# exempt: >-
#   Helper for spike_pyodide.py, not a pipeline step: it runs one Python function inside Pyodide
#   under Node, or in a fresh CPython process, and reports timings, downloads and errors.
# ///
"""Run one self-contained Python function in Pyodide (Node) or CPython, fresh each time.

The function travels to Pyodide as source text, so it must import what it needs inside its
body. Every run gets a fresh package cache, so the download log is the full install cost.
"""

import inspect
import json
import os
import platform
import signal
import subprocess
import sys
import time
from collections.abc import Callable, Mapping, Sequence
from importlib import metadata
from pathlib import Path
from typing import Any
from urllib.parse import unquote, urlparse

RUNNER_JS = r"""
import fs from "node:fs";
const cfg = JSON.parse(fs.readFileSync(process.argv[2], "utf8"));
const { loadPyodide } = await import(cfg.pyodide_mjs);
const out = { phases_s: {}, fetches: [] };
const save = () => fs.writeFileSync(cfg.out_json, JSON.stringify(out));
let phase = "load";
const realFetch = globalThis.fetch;
// Pyodide and micropip both download through fetch, so every wheel's bytes are counted here.
globalThis.fetch = async (input, init) => {
  const response = await realFetch(input, init);
  const body = await response.clone().arrayBuffer();
  const url = typeof input === "string" ? input : input.url;
  out.fetches.push({ phase, url, bytes: body.byteLength });
  return response;
};
async function timed(name, fn) {
  phase = name;
  const t0 = performance.now();
  try {
    return await fn();
  } finally {
    out.phases_s[name] = (performance.now() - t0) / 1000;
    save();
  }
}
const py = await timed("load", () => loadPyodide({ packageCacheDir: cfg.cache_dir }));
py.globals.set("cfg_json", JSON.stringify(cfg.python));
out.install = JSON.parse(await timed("install", async () => {
  await py.loadPackage(["micropip", ...cfg.python.loaded], { messageCallback: () => {} });
  return py.runPythonAsync(cfg.install_py);
}));
if (out.install.ready) {
  for (const [dest, src] of Object.entries(cfg.inputs)) py.FS.writeFile(dest, fs.readFileSync(src));
  out.run = JSON.parse(await timed("run", () => py.runPythonAsync(cfg.run_py)));
  if (out.run.status === "ok") fs.writeFileSync(cfg.results, py.FS.readFile(cfg.python.output));
}
save();
"""

INSTALL_PY = """
import json
from importlib import metadata
import micropip
cfg = json.loads(cfg_json)
steps = {}
def report(name, error=None):
    if error is None:
        return {"status": "installed", "version": metadata.version(name)}
    return {"status": "failed", "error": f"{type(error).__name__}: {error}"}
for requirement in cfg["requirements"]:
    name = requirement.split("==")[0]
    try:
        await micropip.install(requirement, keep_going=True)
        steps[name] = report(name)
    except Exception as error:
        steps[name] = report(name, error)
for name in cfg["loaded"]:
    try:
        __import__(name)
        steps[name] = report(name)
    except Exception as error:
        steps[name] = report(name, error)
if cfg["absent"]:
    try:
        __import__(cfg["absent"])
        steps["absent"] = {"status": "importable"}
    except ImportError as error:
        steps["absent"] = {"status": "not importable", "error": f"{type(error).__name__}: {error}"}
steps["distributions"] = {k: [v.version, v.source] for k, v in micropip.list().items()}
needed = [r.split("==")[0] for r in cfg["requirements"]] + cfg["loaded"]
steps["ready"] = all(steps[n]["status"] == "installed" for n in needed) and (
    steps.get("absent", {}).get("status") != "importable"
)
json.dumps(steps)
"""

RUN_PY = """
import json
import time
import traceback
cfg = json.loads(cfg_json)
namespace = {}
exec(cfg["source"], namespace)
start = time.perf_counter()
try:
    namespace[cfg["function"]](*cfg["args"])
    result = {"status": "ok"}
except Exception as error:
    result = {"status": "error", "error": f"{type(error).__name__}: {error}"}
    result["traceback"] = traceback.format_exc()
result["python_s"] = time.perf_counter() - start
json.dumps(result)
"""


def timed_process(cmd: Sequence[str], log: Path, timeout_s: float) -> dict[str, Any]:
    """Run a command under ``/usr/bin/time`` and kill its whole process group at the timeout.

    Args:
        cmd: The command.
        log: Where its stdout and stderr go, the time report included.
        timeout_s: Wall-clock limit.

    Returns:
        ``status`` (``exited`` or ``timeout``), ``returncode``, ``wall_s`` and
        ``peak_rss_bytes``: ``/usr/bin/time -l`` reports bytes on macOS, ``-v`` KiB on Linux.
    """
    flag = "-l" if sys.platform == "darwin" else "-v"
    start = time.perf_counter()
    with log.open("w") as sink:
        proc = subprocess.Popen(
            ["/usr/bin/time", flag, *cmd], stdout=sink, stderr=sink, start_new_session=True
        )
        try:
            code = proc.wait(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            # Killing /usr/bin/time alone would orphan node; the whole session goes instead.
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
            return {"status": "timeout", "wall_s": time.perf_counter() - start}
    result: dict[str, Any] = {"status": "exited", "returncode": code}
    result["wall_s"] = time.perf_counter() - start
    for line in log.read_text().splitlines():
        if line.strip().endswith("maximum resident set size"):
            result["peak_rss_bytes"] = int(line.split()[0])
        elif "Maximum resident set size (kbytes)" in line:
            result["peak_rss_bytes"] = 1024 * int(line.rsplit(":", 1)[1])
    return result


def in_pyodide(job: Mapping[str, Any], run_dir: Path, node_dir: Path, timeout_s: float) -> dict:
    """Run ``job["function"]`` in a fresh Pyodide under Node, with a fresh package cache.

    Args:
        job: ``function``; ``args`` (paths inside Pyodide); ``inputs`` ({Pyodide path: host
            path}); ``output`` (the Pyodide path the function writes); ``requirements`` for
            micropip; ``loaded``, Pyodide packages to load and import; ``absent``, a module
            that must not import, or "".
        run_dir: A new directory for this run.
        node_dir: The npm prefix holding ``node_modules/pyodide``.
        timeout_s: Wall-clock limit for the whole process.

    Returns:
        ``phases_s``, ``fetches``, ``install``, ``run`` as the runner saved them (partial after
        a timeout), plus ``process`` and ``results``, the host path of the function's output.
    """
    run_dir.mkdir()
    python = {k: job[k] for k in ("args", "requirements", "loaded", "absent", "output")}
    python |= {"source": inspect.getsource(job["function"]), "function": job["function"].__name__}
    config = {
        "pyodide_mjs": str(node_dir / "node_modules" / "pyodide" / "pyodide.mjs"),
        "cache_dir": str(run_dir / "cache"),
        "inputs": {k: str(v) for k, v in job["inputs"].items()},
        "install_py": INSTALL_PY,
        "run_py": RUN_PY,
        "python": python,
        "out_json": str(run_dir / "out.json"),
        "results": str(run_dir / "results.csv"),
    }
    (run_dir / "config.json").write_text(json.dumps(config))
    (run_dir / "runner.mjs").write_text(RUNNER_JS)
    cmd = ["node", str(run_dir / "runner.mjs"), str(run_dir / "config.json")]
    process = timed_process(cmd, run_dir / "process.log", timeout_s)
    out = run_dir / "out.json"
    reported = json.loads(out.read_text()) if out.exists() else {}
    return {**reported, "process": process, "results": run_dir / "results.csv"}


def in_cpython(
    script: Path,
    function: Callable[..., None],
    args: Sequence[Any],
    run_dir: Path,
    timeout_s: float,
) -> dict:
    """Run a function in a fresh CPython process: ``<script> --reference <config>``.

    Args:
        script: A script whose ``--reference`` flag calls ``reference()`` with its functions.
        function: The function, looked up by name in the child.
        args: Its arguments, host paths included.
        run_dir: A new directory for this run.
        timeout_s: Wall-clock limit.

    Returns:
        What ``reference()`` wrote (status, error, python_s, versions), plus ``process``.
    """
    run_dir.mkdir()
    config = run_dir / "config.json"
    config.write_text(json.dumps({"function": function.__name__, "args": list(args)}))
    cmd = [sys.executable, str(script), "--reference", str(config)]
    process = timed_process(cmd, run_dir / "process.log", timeout_s)
    out = run_dir / "out.json"
    return {**(json.loads(out.read_text()) if out.exists() else {}), "process": process}


def reference(
    config: Path, functions: Sequence[Callable[..., None]], report: Sequence[str]
) -> None:
    """The child side of ``in_cpython``: run one function and write ``out.json`` beside config.

    Args:
        config: JSON with ``function`` and ``args``.
        functions: The functions the config may name.
        report: Distributions whose installed versions go into ``out.json``.
    """
    spec = json.loads(config.read_text())
    function = {f.__name__: f for f in functions}[spec["function"]]
    start = time.perf_counter()
    try:
        function(*spec["args"])
        result: dict[str, Any] = {"status": "ok"}
    except Exception as error:  # the record carries the error, as the Pyodide arm's does
        result = {"status": "error", "error": f"{type(error).__name__}: {error}"}
    result["python_s"] = time.perf_counter() - start
    # conda installs typing_extensions under that name and PyPI as typing-extensions (PEP 503).
    installed = {d.metadata["Name"].lower().replace("_", "-") for d in metadata.distributions()}
    result["versions"] = {d: metadata.version(d) for d in report if d in installed}
    result["versions"]["python"] = platform.python_version()
    (config.parent / "out.json").write_text(json.dumps(result))


def wheels(fetches: Sequence[Mapping[str, Any]]) -> dict[str, int]:
    """Bytes per package file downloaded while installing; index pages and metadata excluded.

    Args:
        fetches: The runner's fetch log.

    Returns:
        File name to bytes.
    """
    files = {}
    for fetch in fetches:
        name = unquote(Path(urlparse(fetch["url"]).path).name)
        if fetch["phase"] == "install" and name.endswith((".whl", ".zip", ".tar")):
            files[name] = int(fetch["bytes"])
    return files
