"""M3.14: one probed suite run as tables. Arguments: run folder, label, command, output .md.

The run folder holds `time.txt` (`/usr/bin/time -l`), `out.txt`, the load files, `bt/`
(the basetemp, with `session_runs.log`) and `probe/` (m314_probe.py's files).
"""

import json
import re
import sys
from collections import defaultdict
from pathlib import Path

run, label, command, out = Path(sys.argv[1]), sys.argv[2], sys.argv[3], Path(sys.argv[4])
probe = run / "probe"
recs = [json.loads(x) for f in probe.glob("*.jsonl") for x in f.read_text().splitlines()]
workers = any(r["worker"] != "main" for r in recs)
recs = [r for r in recs if not (workers and r["worker"] == "main")]
fix, scopes = {}, {}
for f in probe.glob("fixtures_*.json"):
    fix.update(json.loads(f.read_text()))
for f in probe.glob("scopes_*.json"):
    scopes.update(json.loads(f.read_text()))
# Before the scope probe existed: conftest's session fixtures, by name.
CONFTEST = {
    "process_root",
    "process_rerun",
    "process_run",
    "quant_built",
    "quant_run",
    "wheel_site",
    "setup_run",
    "gather_run",
    "pipeline",
    "de_run",
    "function_run",
    "report_run",
    "qc_run",
    "clinical_run",
    "proteo_run",
    "genecard_run",
    "ncrna_run",
    "guide_run",
    "run_run",
}


BUILTIN = {"tmp_path_factory", "worker_id", "testrun_uid", "pytestconfig"}


def bare(nodeid: str) -> str:
    """A node ID without the ``@group`` suffix xdist's loadgroup adds."""
    return nodeid.split("@", 1)[0]


def scope_of(nodeid: str) -> dict:
    """The scopes of the fixtures a test reaches, beyond function scope."""
    return scopes.get(nodeid) or scopes.get(bare(nodeid)) or {}


def shared(nodeid: str) -> str:
    """The session fixtures a test reaches, then its module ones."""
    if scope_of(nodeid):
        got = scope_of(nodeid)
        names = sorted(n for n, s in got.items() if s == "session" and n not in BUILTIN)
        mods = sorted(n for n, s in got.items() if s == "module")
    else:
        names = sorted(CONFTEST & set(fix.get(nodeid, [])))
        mods = []
    param = re.search(r"\[(\w+)\]$", bare(nodeid))
    if param and param[1] in CONFTEST | {"across_root", "cross_root", "study_root", "made"}:
        names.append(f"{param[1]} (by request)")
    elif param and f"{param[1]}_run" in CONFTEST:
        names.append(f"{param[1]}_run (by request)")
    text = ", ".join(names) or "none"
    return text + (f"; module: {', '.join(mods)}" if mods else "")


time_txt = (run / "time.txt").read_text()
real, user, system = (
    float(x) for x in re.search(r"([\d.]+) real\s+([\d.]+) user\s+([\d.]+) sys", time_txt).groups()
)
rss = int(re.search(r"(\d+)\s+maximum resident set size", time_txt)[1])
summary = [
    x for x in (run / "out.txt").read_text().splitlines() if re.search(r"\d+ (passed|failed)", x)
][-1].strip("= ")
loads = [
    float(re.search(r"load averages?: ([\d.]+)", x)[1])
    for x in (run / "load_during.txt").read_text().splitlines()
    if "load" in x
]
before = (run / "load_before.txt").read_text().strip().split("load averages: ")[1]
after = (run / "load_after.txt").read_text().strip().split("load averages: ")[1]
log = run / "bt" / "session_runs.log"
built = log.read_text().split() if log.is_file() else []
builds = defaultdict(list)
waits = defaultdict(list)
for r in recs:
    if r["kind"] == "build":
        builds[r["name"]].append((r["worker"], r["wall"]))
    if r["kind"] == "once":
        waits[r["name"]].append((r["worker"], r["wall"]))
phases = sorted(
    (r for r in recs if r["kind"] == "phase" and r["when"] != "teardown"), key=lambda r: -r["dur"]
)
tests = [r for r in recs if r["kind"] == "test"]
cpu = {bare(r["nodeid"]): r["cpu"] for r in tests}
per_file, cpu_file = defaultdict(float), defaultdict(float)
for r in recs:
    if r["kind"] == "phase":
        per_file[bare(r["nodeid"]).split("::")[0]] += r["dur"]
for r in tests:
    cpu_file[bare(r["nodeid"]).split("::")[0]] += r["cpu"]
busy = defaultdict(float)
for r in recs:
    if r["kind"] == "phase":
        busy[r["worker"]] += r["dur"]

lines = [
    f"# {label}",
    "",
    f"Command: `{command}`",
    "",
    "| Measure | Value |",
    "| --- | --- |",
    f"| Result | {summary} |",
    f"| Wall time (`/usr/bin/time -l`) | {real:,.0f} s |",
    f"| CPU time, user + sys | {user:,.0f} + {system:,.0f} = {user + system:,.0f} s |",
    f"| Largest resident set | {rss:,} B |",
    f"| Load average, 1/5/15 min, at start | {before} |",
    f"| Load average, 1 min, during (every 30 s) | {min(loads):.1f} to {max(loads):.1f} |",
    f"| Load average, 1/5/15 min, at end | {after} |",
    f"| Worker busy time, least to most | {min(busy.values()):,.0f} to "
    f"{max(busy.values()):,.0f} s over {len(busy)} workers |",
    "",
    "## Session runs built through `conftest.once`",
    "",
    "Builds counts the run's `session_runs.log`, which `built_twice` reads; wall is the "
    "build on the worker that took the lock; the others waited on it.",
    "",
    "| Name | Builds | Built on | Build wall (s) | Workers that asked | Longest wait (s) |",
    "| --- | --- | --- | --- | --- | --- |",
]
for name in sorted(builds, key=lambda n: -max(w for _, w in builds[n])):
    (worker, wall), *_ = builds[name]
    lines.append(
        f"| {name} | {built.count(name)} | {worker} | {wall:.1f} | "
        f"{len(waits[name])} | {max(w for _, w in waits[name]):.1f} |"
    )
lines += [
    "",
    "## The 40 slowest setups and calls",
    "",
    "A setup's time includes building or waiting for the session runs it is the first on "
    "its worker to ask for. CPU is the whole test's, self and reaped children, where "
    "measured.",
    "",
    "| Wall (s) | Phase | Worker | Test | CPU of the test (s) | Shared fixtures it reaches |",
    "| --- | --- | --- | --- | --- | --- |",
]
for r in phases[:40]:
    c = cpu.get(bare(r["nodeid"]))
    lines.append(
        f"| {r['dur']:.1f} | {r['when']} | {r['worker']} | `{bare(r['nodeid'])}` | "
        f"{'' if c is None else f'{c:.1f}'} | {shared(r['nodeid'])} |"
    )
# Module fixtures: built once per (worker, requesting file), read from each test's record.
module_builds = defaultdict(set)
for r in recs:
    if r["kind"] == "phase" and r["when"] == "setup":
        for name, scope in scope_of(r["nodeid"]).items():
            if scope == "module":
                module_builds[name].add((r["worker"], bare(r["nodeid"]).split("::")[0]))
lines += [
    "",
    "## Module fixtures built more than once",
    "",
    "A module fixture is built once per worker and per file that requests it.",
    "",
    "| Fixture | Builds | Files that request it, with the workers |",
    "| --- | --- | --- |",
]
for name, where in sorted(module_builds.items(), key=lambda kv: (-len(kv[1]), kv[0])):
    if len(where) > 1:
        files = sorted(f"{f.removeprefix('tests/')} ({w})" for w, f in where)
        lines.append(f"| {name} | {len(where)} | {'; '.join(files)} |")
once_only = sum(len(w) == 1 for w in module_builds.values())
lines += ["", f"Module fixtures built once: {once_only}."]
lines += [
    "",
    "## Files by summed setup and call time",
    "",
    "| File | Summed wall (s) | CPU (s) |",
    "| --- | --- | --- |",
]
for f, s in sorted(per_file.items(), key=lambda kv: -kv[1])[:20]:
    lines.append(f"| `{f}` | {s:.1f} | {cpu_file[f]:.1f} |" if tests else f"| `{f}` | {s:.1f} | |")
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("\n".join(lines) + "\n")
print(out, len(lines))
