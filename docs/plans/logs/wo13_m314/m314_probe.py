"""M3.14's timing probe, a pytest plugin: load with ``-p m314_probe``, ``M314_OUT`` set.

Each process (the controller and every xdist worker) appends JSON lines to
``$M314_OUT/<worker>.jsonl``: each phase's duration, each test's wall and CPU (self and reaped
children), and each ``conftest.once`` call, the build timed apart from the wait. Collection
writes each test's fixtures and the scope of every one it reaches.
"""

import json
import os
import resource
import sys
import time
from pathlib import Path

import pytest

OUT = Path(os.environ["M314_OUT"])
WORKER = os.environ.get("PYTEST_XDIST_WORKER", "main")


def _write(record: dict) -> None:
    with (OUT / f"{WORKER}.jsonl").open("a") as handle:
        handle.write(json.dumps(record) + "\n")


def _cpu() -> float:
    kinds = (resource.RUSAGE_SELF, resource.RUSAGE_CHILDREN)
    own, children = (resource.getrusage(kind) for kind in kinds)
    return own.ru_utime + own.ru_stime + children.ru_utime + children.ru_stime


def pytest_sessionstart() -> None:
    """Wrap ``conftest.once`` so each build is timed apart from the wait for it."""
    # Fixtures look `once` up in conftest's namespace when they run, so the wrapper is seen.
    conftest = sys.modules.get("conftest")
    if conftest is None or getattr(conftest, "_m314", False):
        return
    original = conftest.once

    def timed(factory, name, build):
        def built(folder):
            started = time.perf_counter()
            try:
                return build(folder)
            finally:
                wall = time.perf_counter() - started
                _write({"kind": "build", "name": name, "wall": wall, "worker": WORKER})

        started = time.perf_counter()
        value = original(factory, name, built)
        wall = time.perf_counter() - started
        _write({"kind": "once", "name": name, "wall": wall, "worker": WORKER})
        return value

    conftest.once = timed
    conftest._m314 = True


def pytest_runtest_logreport(report: pytest.TestReport) -> None:
    """One line per setup, call and teardown."""
    _write({"kind": "phase", "nodeid": report.nodeid, "when": report.when,
            "dur": report.duration, "outcome": report.outcome, "worker": WORKER})  # fmt: skip


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_protocol(item: pytest.Item):  # noqa: ANN201 - a hook wrapper's generator
    """One line per test: its wall time and the CPU of this process and its reaped children."""
    started, cpu = time.perf_counter(), _cpu()
    yield
    _write({"kind": "test", "nodeid": item.nodeid, "wall": time.perf_counter() - started,
            "cpu": _cpu() - cpu, "worker": WORKER})  # fmt: skip


def pytest_collection_modifyitems(items: list[pytest.Item]) -> None:
    """Each test's fixture names, and the scope of each fixture it reaches."""
    names = {i.nodeid: list(i.fixturenames) for i in items}
    (OUT / f"fixtures_{WORKER}.json").write_text(json.dumps(names))
    scopes = {}
    for item in items:
        defs = item._fixtureinfo.name2fixturedefs  # noqa: SLF001 - no public view
        # A module fixture imported into another file is another FixtureDef: keyed by function.
        scopes[item.nodeid] = {
            (f"{d[-1].func.__module__}.{n}" if d[-1].scope == "module" else n): d[-1].scope
            for n, d in defs.items()
            if d and d[-1].scope != "function"
        }
    (OUT / f"scopes_{WORKER}.json").write_text(json.dumps(scopes))
