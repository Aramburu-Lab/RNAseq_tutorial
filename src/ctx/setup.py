"""M0, environment report: what this machine and environment are, and what a run costs.

Downloads nothing (Phase 1). Every fact on the page is measured in this run, not assumed.
"""

import hashlib
import json
import os
import platform
import urllib.request
from collections.abc import Callable, Sequence
from datetime import UTC, date, datetime
from pathlib import Path
from typing import Any

import numpy
import polars as pl
import threadpoolctl
from matplotlib import font_manager

from ctx import runtime
from ctx.dataset import Resource, load
from ctx.io.naming import MODULES, stem
from ctx.io.readme_companion import write_table
from ctx.report.builder import ModulePage, ParameterCard, Reproduction, SelfCheck, write
from ctx.viz.style import DATA_COLOURS, FONT_CHAIN, PROJECTOR_UNSAFE, THEMES, contrast

MODULE = "setup"
QUESTION = "Is this computer ready for the course, and what will a run cost it?"
CONCEPT = (
    "An analysis is only reproducible if you can say exactly which software ran it and on what. "
    "This page records the versions every later module relies on, how many threads and how "
    "much memory a run may use (the profile), and how long this run took and how much memory "
    "it peaked at. Every later report carries the same reproduction block."
)
PRACTICE = (
    "Keep this page with your results. When a number changes between two runs, compare the "
    "versions tables first: a changed library is the commonest reason."
)
THINK = (
    SelfCheck(
        "Why does the page list the peak memory of the largest child process separately?",
        "Tools such as aligners run as separate processes. Their memory is not part of this "
        "process's peak, and a laptop runs out of memory on the total, not on either half.",
    ),
    SelfCheck(
        "Why might a colleague get a different font on the same figure?",
        "The font chain takes the first installed family. A machine without Helvetica or Arial "
        "draws with the next one, which changes text width but not the data.",
    ),
)
MISTAKES = (
    "Comparing timings measured under different profiles: hpc and laptop use different "
    "thread counts.",
    "Reading a version from memory instead of from this table.",
)
VERSION_COLUMNS = {
    "package": ("none", "Python distribution imported by this run, or python/ctx."),
    "version": ("none", "Version string as installed."),
}


def _memory_bytes() -> int | None:
    """Physical memory, where the OS reports it."""
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (ValueError, OSError):
        return None


def _resolved_font() -> str:
    """The first family in the font chain that matplotlib can actually find here."""
    for family in FONT_CHAIN:
        try:
            font_manager.findfont(
                font_manager.FontProperties(family=family), fallback_to_default=False
            )
        except ValueError:
            continue
        return family
    return "none of the chain; matplotlib's default"


def _facts(settings: runtime.Settings) -> tuple[str, ...]:
    """The measured state of this machine and run, one sentence each."""
    numpy.zeros(1)  # loads BLAS, so its thread pool can be read
    blas = [p["num_threads"] for p in threadpoolctl.threadpool_info() if p["user_api"] == "blas"]
    memory = _memory_bytes()
    budget = runtime.budget(MODULE)
    contrasts = ", ".join(
        f"{c} {contrast(c, THEMES['white'].background):.2f}:1"
        + (" (not for projection)" if c in PROJECTOR_UNSAFE else "")
        for c in DATA_COLOURS
    )
    return (
        f"System: {platform.system()} {platform.release()}, {platform.machine()}; "
        f"Python {platform.python_version()}.",
        f"CPU cores reported: {os.cpu_count()}; physical memory: "
        + (f"{memory / 2**30:.1f} GiB." if memory else "not reported by this OS."),
        f"Profile: {settings.profile}; threads allowed: {settings.threads}; polars thread pool: "
        f"{pl.thread_pool_size()}; BLAS thread pools: {blas}.",
        f"Output root: {settings.root}.",
        f"Budget for this module: Track B, under {budget.minutes} min and {budget.ram_gb} GB.",
        f"Figure font in use: {_resolved_font()} (chain: {', '.join(FONT_CHAIN)}).",
        f"Contrast of each data colour against white: {contrasts}.",
    )


MANIFEST = Path("resources") / "manifest.json"
CHUNK = 2**20


class FetchError(RuntimeError):
    """A resource could not be fetched from any of its mirrors with the declared checksum."""


def _sha256(path: Path) -> str:
    """SHA-256 of a file, read in chunks."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(CHUNK):
            digest.update(block)
    return digest.hexdigest()


def _download(url: str, part: Path, opener: Callable[..., Any]) -> None:
    """Stream one URL to ``part``."""
    with opener(url, timeout=120) as response, part.open("wb") as out:
        while block := response.read(CHUNK):
            out.write(block)


def fetch(
    root: Path,
    resources: Sequence[Resource],
    opener: Callable[..., Any] = urllib.request.urlopen,
) -> list[dict[str, Any]]:
    """Fetch each resource under the root, checked by size and SHA-256 before it is kept.

    Each mirror is tried in order: download to ``<dest>.part``, check, then rename, so a
    file at ``dest`` is always a checked one. A file already there with the right hash is
    not fetched again. Every fetched or confirmed file is recorded in
    ``resources/manifest.json``.

    Args:
        root: The project root.
        resources: The dataset file's resources.
        opener: How to open a URL; tests pass their own.

    Returns:
        One manifest entry per resource.

    Raises:
        FetchError: If a destination leaves the root, or no mirror gives the declared bytes.
    """
    manifest_path = root / MANIFEST
    manifest = json.loads(manifest_path.read_text()) if manifest_path.is_file() else {}
    entries = []
    for resource in resources:
        dest = (root / resource.dest).resolve()
        if not dest.is_relative_to(root.resolve()):
            raise FetchError(f"{resource.name}: {resource.dest} is outside {root}")
        entry = manifest.get(resource.dest)
        if not (dest.is_file() and _sha256(dest) == resource.sha256 and entry):
            entry = _fetch_one(resource, dest, opener)
        manifest[resource.dest] = entry
        entries.append(entry)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return entries


def _fetch_one(resource: Resource, dest: Path, opener: Callable[..., Any]) -> dict[str, Any]:
    """Try each mirror until one gives the declared size and hash; leave no partial file."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    failures = []
    for url in resource.urls:
        try:
            _download(url, part, opener)
            size, digest = part.stat().st_size, _sha256(part)
            if (size, digest) == (resource.bytes, resource.sha256):
                part.replace(dest)
                return {
                    "name": resource.name, "url": url, "bytes": size, "sha256": digest,
                    "retrieved": datetime.now(UTC).isoformat(timespec="seconds"),
                }  # fmt: skip
            failures.append(f"{url}: {size} bytes, sha256 {digest}")
        except OSError as error:
            failures.append(f"{url}: {error}")
        finally:
            part.unlink(missing_ok=True)
    expected = f"{resource.bytes} bytes, sha256 {resource.sha256}"
    raise FetchError(f"{resource.name}: expected {expected}; got\n  " + "\n  ".join(failures))


def run(
    settings: runtime.Settings,
    cards: tuple[ParameterCard, ...],
    command: str,
    dataset_file: Path | None = None,
) -> list[Path]:
    """Write the versions table and the environment report under ``outputs/00_setup``.

    Args:
        settings: The resolved runtime.
        cards: Parameter cards for every option of this command.
        command: The command line, as the reader should retype it.
        dataset_file: A dataset file whose resources to fetch; unset fetches nothing.

    Returns:
        The written files: the report, then the versions table files.

    Raises:
        FetchError: Before any output is written, if a resource cannot be fetched.
    """
    fetched = fetch(settings.root, load(dataset_file).resources) if dataset_file else []
    day = date.today()
    number = dict(MODULES)[MODULE]
    out = settings.root / "outputs" / f"{number:02d}_{MODULE}"
    installed = runtime.versions()
    table = pl.DataFrame({"package": list(installed), "version": list(installed.values())})
    files = write_table(
        table, out, stem(MODULE, "versions", 1, day=day), columns=VERSION_COLUMNS, data_releases={}
    )
    facts = _facts(settings) + tuple(
        f"Resource {e['name']}: {e['bytes']:,} bytes, sha256 {e['sha256']}, from {e['url']} "
        f"({e['retrieved']}); recorded in {MANIFEST}."
        for e in fetched
    )
    # Measured before the page is rendered: the page cannot contain its own rendering time.
    usage = runtime.measure(settings.started)
    budget = runtime.budget(MODULE)
    page = ModulePage(
        title="M0 Setup: this machine and this environment",
        question=QUESTION,
        concept=CONCEPT,
        parameters=cards,
        results=facts,
        practice=PRACTICE,
        think=THINK,
        mistakes=MISTAKES,
        reproduce=Reproduction(
            command=command,
            wall_seconds=usage.wall_seconds,
            peak_rss_self=usage.peak_rss_self,
            peak_rss_children=usage.peak_rss_children,
            profile=f"{settings.profile}, {settings.threads} threads",
            budget=f"Track B: under {budget.minutes} min and {budget.ram_gb} GB",
            versions=installed,
        ),
    )
    report = write(page, out / f"{stem(MODULE, 'environment', 1, day=day)}.html")
    return [report, *files]
