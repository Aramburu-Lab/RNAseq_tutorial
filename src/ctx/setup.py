"""M0, environment report: what this machine and environment are, and what a run costs.

Downloads nothing (Phase 1). Every fact on the page is measured in this run, not assumed.
"""

import functools
import hashlib
import json
import os
import platform
import ssl
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
from ctx.proteo import pdc
from ctx.report import text
from ctx.report.builder import (
    KeyNumber,
    ModulePage,
    ParameterCard,
    Reproduction,
    ResultTable,
    write,
)
from ctx.viz.style import DATA_COLOURS, FONT_CHAIN, PROJECTOR_UNSAFE, THEMES, contrast

MODULE = "setup"
# WCAG 2.1 success criterion 1.4.11: non-text contrast for graphics.
GRAPHICS_CONTRAST = 3.0
MINIMUM = f"{GRAPHICS_CONTRAST:g}:1"
PAGE = MODULE
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
    return text.word(PAGE, "machine.font_missing")


def _memory_shown() -> str:
    """Physical memory as the page prints it, or that the OS did not report it."""
    memory = _memory_bytes()
    if memory:
        return text.fill(PAGE, "keys.memory_value", gib=memory / 2**30)
    return text.word(PAGE, "keys.memory_missing")


def _keys(settings: runtime.Settings) -> tuple[KeyNumber, ...]:
    """The four numbers a reader checks first, measured now."""
    word = functools.partial(text.word, PAGE)
    return (
        KeyNumber(word("keys.cores"), f"{os.cpu_count()}", word("keys.cores_detail")),
        KeyNumber(word("keys.memory"), _memory_shown()),
        KeyNumber(
            word("keys.threads"),
            f"{settings.threads}",
            text.fill(PAGE, "keys.threads_detail", profile=settings.profile),
        ),
        KeyNumber(word("keys.python"), platform.python_version()),
    )


def _colour_use(colour: str, ratio: float) -> str:
    """What a data colour's contrast on white allows, against WCAG's 3:1 for graphics."""
    if colour in PROJECTOR_UNSAFE:
        return text.word(PAGE, "colours.unsafe")
    if ratio >= GRAPHICS_CONTRAST:
        return text.fill(PAGE, "colours.meets", minimum=MINIMUM)
    return text.fill(PAGE, "colours.below", minimum=MINIMUM)


def _ratios() -> dict[str, float]:
    """Each data colour's contrast on the white theme's background."""
    return {c: contrast(c, THEMES["white"].background) for c in DATA_COLOURS}


def _skeleton() -> dict[str, Any]:
    """The page's fixed sections, Method and Interpretation filled from this run's numbers."""
    track_a = runtime.BUDGETS["A"]
    interpretation = text.fill(
        PAGE, "interpretation", minutes=track_a.minutes, ram_gb=track_a.ram_gb,
        laptop_cores=runtime.LAPTOP_CORES, cores=os.cpu_count(), memory=_memory_shown(),
        n_below=sum(r < GRAPHICS_CONTRAST for r in _ratios().values()),
        n_colours=len(DATA_COLOURS), minimum=MINIMUM,
    )  # fmt: skip
    return text.skeleton(PAGE) | {
        "method": text.fill(PAGE, "method", laptop_cores=runtime.LAPTOP_CORES),
        "interpretation": interpretation,
    }


def _facts(settings: runtime.Settings) -> tuple[ResultTable, ...]:
    """The measured state of this machine and run, as tables."""
    numpy.zeros(1)  # loads BLAS, so its thread pool can be read
    blas = [p["num_threads"] for p in threadpoolctl.threadpool_info() if p["user_api"] == "blas"]
    word = functools.partial(text.word, PAGE)
    system = text.fill(
        PAGE,
        "machine.system_value",
        system=platform.system(),
        release=platform.release(),
        machine=platform.machine(),
    )
    pools = text.fill(PAGE, "machine.thread_pools_value", polars=pl.thread_pool_size(), blas=blas)
    font = text.fill(PAGE, "machine.font_value", font=_resolved_font(), chain=", ".join(FONT_CHAIN))
    machine = ResultTable(
        word("machine.title"),
        text.words(PAGE, "machine.columns"),
        (
            (word("machine.system"), system),
            (word("machine.thread_pools"), pools),
            (word("machine.root"), runtime.shown_path(settings.root)),
            (word("machine.font"), font),
        ),
    )  # fmt: skip
    ratios = _ratios()
    colours = ResultTable(
        word("colours.title"),
        text.words(PAGE, "colours.columns"),
        tuple(
            (c, text.fill(PAGE, "colours.ratio", ratio=r), _colour_use(c, r))
            for c, r in ratios.items()
        ),
        numeric=(1,),
        swatch=0,
    )
    return machine, colours


def _resources(fetched: list[dict[str, Any]]) -> tuple[ResultTable, ...]:
    """What setup fetched, each file with its checksum; nothing when it fetched nothing."""
    if not fetched:
        return ()
    rows = tuple(
        (e["name"], f"{e['bytes']:,}", text.fill(PAGE, "resources.checksum", sha256=e["sha256"]),
         e["url"], e["retrieved"])
        for e in fetched
    )  # fmt: skip
    columns = text.words(PAGE, "resources.columns")
    title = text.fill(PAGE, "resources.title", manifest=MANIFEST)
    return (ResultTable(title, columns, rows, numeric=(1,)),)


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


def _md5(path: Path) -> str:
    """MD5 of a file, read in chunks: the checksum PDC publishes."""
    digest = hashlib.md5()  # noqa: S324 - compared with the publisher's MD5, not for security
    with path.open("rb") as handle:
        while block := handle.read(CHUNK):
            digest.update(block)
    return digest.hexdigest()


def tls_context(ciphers: str) -> ssl.SSLContext:
    """Python's default verifying context with ``ciphers`` added to its TLS 1.2 list.

    The default list holds only forward-secret key exchange; a server offering only RSA key
    exchange (miRTarBase's, 2026-10-08) fails the handshake on it. Appending, not replacing,
    keeps every default cipher preferred, and verification is the default context's.

    Args:
        ciphers: An OpenSSL cipher string, as the resource's ``tls_ciphers`` declares it.

    Returns:
        The context to open that resource's URLs with.
    """
    context = ssl.create_default_context()
    default = [c["name"] for c in context.get_ciphers() if c["protocol"] != "TLSv1.3"]
    context.set_ciphers(":".join([*default, ciphers]))
    return context


def _download(url: str, part: Path, opener: Callable[..., Any], ciphers: str = "") -> None:
    """Stream one URL to ``part``, with the resource's own TLS ciphers if it declares any."""
    extra = {"context": tls_context(ciphers)} if ciphers else {}
    with opener(url, timeout=120, **extra) as response, part.open("wb") as out:
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
    not fetched again: one placed there by hand, as GTEx V11's 24 GB were, is registered with
    its hash measured on arrival and ``found`` true, never downloaded. The destination is
    checked against the root by its path, not its link target, so an input folder linked into
    a scratch root counts as inside it. Every fetched or confirmed file is recorded in
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
        base = Path(os.path.normpath(root.absolute()))
        dest = Path(os.path.normpath(base / resource.dest))
        if not dest.is_relative_to(base):
            raise FetchError(f"{resource.name}: {resource.dest} is outside {root}")
        entry = manifest.get(resource.dest)
        present = dest.is_file() and _sha256(dest) == resource.sha256
        if present and not entry:
            entry = _found(resource, dest)
        elif not (present and entry):
            entry = _fetch_one(resource, dest, opener)
        manifest[resource.dest] = entry
        entries.append(entry)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return entries


def _found(resource: Resource, dest: Path) -> dict[str, Any]:
    """The manifest entry of a file already in place with its declared size and hash."""
    return {"name": resource.name, "url": resource.urls[0], "bytes": dest.stat().st_size,
            "sha256": resource.sha256, "found": True,
            "retrieved": datetime.now(UTC).isoformat(timespec="seconds")}  # fmt: skip


def _fetch_one(resource: Resource, dest: Path, opener: Callable[..., Any]) -> dict[str, Any]:
    """Try each mirror until one gives the declared size and hash; leave no partial file."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    failures = []
    for url in resource.urls:
        try:
            # A PDC URL names a file whose signed URL expires: resolved now, MD5 checked first.
            signed = None
            if url.startswith(pdc.SCHEME):
                print(text.fill(PAGE, "resources.dua", name=resource.name, url=pdc.DUA_URL))
                signed = pdc.resolve(url, opener)
            _download(signed.url if signed else url, part, opener, resource.tls_ciphers)
            size, digest = part.stat().st_size, _sha256(part)
            if signed and _md5(part) != signed.md5:
                failures.append(f"{url}: MD5 {_md5(part)}, PDC lists {signed.md5}")
                continue
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
    facts = _facts(settings) + _resources(fetched)
    # Measured before the page is rendered: the page cannot contain its own rendering time.
    usage = runtime.measure(settings.started)
    budget = runtime.budget(MODULE)
    page = ModulePage(
        title=text.word(PAGE, "title"),
        **_skeleton(),
        parameters=cards,
        results=facts,
        key_numbers=_keys(settings),
        reproduce=Reproduction(
            command=command,
            wall_seconds=usage.wall_seconds,
            peak_rss_self=usage.peak_rss_self,
            peak_rss_children=usage.peak_rss_children,
            profile=f"{settings.profile}, {settings.threads} threads",
            budget=text.fill(PAGE, "budget", minutes=budget.minutes, ram_gb=budget.ram_gb),
            versions=installed,
        ),
    )
    report = write(page, out / f"{stem(MODULE, 'environment', 1, day=day)}.html")
    return [report, *files]
