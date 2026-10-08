"""Per-gene API responses, cached under ``resources/cache/`` by ``ctx genecard`` alone.

One response is ``resources/cache/<provider>/<release>/<ENSG>/<endpoint>.gz``, its body as
received, gzipped. ``resources/cache/manifest.json`` records each one's URL, UTC retrieval
time, HTTP status, bytes and the SHA-256 of the body measured on arrival; ``ctx setup``'s own
``resources/manifest.json`` is never touched (TODO.md M9). A cached body is checked against
its hash on every read. A refetch whose body differs under the same key stops the run;
``--refresh`` first moves the gene's folder to ``resources/cache/superseded/<date>/``.
A failed call is never cached: the card shows that panel as unavailable, naming the error.
"""

import gzip
import hashlib
import json
import shutil
import urllib.error
import urllib.request
from collections.abc import Callable
from dataclasses import dataclass, field
from datetime import UTC, date, datetime
from pathlib import Path
from typing import Any

from ctx.genecard import GenecardError

CACHE = Path("resources") / "cache"
MANIFEST = CACHE / "manifest.json"
SUPERSEDED = CACHE / "superseded"
TIMEOUT = 60
# HPA's and Ensembl's servers answer urllib's default agent; a named one says who asks.
AGENT = "ctx-genecard (Cancer Transcriptomics Teaching Suite)"


@dataclass(frozen=True)
class Response:
    """One cached or fetched response; ``body`` is None when the call failed or was not made."""

    url: str
    body: bytes | None
    retrieved: str = ""
    error: str = ""


def _open(url: str) -> Any:
    """Urllib with the suite's agent and timeout: the default opener."""
    request = urllib.request.Request(url, headers={"User-Agent": AGENT})  # noqa: S310 - https
    return urllib.request.urlopen(request, timeout=TIMEOUT)  # noqa: S310 - https URLs only


@dataclass
class Cache:
    """The response cache of one root.

    Attributes:
        root: The ctx root.
        offline: Read the cache only; a response not there is unavailable, never fetched.
        refresh: Move each gene's cached folder aside before its first fetch of this run.
        opener: How to open a URL; tests pass their own.
        today: The date that names a superseded folder.
    """

    root: Path
    offline: bool = False
    refresh: bool = False
    opener: Callable[[str], Any] = _open
    today: date = field(default_factory=date.today)
    _moved: set[Path] = field(default_factory=set)

    def _manifest(self) -> dict[str, Any]:
        """The cache manifest, empty before the first response."""
        path = self.root / MANIFEST
        return json.loads(path.read_text()) if path.is_file() else {}

    def _record(self, key: str, entry: dict[str, Any] | None) -> None:
        """Write one entry, or drop it, and save the manifest at once: a crash loses none."""
        manifest = self._manifest()
        if entry is None:
            manifest.pop(key, None)
        else:
            manifest[key] = entry
        path = self.root / MANIFEST
        path.parent.mkdir(parents=True, exist_ok=True)
        part = path.with_name(path.name + ".part")
        part.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        part.replace(path)

    def _supersede(self, folder: Path) -> None:
        """Move one gene's folder to ``superseded/<date>/``, its manifest entries with it."""
        if folder in self._moved or not folder.is_dir():
            self._moved.add(folder)
            return
        relative = folder.relative_to(self.root / CACHE)
        target = self.root / SUPERSEDED / self.today.isoformat() / relative
        if target.exists():
            raise GenecardError(f"{target} exists: --refresh twice in one day moves nothing")
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(folder, target)
        for key in [k for k in self._manifest() if k.startswith(f"{relative.as_posix()}/")]:
            self._record(key, None)
        self._moved.add(folder)

    def get(self, provider: str, release: str, gene: str, endpoint: str, url: str) -> Response:
        """The response for one gene and endpoint: from the cache, else fetched and cached.

        Raises:
            GenecardError: If a cached body no longer matches its recorded SHA-256, or a
                refetch gives another body under a key the manifest already holds.
        """
        key = f"{provider}/{release}/{gene}/{endpoint}.gz"
        path = self.root / CACHE / key
        if self.refresh and not self.offline:
            self._supersede(path.parent)
        entry = self._manifest().get(key)
        if path.is_file() and entry:
            body = gzip.decompress(path.read_bytes())
            if hashlib.sha256(body).hexdigest() != entry["sha256"]:
                raise GenecardError(f"{path}: body differs from its SHA-256 in {MANIFEST}")
            return Response(url, body, entry["retrieved"])
        if self.offline:
            return Response(url, None, error="not cached, and --offline makes no request")
        return self._fetch(key, path, url, entry)

    def _fetch(self, key: str, path: Path, url: str, entry: dict[str, Any] | None) -> Response:
        """Fetch one URL, check it against an earlier record of the key, and cache it."""
        try:
            with self.opener(url) as answer:
                body, status = answer.read(), getattr(answer, "status", 200)
        except (urllib.error.URLError, OSError) as error:
            return Response(url, None, error=str(error))
        digest = hashlib.sha256(body).hexdigest()
        if entry and entry["sha256"] != digest:
            when = entry["retrieved"]
            raise GenecardError(f"{key}: the provider now sends another body than the one "
                                f"cached {when}; --refresh keeps the old aside")  # fmt: skip
        retrieved = datetime.now(UTC).isoformat(timespec="seconds")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(gzip.compress(body, mtime=0))
        self._record(key, {"url": url, "retrieved": retrieved, "status": status,
                           "bytes": len(body), "sha256": digest})  # fmt: skip
        return Response(url, body, retrieved)
