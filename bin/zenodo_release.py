#!/usr/bin/env python3
# /// codebase-meta
# exempt: Template tooling (D-020), run by hand to archive a release on Zenodo.
# ///
"""Archive one release on Zenodo under one DOI: source, images, models, checksums (D-020).

Zenodo's GitHub integration archives the source and nothing else, and CI cannot see models
that live on HPC2N. So this runs where the files are, at release time, after the tag passes
`release-check`. It uploads `git archive` of the tag, every FILE given, conf/references.sha256
and a SHA256SUMS over all of them, with metadata read from CITATION.cff, the one place it lives.

The first release creates a record; later ones add a version to the concept DOI in
CITATION.cff. Without --publish the deposit stays a draft, published from the web page it
prints after review: publishing mints a DOI that can never be deleted.

Usage: ZENODO_TOKEN=... python3 bin/zenodo_release.py vX.Y.Z [FILE ...] [--publish] [--sandbox]
       --sandbox targets sandbox.zenodo.org, whose token and DOIs are separate and disposable.
"""

import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
from check_release import field, problems  # noqa: E402


def api(method: str, url: str, token: str, body=None, path: Path | None = None) -> dict:
    """One Zenodo call; JSON in and out, or a file streamed up. Raises on any HTTP error."""
    headers = {"Authorization": f"Bearer {token}"}
    data = None
    if path is not None:
        headers |= {
            "Content-Type": "application/octet-stream",
            "Content-Length": str(path.stat().st_size),
        }
        data = path.open("rb")
    elif body is not None:
        headers["Content-Type"] = "application/json"
        data = json.dumps(body).encode()
    with urllib.request.urlopen(urllib.request.Request(url, data, headers, method=method)) as r:
        raw = r.read()
    return json.loads(raw) if raw else {}


def metadata(cff: str, version: str) -> dict:
    """Zenodo metadata from CITATION.cff, so the two can never disagree."""
    creators = []
    for block in re.split(r"\n\s*-\s+(?=given-names|family-names|name)", cff.split("authors:")[1]):
        family, given = field_in(block, "family-names"), field_in(block, "given-names")
        if family or given:
            creator = {"name": f"{family}, {given}".strip(", ")}
            if orcid := field_in(block, "orcid"):
                creator["orcid"] = orcid.rsplit("/", 1)[-1]
            creators.append(creator)
    title = field(cff, "title")
    return {
        "upload_type": "software",
        "title": f"{title} {version}",
        "creators": creators,
        "description": field(cff, "abstract")
        or f"Source code, container images, models and checksums of {title} {version}.",
        "version": version,
        "license": field(cff, "license").lower(),
        "publication_date": field(cff, "date-released"),
    }


def field_in(block: str, key: str) -> str:
    """A key from one indented author entry."""
    hit = re.search(rf"^\s*-?\s*{key}:\s*(.+?)\s*$", block, re.M)
    return hit.group(1).strip("'\"") if hit else ""


def bundle(tag: str, files: list[Path], work: Path) -> list[Path]:
    """The source archive, the given files, reference checksums and a SHA256SUMS of all."""
    slug = ROOT.name
    source = work / f"{slug}-{tag[1:]}.tar.gz"
    subprocess.run(
        [
            "git",
            "archive",
            "--format=tar.gz",
            f"--prefix={slug}-{tag[1:]}/",
            "-o",
            str(source),
            tag,
        ],
        cwd=ROOT,
        check=True,
    )
    out = [source, *files]
    refs = ROOT / "conf" / "references.sha256"
    if refs.is_file():
        out.append(refs)
    sums = work / "SHA256SUMS"
    lines = []
    for p in out:
        digest = hashlib.sha256()
        with p.open("rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                digest.update(chunk)
        lines.append(f"{digest.hexdigest()}  {p.name}")
    sums.write_text("\n".join(lines) + "\n")
    return [*out, sums]


def draft(base: str, token: str, concept: str) -> dict:
    """A new deposit, or a new version of the concept with the old version's files removed."""
    if not concept:
        return api("POST", f"{base}/api/deposit/depositions", token, {})
    latest = api("GET", f"{base}/api/records/{concept}", token)["id"]
    made = api("POST", f"{base}/api/deposit/depositions/{latest}/actions/newversion", token)
    dep = api("GET", made["links"]["latest_draft"], token)
    for f in dep.get("files", []):
        api("DELETE", f"{base}/api/deposit/depositions/{dep['id']}/files/{f['id']}", token)
    return dep


def main() -> int:
    """Upload, set metadata, publish only when asked."""
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not args:
        print(__doc__)
        return 1
    tag, files = args[0], [Path(f).resolve() for f in args[1:]]
    if bad := problems(ROOT, tag):
        print("\n".join(f"release-check {tag}: {p}" for p in bad))
        return 1
    if missing := [str(f) for f in files if not f.is_file()]:
        print(f"zenodo: not found: {', '.join(missing)}")
        return 1
    token = os.environ.get("ZENODO_TOKEN", "")
    if not token:
        print("zenodo: set ZENODO_TOKEN (bin/secrets.env), never on the command line")
        return 1
    sandbox = "--sandbox" in sys.argv
    base = os.environ.get("ZENODO_URL") or f"https://{'sandbox.' if sandbox else ''}zenodo.org"
    cff = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    concept = field(cff, "doi").rsplit(".", 1)[-1] if field(cff, "doi") else ""
    with tempfile.TemporaryDirectory() as tmp:
        dep = draft(base, token, concept)
        for p in bundle(tag, files, Path(tmp)):
            print(f"zenodo: uploading {p.name} ({p.stat().st_size:,} bytes)")
            api("PUT", f"{dep['links']['bucket']}/{p.name}", token, path=p)
    api(
        "PUT",
        f"{base}/api/deposit/depositions/{dep['id']}",
        token,
        {"metadata": metadata(cff, tag[1:])},
    )
    if "--publish" not in sys.argv:
        # Re-running would open a second draft; the review and the Publish button are one page.
        print(f"zenodo: draft ready; review it and press Publish there: {dep['links']['html']}")
        return 0
    done = api("POST", f"{base}/api/deposit/depositions/{dep['id']}/actions/publish", token)
    print(f"zenodo: published {done['doi']}; concept DOI {done.get('conceptdoi', '')}")
    if not concept:
        print(f"zenodo: now add to CITATION.cff and commit:  doi: {done.get('conceptdoi')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
