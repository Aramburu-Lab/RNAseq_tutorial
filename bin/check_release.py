#!/usr/bin/env python3
# /// codebase-meta
# exempt: Template tooling (D-020) that checks a release tag and writes nothing.
# ///
"""Refuse a release tag that the licence and citation metadata do not back up (D-020).

A tag is how every consumer pins this code, and CITATION.cff is what GitHub and Zenodo show
when someone cites it. They drift the moment a tag is cut without editing the file, so CI
runs this on every `v*` tag and `pixi run release-check vX.Y.Z` runs it before pushing one.

Checks: LICENSE exists; the tag is vMAJOR.MINOR.PATCH; CITATION.cff carries that version, a
date-released and the repository origin really points at; conf/references.sha256 names its
references or says none; and once an earlier tag's CITATION.cff carried a doi, this one does.
Earlier tags alone never demand one: riboqc had six and no Zenodo deposit, so a rule keyed on
tags made its next release impossible.

Usage: python3 bin/check_release.py [tag] [project-root]   (tag defaults to HEAD's exact tag)
"""

import re
import subprocess
import sys
from pathlib import Path

SEMVER = re.compile(r"^v(\d+\.\d+\.\d+)$")


def git(root: Path, *args: str) -> str:
    """Run git in the project and return stdout, empty on failure."""
    r = subprocess.run(["git", *args], cwd=root, capture_output=True, text=True, check=False)
    return r.stdout.strip() if r.returncode == 0 else ""


def field(cff: str, key: str) -> str:
    """A top-level scalar from CITATION.cff, quotes stripped; empty when absent."""
    hit = re.search(rf"^{key}:\s*(.+?)\s*(?:\s#.*)?$", cff, re.M)
    return hit.group(1).strip("'\"") if hit else ""


def problems(root: Path, tag: str) -> list[str]:
    """Everything wrong with releasing this tree as `tag`."""
    out = []
    if not (root / "LICENSE").is_file():
        out.append("no LICENSE file: nobody may legally reuse the code")
    semver = SEMVER.match(tag)
    if not semver:
        return [*out, f"tag {tag!r} is not vMAJOR.MINOR.PATCH"]
    cff_path = root / "CITATION.cff"
    if not cff_path.is_file():
        return [*out, "no CITATION.cff: the release cannot be cited"]
    cff = cff_path.read_text(encoding="utf-8")
    if field(cff, "version") != semver.group(1):
        out.append(
            f"CITATION.cff version is {field(cff, 'version')!r}, the tag says {semver.group(1)}"
        )
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", field(cff, "date-released")):
        out.append("CITATION.cff has no date-released: YYYY-MM-DD")
    refs = root / "conf" / "references.sha256"
    text = refs.read_text(encoding="utf-8") if refs.is_file() else ""
    if not re.search(r"^([0-9a-f]{64}  \S|none: \S)", text, re.M):
        out.append("conf/references.sha256 lists no checksum and no 'none: <why>' line")
    return out + citation_links(root, tag, cff)


def bare(url: str) -> str:
    """host/owner/repo of a git or https URL, lowercased, so two spellings compare equal."""
    url = re.sub(r"^(git@|ssh://git@|https?://)", "", url.strip()).replace(":", "/", 1)
    return re.sub(r"\.git$", "", url).rstrip("/").lower()


def citation_links(root: Path, tag: str, cff: str) -> list[str]:
    """repository-code must be the real origin; a DOI once recorded must not be dropped."""
    out = []
    origin = git(root, "remote", "get-url", "origin")
    if origin and bare(field(cff, "repository-code")) != bare(origin):
        out.append(
            f"CITATION.cff repository-code is {field(cff, 'repository-code')!r} but origin is "
            f"{origin!r}: one of them is wrong"
        )
    earlier = [t for t in git(root, "tag", "--list", "v*").split() if SEMVER.match(t) and t != tag]
    cited = [t for t in earlier if field(git(root, "show", f"{t}:CITATION.cff"), "doi")]
    if cited and not field(cff, "doi").startswith("10."):
        out.append(f"CITATION.cff has no doi, but {cited[-1]} carried one: keep the concept DOI")
    return out


def main() -> int:
    """Print what blocks the release; exit 1 if anything does."""
    args = sys.argv[1:]
    root = Path(args[1]).resolve() if len(args) > 1 else Path(__file__).resolve().parents[1]
    tag = args[0] if args else git(root, "describe", "--tags", "--exact-match")
    if not tag:
        print("release-check: HEAD carries no tag; pass one, e.g. pixi run release-check v0.1.0")
        return 1
    bad = problems(root, tag)
    for p in bad:
        print(f"release-check {tag}: {p}")
    if not bad:
        print(f"release-check {tag}: LICENSE, CITATION.cff and tag agree")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
