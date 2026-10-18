"""From a PRIDE project or an own export to the files ``[proteo.table]`` reads (guide, M7).

``python -m ctx.proteo.prepare files <PXD>`` lists a PRIDE project's files with their bytes, so a
download is chosen with its size in view. ``fetch <PXD> <file> --to <folder>`` downloads one over
HTTPS, checks PRIDE's size and SHA-1, and prints the ``[[resources]]`` entry with the SHA-256
measured, for ``ctx setup`` to register in place. ``samples --format <f> --table <file>`` writes
a starter sample table: one row per sample column the adapter reads, the RNA sample and batch
left for the student to fill.
"""

import argparse
import csv
import hashlib
import json
import sys
import urllib.request
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

from ctx.proteo import ProteoError, adapters, quantities

# PRIDE Archive's REST API, version 3: every file of a project, with its size and SHA-1.
PRIDE_FILES = "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{accession}/files/all"
# PRIDE serves its FTP tree over HTTPS too; HTTPS passes where FTP is blocked.
FTP, HTTPS = "ftp://ftp.pride.ebi.ac.uk/", "https://ftp.pride.ebi.ac.uk/"
TIMEOUT = 120
CHUNK = 1 << 20
COLUMNS = ("protein_sample", "rna_sample", "batch")


def listing(
    accession: str, opener: Callable[..., Any] = urllib.request.urlopen
) -> list[dict[str, Any]]:
    """Each file of a PRIDE project: name, bytes, SHA-1 and its HTTPS URL, by name.

    Raises:
        OSError: If PRIDE does not answer with a list.
    """
    with opener(PRIDE_FILES.format(accession=accession), timeout=TIMEOUT) as response:
        found = json.loads(response.read())
    if not isinstance(found, list):
        raise OSError(f"{accession}: PRIDE answered {str(found)[:200]}")
    out = []
    for entry in found:
        ftp = [
            loc["value"]
            for loc in entry.get("publicFileLocations", [])
            if str(loc.get("value", "")).startswith(FTP)
        ]
        out.append(
            {
                "name": entry["fileName"],
                "bytes": int(entry["fileSizeBytes"]),
                "sha1": entry.get("checksum", ""),
                "url": ftp[0].replace(FTP, HTTPS, 1) if ftp else "",
            }
        )
    return sorted(out, key=lambda e: e["name"])


def fetch(
    entry: dict[str, Any], folder: Path, opener: Callable[..., Any] = urllib.request.urlopen
) -> Path:
    """Download one listed file into ``folder`` unless it is there, checked by size and SHA-1.

    The download goes to ``<name>.part`` and is renamed only once both checks pass.

    Raises:
        OSError: If the file has no HTTPS URL or fails a check; no partial file is left.
    """
    dest: Path = folder / str(entry["name"])
    if dest.is_file() and dest.stat().st_size == entry["bytes"]:
        return dest
    if not entry["url"]:
        raise OSError(f"{entry['name']}: PRIDE lists no FTP location")
    folder.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    sha1 = hashlib.sha1()  # noqa: S324 - PRIDE's own checksum, compared, not trusted for security
    try:
        with opener(entry["url"], timeout=TIMEOUT) as response, part.open("wb") as out:
            while block := response.read(CHUNK):
                sha1.update(block)
                out.write(block)
        size = part.stat().st_size
        if size != entry["bytes"] or (entry["sha1"] and sha1.hexdigest() != entry["sha1"]):
            raise OSError(
                f"{entry['name']}: {size} bytes, SHA-1 {sha1.hexdigest()}; PRIDE "
                f"lists {entry['bytes']} and {entry['sha1']}"
            )
        part.replace(dest)
    finally:
        part.unlink(missing_ok=True)
    return dest


def resource_block(name: str, url: str, path: Path, dest: str) -> str:
    """A ``[[resources]]`` entry for a file in place: its URL, bytes and SHA-256 measured."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(CHUNK):
            digest.update(block)
    return (
        f'[[resources]]\nname = "{name}"\nurls = ["{url}"]\nbytes = {path.stat().st_size}\n'
        f'sha256 = "{digest.hexdigest()}"\ndest = "{dest}"\n'
    )


def sample_columns(fmt: str, path: Path, member: str | None, quantity: str) -> list[str]:
    """The sample columns an adapter reads from an export, in the file's order.

    Raises:
        ProteoError: If the format is not a tool's own export.
    """
    data = quantities.raw_bytes(path, member)
    if fmt == "maxquant":
        names = adapters.header(data)
        return [c.removeprefix(f"{quantity} ") for c in names if c.startswith(f"{quantity} ")]
    if fmt == "diann":
        return [
            adapters.run_name(c)
            for c in adapters.header(data)
            if c not in adapters.DIANN_DESCRIPTORS
        ]
    raise ProteoError(f"--format {fmt!r}: a plain table names its own sample columns")


def parser() -> argparse.ArgumentParser:
    """The three subcommands, as the guide's M7 part types them."""
    top = argparse.ArgumentParser(prog="python -m ctx.proteo.prepare", description=__doc__)
    sub = top.add_subparsers(dest="command", required=True)
    files = sub.add_parser("files", help="list a PRIDE project's files and bytes")
    files.add_argument("accession")
    files.add_argument("--match", default="", help="keep names containing this text")
    get = sub.add_parser("fetch", help="download one file and print its [[resources]] entry")
    get.add_argument("accession")
    get.add_argument("file")
    get.add_argument("--root", type=Path, required=True, help="the ctx root")
    get.add_argument("--slug", required=True, help="the dataset's slug")
    get.add_argument("--name", default="protein_table", help="the resource's name")
    table = sub.add_parser("samples", help="write a starter sample table")
    table.add_argument("--format", choices=("maxquant", "diann"), required=True)
    table.add_argument("--table", type=Path, required=True, help="the export or its archive")
    table.add_argument("--member", default=None, help="the export's name inside an archive")
    table.add_argument("--quantity", default="LFQ intensity", help="MaxQuant's column prefix")
    table.add_argument("--out", type=Path, required=True, help="the sample table to write")
    return top


def main(argv: Sequence[str] | None = None) -> int:
    """Run one subcommand; a problem prints and returns 1."""
    args = parser().parse_args(argv)
    try:
        if args.command == "files":
            for e in listing(args.accession):
                if args.match in e["name"]:
                    print(f"{e['bytes']:>15,}  {e['name']}")
        elif args.command == "fetch":
            entry = next((e for e in listing(args.accession) if e["name"] == args.file), None)
            if entry is None:
                raise OSError(f"{args.accession} lists no file {args.file}")
            dest = f"data/raw/{args.slug}/{args.file}"
            path = fetch(entry, args.root / Path(dest).parent)
            print(resource_block(args.name, entry["url"], path, dest), end="")
        else:
            if args.out.exists():
                raise OSError(f"{args.out} exists; a starter never overwrites a table")
            names = sample_columns(args.format, args.table, args.member, args.quantity)
            with args.out.open("w", newline="") as handle:
                csv.writer(handle).writerows([COLUMNS, *([n, "", ""] for n in names)])
            print(f"{args.out}: {len(names)} samples; fill rna_sample and batch")
    except (OSError, ProteoError) as error:
        print(error, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
