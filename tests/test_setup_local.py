"""ctx setup's two R3.4 routes: a resource copied from a local folder, and a zip member kept."""

import hashlib
import io
import json
import os
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

from ctx.dataset import DatasetError, Resource, load
from ctx.report.builder import unpacked
from ctx.setup import MANIFEST, FetchError, fetch

CTX = Path(sys.executable).parent / "ctx"
BODY = b"Gene_id\tGene_name\tTMT1\n"
DIGEST = hashlib.sha256(BODY).hexdigest()
# Port 9 (discard) on loopback: any attempt to download fails, so a pass proves no download.
UNREACHABLE = "http://127.0.0.1:9/supplementaryFiles"


def resource(**extra: str) -> Resource:
    return Resource("table", (UNREACHABLE,), len(BODY), DIGEST, "data/raw/s2.xlsx", **extra)


def no_network(url: str, **_: object) -> io.BytesIO:
    raise AssertionError(f"opened {url}")


def test_copied_from_folder_never_moved(tmp_path):
    # Mutations: the file moved out of the folder; the manifest silent on where it came from.
    folder = tmp_path / "copies" / "deeper"
    folder.mkdir(parents=True)
    source = folder / "renamed.bin"
    source.write_bytes(BODY)
    (entry,) = fetch(tmp_path / "root", [resource()], no_network, from_dir=tmp_path / "copies")
    assert (tmp_path / "root" / "data" / "raw" / "s2.xlsx").read_bytes() == BODY
    assert source.read_bytes() == BODY, "the source must stay where it was"
    manifest = json.loads((tmp_path / "root" / MANIFEST).read_text())
    assert manifest["data/raw/s2.xlsx"]["copied_from"] == str(source.resolve())
    assert entry["sha256"] == DIGEST and entry["url"] == UNREACHABLE


def test_same_named_mismatch_refused(tmp_path):
    # Mutation: a damaged copy passed over, so setup downloads, or kept.
    folder = tmp_path / "copies"
    folder.mkdir()
    (folder / "s2.xlsx").write_bytes(BODY.replace(b"TMT", b"XXX"))
    with pytest.raises(FetchError, match="s2.xlsx has"):
        fetch(tmp_path / "root", [resource()], no_network, from_dir=folder)
    assert not (tmp_path / "root" / "data" / "raw" / "s2.xlsx").exists()


def test_folder_without_the_file_downloads(tmp_path):
    folder = tmp_path / "copies"
    folder.mkdir()
    (folder / "other.txt").write_bytes(b"unrelated")
    with pytest.raises(AssertionError, match="opened"):
        fetch(tmp_path / "root", [resource()], no_network, from_dir=folder)


def zipped(member: str, body: bytes) -> bytes:
    # A fixed entry time: a zip stores when it was written, so cases built at collection
    # differed between xdist workers, which then refused to run (2026-10-21).
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as bundle:
        for name, data in (("other.gif", b"GIF89a"), (member, body)):
            bundle.writestr(zipfile.ZipInfo(name, date_time=(2026, 10, 21, 0, 0, 0)), data)
    return buffer.getvalue()


def test_member_kept_and_checked(tmp_path):
    # Mutation: the archive kept at dest, or checked against the member's hash whole.
    archive = zipped("s2.xlsx", BODY)
    (entry,) = fetch(
        tmp_path / "root", [resource(member="s2.xlsx")], lambda *a, **k: io.BytesIO(archive)
    )
    assert (tmp_path / "root" / "data" / "raw" / "s2.xlsx").read_bytes() == BODY
    assert entry["member"] == "s2.xlsx"
    assert entry["archive_sha256"] == hashlib.sha256(archive).hexdigest()


@pytest.mark.parametrize(
    ("archive", "match"),
    [(zipped("elsewhere.xlsx", BODY), "no member"), (BODY, "no member"),
     (zipped("s2.xlsx", BODY + b"x"), "sha256")],
    ids=["member missing", "not a zip", "member's hash wrong"],
)  # fmt: skip
def test_member_refused(tmp_path, archive, match):
    with pytest.raises(FetchError, match=match):
        fetch(tmp_path / "root", [resource(member="s2.xlsx")], lambda *a, **k: io.BytesIO(archive))
    assert not list((tmp_path / "root" / "data" / "raw").iterdir()), "no partial file may remain"


def spec_text(**extra: object) -> str:
    lines = "".join(f"{k} = {json.dumps(v)}\n" for k, v in extra.items())
    return (
        f'[dataset]\nslug = "local"\n\n[[resources]]\nname = "table"\nurls = ["{UNREACHABLE}"]\n'
        f'bytes = {len(BODY)}\nsha256 = "{DIGEST}"\ndest = "data/raw/s2.xlsx"\n{lines}\n'
        "[guard]\nids = []\npatterns = []\n"
    )


def test_loader_reads_member(tmp_path):
    spec = tmp_path / "local.toml"
    spec.write_text(spec_text(member="s2.xlsx", archive_sha256="0" * 64))
    (loaded,) = load(spec).resources
    assert (loaded.member, loaded.archive_sha256) == ("s2.xlsx", "0" * 64)


@pytest.mark.parametrize(
    ("extra", "match"),
    [
        ({"archive_sha256": "0" * 64}, "archive of a member"),
        ({"member": 3}, "member: expected str"),
    ],
)
def test_loader_refuses(tmp_path, extra, match):
    spec = tmp_path / "local.toml"
    spec.write_text(spec_text(**extra))
    with pytest.raises(DatasetError, match=match):
        load(spec)


def run_setup(root: Path, *args: str) -> subprocess.CompletedProcess[str]:
    env = {k: v for k, v in os.environ.items() if not k.startswith(("SLURM_", "CTX_"))}
    return subprocess.run(
        [str(CTX), "--root", str(root), "setup", *args],
        cwd=root, env=env, capture_output=True, text=True, check=False,
    )  # fmt: skip


def test_cli_from_dir(tmp_path):
    spec = tmp_path / "local.toml"
    spec.write_text(spec_text())
    folder = tmp_path / "copies"
    folder.mkdir()
    (folder / "copy.xlsx").write_bytes(BODY)
    root = tmp_path / "root"
    root.mkdir()
    done = run_setup(root, "--dataset-file", str(spec), "--from-dir", str(folder))
    assert done.returncode == 0, done.stderr
    assert (root / "data" / "raw" / "s2.xlsx").read_bytes() == BODY
    (page,) = (root / "outputs" / "00_setup").glob("*.html")
    assert "--from-dir" in unpacked(page.read_text())
    alone = run_setup(root, "--from-dir", str(folder))
    assert alone.returncode == 1 and "give --dataset-file" in alone.stderr
