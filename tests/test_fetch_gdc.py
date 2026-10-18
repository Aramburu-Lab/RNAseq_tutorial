"""R1: scripts/bash/fetch_gdc.sbatch against a made-up GDC API on this machine, read by gather.

The API serves the folders tests/gdc_fixture.py and tests/gdc_mirna_fixture.py build, so a
fetched folder passes ``ctx.gather.gdc``'s checks only if the script lays it out as gather reads
one. No request leaves 127.0.0.1. Each test names the mutation it catches.
"""

import hashlib
import json
import shutil
import threading
from collections.abc import Iterator
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import gdc_mirna_fixture as mirna_fixture
import pytest
from gdc_fixture import RELEASE, aliquot, build
from test_cluster_scripts import BASH, _run
from test_cluster_scripts import site as _cluster_site

from ctx.dataset import load
from ctx.gather import gdc, gdc_mirna

SCRIPT = "fetch_gdc.sbatch"
CDR_ID = "1b5f413e-a8d1-4d10-92eb-7c4ae739ed81"
CDR_NAME = "TCGA-CDR-SupplementalTableS1.xlsx"
CDR_BODY = b"a made-up workbook"
PANCAN = "PanCan-General_Open_GDC-Manifest_2.txt"
# The STAR aliquots: two cases with a tumour and a normal, one tumour alone.
STAR = (aliquot(1), aliquot(1, "11"), aliquot(2), aliquot(2, "11"))
KINDS = {"Gene Expression Quantification": "star", "miRNA Expression Quantification": "mirna",
         "Isoform Expression Quantification": "isoform"}  # fmt: skip


@pytest.fixture
def site(tmp_path):
    """test_cluster_scripts' site, its clone holding this script too."""
    root, env = _cluster_site.__wrapped__(tmp_path)
    shutil.copy2(BASH / SCRIPT, Path(env["CTX_REPO"]) / "scripts" / "bash" / SCRIPT)
    return root, env


class Api:
    """The folders a made-up API serves, and every request it answered."""

    def __init__(self, source: Path) -> None:
        """Serve the fixture folders under ``source``: star, mirna and isoform."""
        self.source = source
        self.requests: list[str] = []
        # File IDs served with a byte added, as a transfer that went wrong would deliver them.
        self.corrupt: set[str] = set()
        # Records whose pagination says GDC holds more hits than it returned.
        self.short: set[str] = set()
        self.files = {p.parent.name: p for p in source.rglob("*") if p.is_file()
                      and len(p.parent.name) == 36}  # fmt: skip
        self.files[CDR_ID] = source / CDR_NAME
        (source / CDR_NAME).write_bytes(CDR_BODY)

    def answer(self, method: str, path: str, body: bytes) -> bytes:  # noqa: PLR0911 - one return per endpoint
        """What GDC's endpoint would return: records with pagination, files by ID."""
        self.requests.append(f"{method} {path}")
        if path == "/status":
            return json.dumps({"data_release": RELEASE, "status": "OK"}).encode()
        if path.startswith("/data/"):
            file_id = path.removeprefix("/data/")
            return self.files[file_id].read_bytes() + (b"x" if file_id in self.corrupt else b"")
        if path.endswith(PANCAN):
            md5 = hashlib.md5(CDR_BODY).hexdigest()  # noqa: S324 - GDC's checksum
            row = f"{CDR_ID}\t{CDR_NAME}\t{md5}\t{len(CDR_BODY)}"
            return f"id\tfilename\tmd5\tsize\n{row}\n".encode()
        query = json.loads(body)
        if path in ("/cases", "/annotations"):
            name = gdc.CASES if path == "/cases" else gdc.ANNOTATIONS
            paged = self._paged(json.loads((self.source / "star" / name).read_text()))
            if path in self.short:
                record = json.loads(paged)
                record["data"]["pagination"]["total"] += 1
                return json.dumps(record).encode()
            return paged
        types = [c["content"]["value"][0] for c in query["filters"]["content"]
                 if c["content"]["field"] == "data_type"]  # fmt: skip
        if types == ["Aligned Reads"]:
            return self._paged({"data": {"hits": []}})
        kind = KINDS[types[0]]
        folder = self.source / kind
        if kind == "star":
            metadata, manifest = gdc.METADATA, gdc.MANIFEST
        else:
            metadata, manifest, _ = gdc_mirna.KINDS[kind]
        if "return_type" in query:
            # GDC ends a manifest's lines with LF and a TSV's with CRLF (both read 2026-10-18).
            return (folder / manifest).read_bytes().replace(b"\r\n", b"\n")
        return (folder / metadata).read_bytes()

    @staticmethod
    def _paged(record: dict) -> bytes:
        """A record with the pagination GDC's JSON carries, all hits on one page."""
        hits = record["data"]["hits"]
        record["data"]["pagination"] = {"total": len(hits), "count": len(hits)}
        record["warnings"] = {}
        return json.dumps(record).encode()


@pytest.fixture
def api(tmp_path) -> Iterator[tuple[Api, str]]:
    """The made-up API on a free port of 127.0.0.1, and its base URL."""
    source = tmp_path / "gdc_source"
    build(source / "star", STAR)
    small = tuple(mirna_fixture.mirna_aliquot(a) for a in STAR)
    for kind in (gdc_mirna.MIRNA, gdc_mirna.ISOFORM):
        mirna_fixture.build(source / kind, kind, small)
    served = Api(source)

    class Handler(BaseHTTPRequestHandler):
        """GET and POST, answered from ``served``."""

        def _reply(self, method: str) -> None:
            """Answer one request, or 404 for a path the API does not hold."""
            body = self.rfile.read(int(self.headers.get("Content-Length") or 0))
            try:
                data = served.answer(method, self.path, body)
            except (KeyError, OSError):
                self.send_response(404)
                self.end_headers()
                return
            self.send_response(200)
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)

        def do_GET(self) -> None:  # noqa: N802 - http.server's name
            """A GET."""
            self._reply("GET")

        def do_POST(self) -> None:  # noqa: N802 - http.server's name
            """A POST."""
            self._reply("POST")

        def log_message(self, *args: object) -> None:
            """Silent."""

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield served, f"http://127.0.0.1:{server.server_address[1]}"
    server.shutdown()


def _fetch(env, api_url, folder, cwd, **given):
    return _run([Path(env["CTX_REPO"]) / "scripts" / "bash" / SCRIPT, "FAKE-BRCA", folder], cwd,
                env, GDC_API=api_url, GDC_WEB=f"{api_url}/web", CTX_SITE="laptop",
                **given)  # fmt: skip


def _spec(tmp_path, folder):
    path = mirna_fixture.spec_file(tmp_path / "fake.toml")
    path.write_text(path.read_text().replace('normal = "Solid Tissue Normal"\n',
                                             f'normal = "Solid Tissue Normal"\ncdr_file = '
                                             f'"{CDR_NAME}"\ncdr_type = "BRCA"\n'))  # fmt: skip
    assert folder.is_dir()
    return load(path).gdc


def test_fetched_folders_pass_gathers_checks(site, api, tmp_path):
    """The three folders the script writes pass every check gather's readers make, the CDR
    beside the STAR metadata with its retrieval row, and the end check passes. Catches a file
    name, query field or log row gather reads that the script writes otherwise, and the CDR
    left in a file_id folder."""
    root, env = site
    served, url = api
    folder = root / "home" / "ctx_fake" / "data" / "raw" / "fake"
    done = _fetch(env, url, folder, root, CDR="1")
    assert done.returncode == 0, done.stdout[-3000:] + done.stderr[-2000:]
    assert done.stdout.rstrip().endswith("every check passed")
    spec = _spec(tmp_path, folder)
    star = gdc.read_folder(folder, spec)
    gdc.verify(star)
    assert star.files.height == len(STAR)
    for kind in (gdc_mirna.MIRNA, gdc_mirna.ISOFORM):
        gdc.verify(gdc_mirna.read_kind(Path(f"{folder}_{kind}"), spec, kind))
    assert (folder / CDR_NAME).read_bytes() == CDR_BODY
    assert gdc.library_selection(folder, spec) is not None
    logged = (folder / gdc.RETRIEVAL).read_text()
    assert f"\t{CDR_NAME}\tGET " in logged and f"\t{gdc.READ_GROUPS}\tPOST " in logged
    assert "PASTE pairs, cases with a primary tumour and a solid tissue normal: 2" in done.stdout


def test_size_first_then_resume_fetches_nothing_twice(site, api):
    """DATA=0 saves the metadata alone and ends on the bytes the manifests list; the download
    then fetches each data file once, and a rerun fetches none and keeps every record. Catches
    a data file fetched by the size run, a record fetched again, and a size off the manifests."""
    root, env = site
    served, url = api
    folder = root / "store" / "fake"
    sized = _fetch(env, url, folder, root, DATA="0", CDR="0")
    assert sized.returncode == 0, sized.stdout[-3000:]
    listed = sum(p.stat().st_size for k, p in served.files.items() if k != CDR_ID)
    assert sized.stdout.rstrip().splitlines()[-1].endswith(f": {listed} bytes")
    assert not any(r.startswith("GET /data/") for r in served.requests)
    done = _fetch(env, url, folder, root, CDR="0")
    assert done.returncode == 0, done.stdout[-3000:]
    data = [r for r in served.requests if r.startswith("GET /data/")]
    assert len(data) == len(set(data)) == 3 * len(STAR)
    posts = len([r for r in served.requests if r.startswith("POST")])
    again = _fetch(env, url, folder, root, CDR="0")
    assert again.returncode == 0, again.stdout[-3000:]
    assert [r for r in served.requests if r.startswith("GET /data/")] == data
    assert len([r for r in served.requests if r.startswith("POST")]) == posts


def test_a_partial_folder_fails_and_resumes(site, api):
    """A file lost from the folder makes the end check fail, naming it, and the same command
    fetches it again; a file that no longer matches its md5 is moved aside, never deleted.
    Catches an end check that passes a partial folder, and a mismatched file overwritten."""
    root, env = site
    _, url = api
    folder = root / "store" / "fake"
    assert _fetch(env, url, folder, root, KINDS="star", CDR="0").returncode == 0
    first, second = sorted(p for p in folder.glob("*/*.tsv"))[:2]
    first.unlink()
    second.write_text("changed")
    served, _ = api
    held = served.files.pop(first.parent.name)
    broken = _fetch(env, url, folder, root, KINDS="star", CDR="0", RETRIES="0")
    served.files[first.parent.name] = held
    assert broken.returncode == 1 and f"missing: {first}" in broken.stdout, broken.stdout[-2000:]
    files = len(STAR) - 1, 0, 1
    assert f"PASTE data files md5 OK, bad, missing: {', '.join(map(str, files))}   FAILS" in (
        broken.stdout)  # fmt: skip
    assert broken.stdout.rstrip().endswith("read the lines marked FAIL")
    fixed = _fetch(env, url, folder, root, KINDS="star", CDR="0")
    assert fixed.returncode == 0, fixed.stdout[-3000:]
    assert [p.read_text() for p in second.parent.glob("*.md5_mismatch_*")] == ["changed"]


def test_refuses_a_folder_of_another_project(site, api):
    """A folder whose saved file query names another project is refused before any request.
    Catches two projects' files mixed in one folder, which gather would read as one."""
    root, env = site
    _, url = api
    folder = root / "store" / "fake"
    assert _fetch(env, url, folder, root, KINDS="star", CDR="0").returncode == 0
    saved = folder / f"{gdc.METADATA}{gdc.QUERY}"
    saved.write_text(saved.read_text().replace("FAKE-BRCA", "FAKE-KIRC"))
    done = _fetch(env, url, folder, root, KINDS="star", CDR="0")
    assert done.returncode == 1 and "holds another query's files" in done.stdout


def test_kebnekaise_refuses_home_and_takes_the_store(site, api):
    """On Kebnekaise, the default site, a folder in home is refused before anything is fetched,
    and a folder on the store is fetched with the fetch environment's Python. Catches the home
    check dropped from this script, and the end check run on another Python."""
    root, env = site
    served, url = api
    home = Path(env["HOME"]) / "fake"
    run = {k: v for k, v in env.items() if k != "CTX_SITE"}
    script = Path(env["CTX_REPO"]) / "scripts" / "bash" / SCRIPT
    refused = _run([script, "FAKE-BRCA", home], root, run, GDC_API=url, KINDS="star", CDR="0")
    assert refused.returncode == 1 and "home folder" in refused.stderr
    assert served.requests == []
    store = Path(env["CTX_STORE"]) / "fake"
    done = _run([script, "FAKE-BRCA", store], root, run, GDC_API=url, KINDS="star", CDR="0")
    assert done.returncode == 0, done.stdout[-3000:] + done.stderr[-2000:]
    assert (root / "log" / "python3.calls").read_text().split() == ["called"]


def test_a_corrupt_transfer_is_never_kept(site, api):
    """A file whose bytes differ from the manifest's md5 on arrival stays a .part, is never
    moved into place, and the end check fails on it; served right, the same command fixes it.
    Catches the md5 check on arrival dropped, which would keep a damaged file as if it were
    GDC's."""
    root, env = site
    served, url = api
    folder = root / "store" / "fake"
    (star_file,) = [k for k, p in served.files.items() if p.parent.parent.name == "star"][:1]
    served.corrupt.add(star_file)
    broken = _fetch(env, url, folder, root, KINDS="star", CDR="0", RETRIES="0")
    assert broken.returncode == 1, broken.stdout[-3000:]
    assert f"FAIL md5 {folder}/{star_file}/" in broken.stdout
    assert not [p for p in (folder / star_file).iterdir() if not p.name.endswith(".part")]
    served.corrupt.clear()
    fixed = _fetch(env, url, folder, root, KINDS="star", CDR="0")
    assert fixed.returncode == 0, fixed.stdout[-3000:]


def test_a_record_cut_short_fails(site, api):
    """A clinical record whose pagination counts more hits than it holds fails the end check,
    named. Catches a record cut at the size limit passed as the whole project's."""
    root, env = site
    served, url = api
    served.short.add("/cases")
    done = _fetch(env, url, root / "store" / "fake", root, KINDS="star", CDR="0")
    assert done.returncode == 1, done.stdout[-3000:]
    cases = len({a.case for a in STAR})
    assert f"PASTE gdc_cases_clinical.json, warnings: {cases} of {cases + 1}; {{}}   FAILS" in (
        done.stdout)  # fmt: skip
