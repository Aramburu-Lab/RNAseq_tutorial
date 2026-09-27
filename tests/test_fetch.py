"""ctx setup's resource fetch, served by http.server on 127.0.0.1; one network test, marked."""

import functools
import hashlib
import json
import os
import subprocess
import sys
import threading
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pytest

from ctx.dataset import Resource, load
from ctx.setup import MANIFEST, FetchError, fetch

ROOT = Path(__file__).resolve().parents[1]
CTX = Path(sys.executable).parent / "ctx"
BODY = b"gene set\tdescription\tA\tB\tC\n"


@pytest.fixture
def server(tmp_path):
    served = tmp_path / "served"
    served.mkdir()
    (served / "good.txt").write_bytes(BODY)
    (served / "tampered.txt").write_bytes(BODY.replace(b"A", b"Z"))
    hits: list[str] = []

    class Counting(SimpleHTTPRequestHandler):
        def do_GET(self):  # noqa: N802 - the http.server hook name
            hits.append(self.path)
            super().do_GET()

        def log_message(self, *args):
            pass

    httpd = ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(Counting, directory=served))
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{httpd.server_port}", hits
    httpd.shutdown()


def resource(base: str, *names: str, dest: str = "resources/sets/good.gmt") -> Resource:
    digest = hashlib.sha256(BODY).hexdigest()
    return Resource("sets", tuple(f"{base}/{n}" for n in names), len(BODY), digest, dest)


def test_fetch_verifies(server, tmp_path):
    base, _ = server
    (entry,) = fetch(tmp_path / "root", [resource(base, "good.txt")])
    dest = tmp_path / "root" / "resources" / "sets" / "good.gmt"
    assert dest.read_bytes() == BODY
    manifest = json.loads((tmp_path / "root" / MANIFEST).read_text())
    assert manifest["resources/sets/good.gmt"]["sha256"] == hashlib.sha256(BODY).hexdigest()
    assert entry["url"] == f"{base}/good.txt"


def test_tampered_refused(server, tmp_path):
    base, _ = server
    with pytest.raises(FetchError, match="tampered.txt"):
        fetch(tmp_path / "root", [resource(base, "tampered.txt")])
    folder = tmp_path / "root" / "resources" / "sets"
    assert not list(folder.iterdir()), "neither the file nor a .part may remain"


def test_mirror_fallback(server, tmp_path):
    base, hits = server
    (entry,) = fetch(tmp_path / "root", [resource(base, "missing.txt", "good.txt")])
    assert entry["url"] == f"{base}/good.txt"
    assert hits == ["/missing.txt", "/good.txt"]


def test_idempotent(server, tmp_path):
    base, hits = server
    fetch(tmp_path / "root", [resource(base, "good.txt")])
    before = len(hits)
    fetch(tmp_path / "root", [resource(base, "good.txt")])
    assert len(hits) == before


def test_dest_outside_root_refused(server, tmp_path):
    base, hits = server
    with pytest.raises(FetchError, match="outside"):
        fetch(tmp_path / "root", [resource(base, "good.txt", dest="../escape.gmt")])
    assert hits == []


def test_cli_setup_fetches(server, tmp_path):
    base, _ = server
    digest = hashlib.sha256(BODY).hexdigest()
    spec = tmp_path / "local.toml"
    spec.write_text(
        f'[dataset]\nslug = "local"\n\n[[resources]]\nname = "sets"\nurls = ["{base}/good.txt"]\n'
        f'bytes = {len(BODY)}\nsha256 = "{digest}"\ndest = "resources/sets/good.gmt"\n\n'
        "[guard]\nids = []\npatterns = []\n"
    )
    root = tmp_path / "root"
    root.mkdir()
    env = {k: v for k, v in os.environ.items() if not k.startswith(("SLURM_", "CTX_"))}
    done = subprocess.run(
        [str(CTX), "--root", str(root), "setup", "--dataset-file", str(spec)],
        cwd=root, env=env, capture_output=True, text=True, check=False,
    )  # fmt: skip
    assert done.returncode == 0, done.stderr
    assert (root / "resources" / "sets" / "good.gmt").read_bytes() == BODY
    (page,) = (root / "outputs" / "00_setup").glob("*.html")
    assert f"sha256 {digest}" in page.read_text()


@pytest.mark.network
def test_fetch_real(tmp_path):
    spec = load(ROOT / "conf" / "datasets" / "tier1.toml")
    entries = fetch(tmp_path, spec.resources)
    assert [e["name"] for e in entries] == [r.name for r in spec.resources]
    for r in spec.resources:
        assert hashlib.sha256((tmp_path / r.dest).read_bytes()).hexdigest() == r.sha256
