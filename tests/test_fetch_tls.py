"""ctx setup's fetch from a server that, like miRTarBase's, offers only RSA key exchange.

A local TLS 1.2 server on 127.0.0.1 is limited to AES128-GCM-SHA256, the suite miRTarBase's
server chose on 2026-10-08, with a certificate made for the test and trusted through
SSL_CERT_FILE, so certificate checking stays on throughout.
Mutations these tests catch: ``_download`` dropping the resource's context (the declared
fetch fails the handshake); ``tls_context`` replacing the default list instead of appending
(the forward-secret suites vanish); ``tls_context`` setting CERT_NONE (the untrusted server
is accepted); the loader ignoring ``tls_ciphers`` (the brca entry declares nothing).
"""

import functools
import hashlib
import shutil
import ssl
import subprocess
import sys
import threading
from dataclasses import replace
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pytest

from ctx.dataset import Resource, load
from ctx.setup import FetchError, fetch, tls_context

ROOT = Path(__file__).resolve().parents[1]
BRCA = ROOT / "conf" / "datasets" / "tcga_brca.toml"
BODY = b"miRTarBase ID,miRNA,Target Gene\nMIRT000002,hsa-miR-20a-5p,HIF1A\n"
ONLY_SUITE = "AES128-GCM-SHA256"


def _certificate(folder: Path, name: str) -> tuple[Path, Path]:
    """A one-day self-signed certificate for 127.0.0.1, made by the env's openssl."""
    openssl = shutil.which("openssl", path=str(Path(sys.executable).parent)) or "openssl"
    key, cert = folder / f"{name}.key", folder / f"{name}.pem"
    subprocess.run(
        [openssl, "req", "-x509", "-newkey", "rsa:2048", "-nodes", "-days", "1",
         "-keyout", str(key), "-out", str(cert), "-subj", "/CN=127.0.0.1",
         "-addext", "subjectAltName=IP:127.0.0.1"],
        check=True, capture_output=True,
    )  # fmt: skip
    return cert, key


@pytest.fixture
def rsa_only(tmp_path, monkeypatch):
    """An HTTPS server offering one RSA key-exchange suite; yields its base URL and a stranger."""
    served = tmp_path / "served"
    served.mkdir()
    (served / "targets.csv").write_bytes(BODY)
    cert, key = _certificate(tmp_path, "server")
    stranger, _ = _certificate(tmp_path, "stranger")
    monkeypatch.setenv("SSL_CERT_FILE", str(cert))
    server_context = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)
    server_context.maximum_version = ssl.TLSVersion.TLSv1_2
    server_context.set_ciphers(ONLY_SUITE)
    server_context.load_cert_chain(cert, key)

    class Quiet(SimpleHTTPRequestHandler):
        def log_message(self, *args):
            pass

    httpd = ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=served))
    httpd.socket = server_context.wrap_socket(httpd.socket, server_side=True)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield f"https://127.0.0.1:{httpd.server_port}", stranger
    httpd.shutdown()


def _declared() -> Resource:
    """The miRTarBase entry exactly as tcga_brca.toml declares it."""
    return load(BRCA).resource("mirtarbase_strong")


def _local(base: str, ciphers: str) -> Resource:
    """The brca entry pointed at the local server, with its body's size and hash."""
    return replace(
        _declared(), urls=(f"{base}/targets.csv",), bytes=len(BODY),
        sha256=hashlib.sha256(BODY).hexdigest(), dest="data/raw/mirna_targets/t.csv",
        tls_ciphers=ciphers,
    )  # fmt: skip


def test_brca_declares_the_rsa_suites():
    assert ONLY_SUITE in _declared().tls_ciphers.split(":")


def test_declared_ciphers_fetch_and_default_fails(rsa_only, tmp_path):
    base, _ = rsa_only
    with pytest.raises(FetchError, match="HANDSHAKE_FAILURE"):
        fetch(tmp_path / "plain", [_local(base, "")])
    (entry,) = fetch(tmp_path / "root", [_local(base, _declared().tls_ciphers)])
    assert (tmp_path / "root" / "data/raw/mirna_targets/t.csv").read_bytes() == BODY
    assert entry["sha256"] == hashlib.sha256(BODY).hexdigest()


def test_certificate_still_checked(rsa_only, tmp_path, monkeypatch):
    base, stranger = rsa_only
    monkeypatch.setenv("SSL_CERT_FILE", str(stranger))
    with pytest.raises(FetchError, match="CERTIFICATE_VERIFY_FAILED"):
        fetch(tmp_path / "root", [_local(base, _declared().tls_ciphers)])


def test_context_appends_to_the_default_list():
    default = {c["name"] for c in ssl.create_default_context().get_ciphers()}
    context = tls_context(_declared().tls_ciphers)
    names = {c["name"] for c in context.get_ciphers()}
    assert default < names
    assert ONLY_SUITE in names - default
    assert context.verify_mode == ssl.CERT_REQUIRED and context.check_hostname


@pytest.mark.network
def test_mirtarbase_real(tmp_path):
    (entry,) = fetch(tmp_path, [_declared()])
    assert entry["bytes"] == _declared().bytes
