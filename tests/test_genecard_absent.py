"""M9's cache records "no such entry" (HTTP 404) as an answer; a 5xx, a timeout or a lost
network is never cached (NC, 2026-10-20).

Each test names the mutation it catches; every one was planted and seen to fail. No test
sends a request beyond the loopback server it starts itself.
"""

import http.server
import json
import re
import threading
import time
import urllib.error
from datetime import date
from pathlib import Path

from genecard_root import ACTB, SLUG, Stub, make, responses

from ctx import runtime
from ctx.genecard import annotation, hpa, junctions
from ctx.genecard import cache as cache_module
from ctx.genecard import run as module
from ctx.genecard.cache import MANIFEST, SUPERSEDED, Cache, absent
from ctx.report import genes
from ctx.report.builder import unpacked
from ctx.setup import fetch

GENE = "ENSG00000106236"
URL = hpa.url("v25", GENE, "xml")
NOT_CACHED = "not cached, and --offline makes no request"


class Failing(Stub):
    """A stub answering chosen URLs with a 500 or a timeout, and the rest as ``Stub`` does."""

    def __init__(self, bodies: dict[str, bytes], errors: dict[str, int | str]) -> None:
        """Serve ``bodies``; answer each URL of ``errors`` with its status or "timeout"."""
        super().__init__(bodies)
        self.errors = errors

    def __call__(self, url: str):  # noqa: ANN204 - Stub's own return
        """The body of ``url``, its planted error, or an HTTP 404."""
        if (error := self.errors.get(url)) == "timeout":
            self.calls.append(url)
            raise urllib.error.URLError(TimeoutError("timed out"))
        if isinstance(error, int):
            self.calls.append(url)
            raise urllib.error.HTTPError(url, error, "Server Error", None, None)  # type: ignore[arg-type]
        return super().__call__(url)


def _ask(root: Path, opener, **kwargs) -> dict[str, cache_module.Response]:  # noqa: ANN001
    cache = Cache(root, opener=opener, **kwargs)
    return {name: cache.get("hpa", "v25", GENE, name, f"https://x/{name}")
            for name in ("missing", "broken", "slow")}  # fmt: skip


def test_404_is_recorded_and_5xx_and_timeouts_are_not(tmp_path):
    # Mutations: a 5xx or a timeout cached; a 404 not cached; the offline wording of a 404
    # differing from the online one; a recorded 404 asked again before --refresh.
    stub = Failing({}, {"https://x/broken": 500, "https://x/slow": "timeout"})
    online = _ask(tmp_path, stub)
    assert all(a.body is None for a in online.values())
    assert online["missing"].error == absent(404, online["missing"].retrieved)
    assert online["broken"].error == "HTTP Error 500: Server Error"
    assert "timed out" in online["slow"].error
    recorded = json.loads((tmp_path / MANIFEST).read_text())
    assert list(recorded) == [f"hpa/v25/{GENE}/missing.gz"]
    entry = recorded[f"hpa/v25/{GENE}/missing.gz"]
    assert (entry["url"], entry["status"], entry["release"]) == ("https://x/missing", 404, "v25")
    assert entry["retrieved"] == online["missing"].retrieved
    offline = _ask(tmp_path, stub, offline=True)
    assert offline["missing"].error == online["missing"].error
    assert offline["broken"].error == offline["slow"].error == NOT_CACHED
    assert len(stub.calls) == 3  # noqa: PLR2004 - offline asks nothing
    again = _ask(tmp_path, stub)
    assert again["missing"].error == online["missing"].error
    assert stub.calls[3:] == ["https://x/broken", "https://x/slow"]


def test_refresh_moves_a_recorded_absence_aside_and_asks_again(tmp_path):
    # Mutation: a recorded 404 kept for ever, or refreshed by deleting it.
    first = Cache(tmp_path, opener=Stub({})).get("hpa", "v25", GENE, "xml", URL)
    assert first.body is None and first.error.startswith("no entry: HTTP 404")
    stub = Stub({URL: b"<entry/>"})
    fresh = Cache(tmp_path, refresh=True, opener=stub, today=date(2026, 10, 20))
    assert fresh.get("hpa", "v25", GENE, "xml", URL).body == b"<entry/>"
    assert stub.calls == [URL]
    assert (tmp_path / SUPERSEDED / "2026-10-20" / "hpa" / "v25" / GENE / "xml.gz").is_file()
    assert json.loads((tmp_path / MANIFEST).read_text())[f"hpa/v25/{GENE}/xml.gz"]["status"] == 200


class _Handler(http.server.BaseHTTPRequestHandler):
    """Loopback answers: 404 with a body, 500, and a reply slower than the timeout."""

    def do_GET(self) -> None:  # noqa: N802 - http.server's name
        """Answer by path."""
        if self.path == "/slow":
            time.sleep(2)
        status = {"/missing": 404, "/broken": 503}.get(self.path, 200)
        body = b"no such gene" if status == 404 else b"down"  # noqa: PLR2004
        self.send_response(status)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, *args: object) -> None:
        """Keep the test log quiet."""


def test_real_urllib_against_a_loopback_server(tmp_path, monkeypatch):
    # Mutation: the 404 caught as a URLError before its status is read (the defect found by
    # the kit v10 build), which only the real urlopen raises the way HPA's server makes it.
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), _Handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    monkeypatch.setattr(cache_module, "TIMEOUT", 0.5)
    base = f"http://127.0.0.1:{server.server_address[1]}"
    try:
        cache = Cache(tmp_path)
        answers = {p: cache.get("t", "v1", GENE, p, f"{base}/{p}") for p in ("missing", "broken",
                                                                            "slow")}  # fmt: skip
    finally:
        server.shutdown()
    assert answers["missing"].error.startswith("no entry: HTTP 404 on ")
    assert answers["broken"].error == "HTTP Error 503: Service Unavailable"
    assert "timed out" in answers["slow"].error
    entry = json.loads((tmp_path / MANIFEST).read_text())
    assert list(entry) == [f"t/v1/{GENE}/missing.gz"]
    assert entry[f"t/v1/{GENE}/missing.gz"]["bytes"] == len(b"no such gene")


def _unavailable(page: Path) -> dict[str, str]:
    text = re.sub(r"<[^>]+>", " ", unpacked(page.read_text()))
    found = re.findall(r"(Human Protein Atlas|Ensembl REST|GTEx V10 junctions) unavailable for "
                       r"ACTB: (.*?)\.\s", text)  # fmt: skip
    return dict(found)


# What a run states about itself, set aside when an offline rerun is compared with the online
# page (K.16): its measured cost, and the --offline its follow-up commands carry so they keep
# the run's network rule.
_COST = re.compile(r'(<td class="(?:wall|rss-self|rss-children)">)[^<]*(</td>)')


def _evidence(html: str) -> str:
    return _COST.sub(r"\1\2", html).replace(" --offline", "")


def test_offline_rerun_prints_the_online_pages_no_entry(tmp_path):
    # Mutations: the offline page saying "not cached" for HPA's 404; a 500 or a timeout
    # shown offline as if cached; any other word of the page differing offline, such as a
    # "cached" beside a retrieval date (K.16); a follow-up command dropping its --offline.
    folder, spec = tmp_path / "root", make(tmp_path / "root")
    from ctx.dataset import load  # noqa: PLC0415

    fetch(folder, load(spec).resources, opener=lambda url, **_: None)
    actb = genes.ensembl(ACTB)
    stub = Failing(responses(), {annotation.lookup_url(actb): 500,
                                 junctions.gene_url(actb): "timeout"})  # fmt: skip
    settings = runtime.Settings("laptop", folder, 1, None, 0.0)
    from test_genecard import _cards  # noqa: PLC0415

    online = module.run(settings, _cards(), "ctx genecard", SLUG,
                        module.Options(spec, genes=("NPTX2", "ACTB")), stub)  # fmt: skip
    said, page = _unavailable(online[0]), unpacked(online[0].read_text())
    offline = module.run(settings, _cards(), "ctx genecard", SLUG,
                         module.Options(spec, genes=("NPTX2", "ACTB"), offline=True))  # fmt: skip
    again, rerun = _unavailable(offline[0]), unpacked(offline[0].read_text())
    no_entry = re.fullmatch(r"xml: no entry: HTTP 404 on \d{4}-\d{2}-\d{2}, cached until "
                            r"--refresh; json: no entry: HTTP 404 on \d{4}-\d{2}-\d{2}, cached "
                            r"until --refresh", said["Human Protein Atlas"])  # fmt: skip
    assert no_entry, said
    assert again["Human Protein Atlas"] == said["Human Protein Atlas"]
    assert said["Ensembl REST"] == "HTTP Error 500: Server Error"
    assert "timed out" in said["GTEx V10 junctions"]
    assert again["Ensembl REST"] == again["GTEx V10 junctions"] == NOT_CACHED
    # The rest of the rerun is the online page, but for the two failures no cache records.
    expected = page.replace(said["Ensembl REST"], NOT_CACHED)
    expected = expected.replace(said["GTEx V10 junctions"], NOT_CACHED)
    assert _evidence(rerun) == _evidence(expected)
    assert rerun.count(" --offline") > expected.count(" --offline")
