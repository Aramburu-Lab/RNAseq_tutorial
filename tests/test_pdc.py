"""A PDC resource through ``ctx setup``'s fetch: a fresh signed URL, PDC's MD5, then SHA-256.

The opener is a stand-in for PDC: it answers the GraphQL POST and serves the signed URL. No
network. Each test names the mutation it catches.
"""

import hashlib
import io
import json

import pytest

from ctx.dataset import Resource
from ctx.proteo import pdc
from ctx.setup import FetchError, fetch

BODY = b"a made-up freeze archive"
URL = "pdc://PDCFAKE/freeze.tar.gz"


def opener(md5: str = hashlib.md5(BODY).hexdigest(), listed: int = 1):  # noqa: S324 - PDC's checksum
    """A stand-in for urlopen: PDC's GraphQL answer, then the file at its signed URL."""
    seen = []

    def open_(target, timeout):
        if hasattr(target, "data"):
            query = json.loads(target.data)["query"]
            seen.append(query)
            files = [{"file_name": "freeze.tar.gz", "md5sum": md5, "file_size": len(BODY),
                      "signedUrl": {"url": "https://signed.example/once"}}] * listed  # fmt: skip
            return io.BytesIO(json.dumps({"data": {"filesPerStudy": files}}).encode())
        seen.append(target)
        assert target == "https://signed.example/once"
        return io.BytesIO(BODY)

    return open_, seen


def resource() -> Resource:
    return Resource("freeze", (URL,), len(BODY), hashlib.sha256(BODY).hexdigest(),
                    "resources/freeze.tar.gz")  # fmt: skip


def test_resolved_fetched_and_recorded_by_its_pdc_url(tmp_path, capsys):
    """Mutations caught: the pdc:// URL handed to urlopen unresolved; the signed URL recorded in
    the manifest, where it expires; the study ID left out of the query; the agreement accepted
    without saying so first (Rytis, 2026-10-08)."""
    open_, seen = opener()
    (entry,) = fetch(tmp_path, [resource()], opener=open_)
    said = capsys.readouterr().out
    assert "accepts the CPTAC data use agreement" in said and pdc.DUA_URL in said
    assert (tmp_path / "resources" / "freeze.tar.gz").read_bytes() == BODY
    assert entry["url"] == URL and 'pdc_study_id: "PDCFAKE"' in seen[0]
    assert seen[1] == "https://signed.example/once"


def test_md5_checked_against_pdc(tmp_path):
    """Mutation caught: PDC's MD5 never compared, so only the declared SHA-256 guards the file."""
    open_, _ = opener(md5="0" * 32)
    with pytest.raises(FetchError, match="PDC lists"):
        fetch(tmp_path, [resource()], opener=open_)
    assert not (tmp_path / "resources" / "freeze.tar.gz").exists()


REFUSED = [("pdc://PDCFAKE", 1, "expected pdc://"), (URL, 2, "2 times"), (URL, 0, "0 times")]


@pytest.mark.parametrize(("url", "listed", "said"), REFUSED)
def test_resolve_refuses(url, listed, said):
    """Mutations caught: a URL without a file name sent to PDC; a file listed twice taken first."""
    open_, _ = opener(listed=listed)
    with pytest.raises(OSError, match=said):
        pdc.resolve(url, open_)
