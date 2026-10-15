"""A Proteomic Data Commons file as a ``ctx setup`` resource: ``pdc://<study>/<file name>``.

PDC serves a file through a signed URL that expires, so a dataset file cannot hold the URL
itself. ``resolve`` asks PDC's GraphQL API for a fresh one, with PDC's MD5 and size, when setup
fetches the file; setup then checks the MD5 against PDC's and the SHA-256 against the dataset
file's (M7 plan, P7.5). The query accepts the CPTAC data use agreement, as the PDC portal's
download does, for whoever runs ``ctx setup`` (Rytis, 2026-10-08); setup says so first.
"""

import json
import urllib.request
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

SCHEME = "pdc://"
GRAPHQL = "https://pdc.cancer.gov/graphql"
# The PDC data categories a resource may sit in: a paper's tables, and the study's own metadata
# workbooks (TMT plex maps, clinical data), which PDC files as Other Metadata.
CATEGORIES = ("Publication Supplementary Material", "Other Metadata")
QUERY = ('{ filesPerStudy(pdc_study_id: "%s", data_category: "%s", acceptDUA: true) '
         "{ file_name md5sum file_size signedUrl { url } } }")  # fmt: skip
# The agreement the query accepts; ctx setup names it before fetching (Rytis, 2026-10-08).
DUA_URL = "https://cptac-data-portal.georgetown.edu/cptac/s?id=73217"
TIMEOUT = 120


@dataclass(frozen=True)
class Signed:
    """A file's fresh download URL with PDC's own MD5 and size."""

    url: str
    md5: str
    size: int


def resolve(url: str, opener: Callable[..., Any] = urllib.request.urlopen) -> Signed:
    """The signed URL, MD5 and size PDC gives now for ``pdc://<study>/<file name>``.

    Every category of ``CATEGORIES`` is asked; the file must be listed once over all of them.

    Raises:
        OSError: If the URL is malformed, PDC answers with errors, or lists the file other
            than once; setup records it as that mirror's failure.
    """
    study, _, name = url.removeprefix(SCHEME).partition("/")
    if not url.startswith(SCHEME) or not study or not name:
        raise OSError(f"{url}: expected {SCHEME}<study>/<file name>")
    files = []
    for category in CATEGORIES:
        body = json.dumps({"query": QUERY % (study, category)}).encode()
        request = urllib.request.Request(  # noqa: S310 - a fixed https endpoint
            GRAPHQL, data=body, headers={"Content-Type": "application/json"}, method="POST"
        )
        with opener(request, timeout=TIMEOUT) as response:
            answer = json.loads(response.read())
        if answer.get("errors"):
            raise OSError(f"{url}: PDC answered {answer['errors']}")
        files += [f for f in (answer.get("data") or {}).get("filesPerStudy") or []
                  if f.get("file_name") == name]  # fmt: skip
    if len(files) != 1 or not (files[0].get("signedUrl") or {}).get("url"):
        raise OSError(f"{url}: PDC lists the file {len(files)} times or without a signed URL")
    found = files[0]
    return Signed(found["signedUrl"]["url"], str(found["md5sum"]), int(found["file_size"]))
