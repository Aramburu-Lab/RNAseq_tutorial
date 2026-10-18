"""wo13-r1: open STAR counts, miRNA and isoform files per candidate project, from the GDC API."""

import json
import sys
import urllib.request

API = "https://api.gdc.cancer.gov"
KINDS = {
    "star": ("Gene Expression Quantification", "STAR - Counts"),
    "mirna": ("miRNA Expression Quantification", "BCGSC miRNA Profiling"),
    "isoform": ("Isoform Expression Quantification", "BCGSC miRNA Profiling"),
}


def query(project: str, data_type: str, workflow: str) -> dict:
    """Files of one kind: count, bytes and sample types, as GDC lists them."""
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [project]}},
            {"op": "in", "content": {"field": "data_type", "value": [data_type]}},
            {"op": "in", "content": {"field": "analysis.workflow_type", "value": [workflow]}},
            {"op": "in", "content": {"field": "access", "value": ["open"]}},
        ],
    }
    body = {"filters": filters, "size": 10000, "fields": "file_size,cases.samples.sample_type"}
    req = urllib.request.Request(
        f"{API}/files", json.dumps(body).encode(), {"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=120) as r:  # noqa: S310 - the GDC API
        hits = json.load(r)["data"]["hits"]
    types: dict[str, int] = {}
    for h in hits:
        t = h["cases"][0]["samples"][0]["sample_type"]
        types[t] = types.get(t, 0) + 1
    return {"files": len(hits), "bytes": sum(h["file_size"] for h in hits), "types": types}


def main() -> None:
    """Print each project's three kinds and their total bytes."""
    with urllib.request.urlopen(f"{API}/status", timeout=60) as r:  # noqa: S310 - the GDC API
        print("status", json.load(r)["data_release"])
    for project in sys.argv[1:]:
        total = 0
        for kind, (dt, wf) in KINDS.items():
            got = query(project, dt, wf)
            total += got["bytes"]
            print(project, kind, got["files"], got["bytes"], json.dumps(got["types"]))
        print(project, "all three kinds", total)


if __name__ == "__main__":
    main()
