"""wo13-u R2.2: scout GEO series whose authors' files are featureCounts or Ensembl-keyed counts.

Reads NCBI E-utilities summaries and each series' suppl folder listing; downloads nothing but
those listings. Usage: python r22_scout.py <term> <max>. No personal data is sent.
"""

import json
import re
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
WORDS = re.compile(r"featurecount|count|raw|matrix", re.I)


def get(url: str) -> str:
    """One GET, text."""
    with urllib.request.urlopen(url, timeout=60) as response:
        return response.read().decode("utf-8", "replace")


def listing(gse: str) -> list[tuple[str, str]]:
    """(file name, size as the FTP index prints it) of a series' suppl folder."""
    stub = gse[:-3] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/suppl/"
    try:
        page = get(url)
    except OSError:
        return []
    return re.findall(r'<a href="([^"?/]+)">[^<]*</a>\s+\S+\s+\S+\s+(\S+)', page)


def main() -> None:
    """Print each series' accession, sample count, title and its count-like files."""
    term, top = sys.argv[1], int(sys.argv[2])
    query = urllib.parse.urlencode({"db": "gds", "retmax": top, "retmode": "json", "term": term})
    ids = json.loads(get(f"{EUTILS}/esearch.fcgi?{query}"))["esearchresult"]["idlist"]
    summary = json.loads(get(f"{EUTILS}/esummary.fcgi?db=gds&retmode=json&id={','.join(ids)}"))
    for uid in ids:
        record = summary["result"][uid]
        gse = "GSE" + record["gse"]
        files = [(n, s) for n, s in listing(gse) if WORDS.search(n)]
        print(f"{gse}\tn={record['n_samples']}\t{record['title'][:90]}")
        for name, size in files:
            print(f"    {name}\t{size}")
        time.sleep(0.4)


if __name__ == "__main__":
    main()
