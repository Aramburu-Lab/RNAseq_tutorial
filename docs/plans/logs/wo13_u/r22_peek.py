"""wo13-u R2.2: the first lines of a gzipped GEO file from its first 32 KiB, nothing more.

Usage: python r22_peek.py <url> [lines]
"""

import sys
import urllib.request
import zlib


def main() -> None:
    """Print the first lines a partial read decompresses to, each cut at 300 characters."""
    url, *rest = sys.argv[1:]
    lines = int(rest[0]) if rest else 3
    request = urllib.request.Request(url, headers={"Range": "bytes=0-32767"})
    with urllib.request.urlopen(request, timeout=60) as response:
        data = response.read()
    text = zlib.decompressobj(31).decompress(data).decode("utf-8", "replace")
    for line in text.splitlines()[:lines]:
        print(line[:300])


if __name__ == "__main__":
    main()
