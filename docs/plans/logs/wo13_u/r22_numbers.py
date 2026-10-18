"""wo13-u R2.2: a proof root's gather checks and DE counts, read from the files it wrote.

Usage: python r22_numbers.py <root> <slug>
"""

import html
import re
import sys
from pathlib import Path

import polars as pl


def newest(folder: Path, pattern: str) -> Path:
    """The newest file of a pattern in a folder."""
    return sorted(folder.glob(pattern))[-1]


def main() -> None:
    """Print the gather page's checks, the left-out reasons and each contrast's calls."""
    root, slug = Path(sys.argv[1]).expanduser(), sys.argv[2]
    gather = root / "outputs" / "01_gather"
    page = newest(gather, f"*_gather_{slug}_overview_v1.html").read_text()
    page = re.sub(r"<script.*?</script>|<style.*?</style>", "", page, flags=re.S)
    shown = re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", page)))
    start = shown.find("Matrix checks")
    print(shown[start : shown.find("Interpretation", start)])
    sheet = pl.read_parquet(newest(gather, f"*_gather_{slug}_samplesheet_v1.parquet"))
    print("sheet:", sheet.height, "samples;", dict(sheet.group_by("condition").len().rows()))
    for results in sorted((root / "outputs" / "05_de").glob(f"*_de_{slug}_*_results_v1.parquet")):
        frame = pl.read_parquet(results)
        calls = dict(frame.group_by("call").len().rows())
        print(results.name, "tested", frame.height, "calls", calls)


if __name__ == "__main__":
    main()
