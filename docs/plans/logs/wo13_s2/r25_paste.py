"""wo13-s2 R2.5: paste each block a `ctx run --check` log prints into the dataset file.

As the guide's paired note asks; a block already in the file is refused, never pasted twice.

Usage: python r25_paste.py <check log> <dataset file>
"""

import sys
from pathlib import Path

MARK = "  add these lines to the dataset file:\n"


def main() -> None:
    """Append every printed block, each ending at its first blank line."""
    log, path = Path(sys.argv[1]).read_text(), Path(sys.argv[2])
    body = path.read_text()
    blocks = [part.split("\n\n", 1)[0].strip() + "\n" for part in log.split(MARK)[1:]]
    for block in blocks:
        if block in body:
            raise SystemExit(f"already in {path}: {block.splitlines()[0]}")
        body = body.rstrip("\n") + "\n\n" + block
    path.write_text(body)
    print(f"{path}: {len(blocks)} blocks pasted")


if __name__ == "__main__":
    main()
