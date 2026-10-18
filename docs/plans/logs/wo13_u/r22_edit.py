"""wo13-u R2.2: the EDIT lines of a copied counts template, as the guide's edit note asks.

Each replacement must match once; the file is the root's own copy, never the template.
Usage: python r22_edit.py <dataset file> <edits.json>, edits a list of [old, new].
"""

import json
import sys
from pathlib import Path


def main() -> None:
    """Apply every replacement or stop at the first that does not match exactly once."""
    path, edits = Path(sys.argv[1]), json.loads(Path(sys.argv[2]).read_text())
    body = path.read_text()
    for old, new in edits:
        if body.count(old) != 1:
            raise SystemExit(f"{old!r} matches {body.count(old)} times in {path}")
        body = body.replace(old, new)
    path.write_text(body)
    print(f"{path}: {len(edits)} edits")


if __name__ == "__main__":
    main()
