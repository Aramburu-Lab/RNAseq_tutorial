"""wo13-u R2.2: copy a proof root's dataset file and table into conf/datasets, headed.

The root's files stay as the proof wrote them; the copies carry a header naming the proof.
Usage: python r22_keep.py <root> <slug> <conf/datasets folder> <header line>...
"""

import shutil
import sys
from pathlib import Path


def main() -> None:
    """Write <slug>.toml with the header lines above the root's text, and copy <slug>.csv."""
    root, slug, out, *header = sys.argv[1:]
    source, folder = Path(root).expanduser(), Path(out)
    lines = "".join(f"# {line}\n" for line in header)
    (folder / f"{slug}.toml").write_text(lines + "\n" + (source / f"{slug}.toml").read_text())
    shutil.copyfile(source / f"{slug}.csv", folder / f"{slug}.csv")
    print(folder / f"{slug}.toml", folder / f"{slug}.csv")


if __name__ == "__main__":
    main()
