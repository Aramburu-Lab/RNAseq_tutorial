"""wo13-y item 5: every real dataset file and template, loaded, one repr line per file.

Usage: python f6a_real_files.py <out>. Run once with PYTHONPATH at an export of main's src and
once at the worktree's; the two outputs must be equal. Prints which ctx it imported, so the
comparison cannot be the worktree against itself.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "tests"))
from dataset_variants import ROOT, sources  # noqa: E402

import ctx  # noqa: E402
from ctx.dataset import DatasetError, load  # noqa: E402

print("ctx from", ctx.__file__)
lines = []
for path in sources():
    try:
        got = "OK " + repr(load(path))
    except DatasetError as error:
        got = f"DatasetError {error}"
    lines.append(f"{path.relative_to(ROOT)}\t{got}")
Path(sys.argv[1]).write_text("\n".join(lines) + "\n")
print(len(lines), "files,", sum(x.split("\t")[1].startswith("OK") for x in lines), "load")
