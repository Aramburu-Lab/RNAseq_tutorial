"""F.5's broken dataset files, one definition: each key line deleted or retyped (F.6).

tests/test_dataset_variants.py loads a sample of them; docs/plans/logs/g6_classify.py loads
every one. The rule is F.5's (docs/plans/logs/f5_equality.py), which ran it before and after
the split; that log keeps its own copy, as a record of what it ran.
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
KEY = re.compile(r"^(\s*[A-Za-z_][\w.\"-]*\s*=\s*)(.+)$")
SWAPS = ("", "7", "true", '"x"')


def sources(root: Path = ROOT) -> list[Path]:
    """Every dataset file and template the variants are made from, in a fixed order."""
    return sorted(
        [*root.glob("conf/datasets/*.toml"), *root.glob("tests/fixtures/datasets/*.toml"),
         *root.glob("conf/dataset_template*.toml")]
    )  # fmt: skip


def variants(text: str) -> list[str]:
    """Each key line deleted or retyped, one at a time, then every line retyped together."""
    lines = text.splitlines()
    out = []
    for i, line in enumerate(lines):
        if not (match := KEY.match(line)):
            continue
        for swap in SWAPS:
            new = match.group(1) + swap if swap else ""
            out.append("\n".join([*lines[:i], new, *lines[i + 1 :]]))
    for swap in SWAPS[1:]:
        out.append("\n".join(
            (m.group(1) + swap) if (m := KEY.match(line)) else line for line in lines
        ))  # fmt: skip
    out.append(text + "\n[bogus]\nx = 1\n[clinical]\ngenes = []\n")
    return out
