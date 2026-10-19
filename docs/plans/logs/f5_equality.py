"""F.5's equality proof: every dataset file and many broken copies, loaded and dumped.

Run once before the split and once after, from the worktree root, and compare the two dumps
byte for byte: ``python docs/plans/logs/f5_equality.py <scratch dir> <dump file>``.
Each valid file is loaded as it is; then every ``key = value`` line of it is deleted, set to
an integer, a boolean and a string in turn, and finally all of them together, so one file
holds many problems and their order is compared too.
"""

import re
import sys
from pathlib import Path

from ctx.dataset import load, loaded

ROOT = Path.cwd()
SOURCES = sorted(
    [*ROOT.glob("conf/datasets/*.toml"), *ROOT.glob("tests/fixtures/datasets/*.toml"),
     *ROOT.glob("conf/dataset_template*.toml")]
)  # fmt: skip
KEY = re.compile(r"^(\s*[A-Za-z_][\w.\"-]*\s*=\s*)(.+)$")
SWAPS = ("", "7", "true", '"x"')


def outcome(path: Path) -> str:
    """The parsed Dataset's repr, or the exception's type and message."""
    try:
        return "OK " + repr(load(path))
    except Exception as error:  # noqa: BLE001 - every outcome is recorded, whatever it is
        return f"{type(error).__name__} {error}"


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


def main(scratch: Path, dump: Path) -> None:
    """Write every outcome, one per line, to ``dump``."""
    scratch.mkdir(parents=True, exist_ok=True)
    rows, broken = [], 0
    for source in SOURCES:
        rows.append(f"{source.relative_to(ROOT)}\t{outcome(source)}")
        for j, text in enumerate(variants(source.read_text())):
            # Beside the source, so [sdrf] file and subset_of resolve as they do there.
            path = source.with_name(f"zz_f5_{source.stem}_{j}.toml")
            path.write_text(text)
            try:
                result = outcome(path).replace(str(path), "<variant>")
            finally:
                path.unlink()
            rows.append(f"{source.relative_to(ROOT)}#{j}\t{result}")
            broken += not result.startswith("OK ")
    rows.append(f"loaded\t{sorted(loaded())}")
    dump.write_text("\n".join(rows) + "\n")
    print(f"{len(SOURCES)} files, {len(rows) - 1 - len(SOURCES)} variants, {broken} not loaded")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
