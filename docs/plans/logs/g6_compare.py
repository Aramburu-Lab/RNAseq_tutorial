"""F.6's comparison of two g6_classify.py dumps, before and after the fix, variant by variant.

``python docs/plans/logs/g6_compare.py <before dump> <after dump>`` prints each class of
change and every variant that changed in a way the brief does not allow: a parsed Dataset
whose repr differs or that is now refused, and a DatasetError message that lost a problem.
A message that keeps every old problem, in order, and adds some is listed with what it adds.
"""

import re
import sys
from collections import Counter
from pathlib import Path

ROW = re.compile(r"^(\S+#\d+)\t", re.MULTILINE)


def rows(dump: Path) -> dict[str, str]:
    """Variant id to its outcome text, which can span lines."""
    parts = ROW.split(dump.read_text())
    return dict(zip(parts[1::2], (p.rstrip("\n") for p in parts[2::2]), strict=True))


def kept_in_order(old: list[str], new: list[str]) -> bool:
    """Every old line is in new, in the same order."""
    found = iter(new)
    return all(line in found for line in old)


def classify(old: str, new: str) -> str:
    """One variant's change, as a class name."""
    kind_old, kind_new = old.split(" ", 1)[0], new.split(" ", 1)[0]
    if old == new:
        return f"{kind_old} equal"
    if kind_old not in ("OK", "DatasetError"):
        return f"{kind_old} to {kind_new}"
    if kind_old == kind_new == "DatasetError" and kept_in_order(old.splitlines(),
                                                               new.splitlines()):  # fmt: skip
        return "DatasetError adds problems"
    return f"{kind_old} to {kind_new}, CHANGED"


def main(before: Path, after: Path) -> None:
    """Print the classes, then each added problem and each disallowed change."""
    old, new = rows(before), rows(after)
    assert old.keys() == new.keys(), "the two dumps hold different variants"
    classes = {key: classify(old[key], new[key]) for key in old}
    for name, n in sorted(Counter(classes.values()).items()):
        print(f"{n}\t{name}")
    for key, name in classes.items():
        if name == "DatasetError adds problems":
            added = [x for x in new[key].splitlines() if x not in old[key].splitlines()]
            print(f"ADDS\t{key}\t{' | '.join(a.strip() for a in added)}")
        elif "CHANGED" in name:
            print(f"CHANGED\t{key}\n  before: {old[key][:300]}\n  after:  {new[key][:300]}")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
