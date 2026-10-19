"""F.6: every F.5 variant (tests/dataset_variants.py) loaded, each non-DatasetError by cause.

Run from the worktree root: ``python docs/plans/logs/g6_classify.py <dump file> <table file>``.
The dump has one line per variant, as f5_equality.py writes them; the table counts each
non-DatasetError by type and the dataset-package line that raised it.
"""

import collections
import sys
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[3] / "tests"))
from dataset_variants import ROOT, sources, variants  # noqa: E402

from ctx.dataset import DatasetError, load  # noqa: E402


def cause(error: BaseException) -> str:
    """The innermost frame inside src/ctx/dataset, as module:line and its code."""
    frames = traceback.extract_tb(error.__traceback__)
    inside = [f for f in frames if "/ctx/dataset/" in f.filename]
    last = (inside or frames)[-1]
    return f"{Path(last.filename).name}:{last.lineno} {last.line}"


def outcome(path: Path, causes: collections.Counter[str]) -> str:
    """The variant's outcome; a non-DatasetError also counts its cause."""
    try:
        return "OK " + repr(load(path))
    except DatasetError as error:
        return f"DatasetError {error}"
    except Exception as error:  # noqa: BLE001 - classified, not hidden
        causes[f"{type(error).__name__}\t{cause(error)}"] += 1
        return f"{type(error).__name__} {error}"


def main(dump: Path, table: Path) -> None:
    """Write every outcome to ``dump`` and the cause counts to ``table``."""
    rows: list[str] = []
    causes: collections.Counter[str] = collections.Counter()
    kinds: collections.Counter[str] = collections.Counter()
    for source in sources():
        for j, text in enumerate(variants(source.read_text())):
            path = source.with_name(f"zz_f5_{source.stem}_{j}.toml")
            path.write_text(text)
            try:
                result = outcome(path, causes).replace(str(path), "<variant>")
            finally:
                path.unlink()
            kinds[result.split(" ", 1)[0]] += 1
            rows.append(f"{source.relative_to(ROOT)}#{j}\t{result}")
    dump.write_text("\n".join(rows) + "\n")
    lines = [f"{k}\t{n}" for k, n in sorted(kinds.items())] + [""]
    lines += [f"{n}\t{c}" for c, n in causes.most_common()]
    table.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
