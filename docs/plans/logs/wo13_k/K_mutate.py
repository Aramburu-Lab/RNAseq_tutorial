"""Plant one mutation per K.9 and K.10 rule in build_kit.py, run the kit tests, restore.

Run through K_mutate.sh. Each mutation must turn at least one test red; the file is restored
byte for byte after each, and the script says if it was not.
"""

import subprocess
import sys
from pathlib import Path

W = Path(__file__).resolve().parents[4]
TARGET = W / "scripts" / "python" / "build_kit.py"
MUTATIONS = {
    "no mapping of a committed dataset file": (
        'mapped = (codebase or {}).get((recorded or {}).get(Path(word).name, ""))',
        "mapped = None",
    ),
    "mapping by path, not by bytes": (
        'mapped = (codebase or {}).get((recorded or {}).get(Path(word).name, ""))',
        'mapped = f"conf/datasets/{Path(word).name}" if Path(word).parent.name == "datasets"'
        " else None",
    ),
    "only the newest commit read": ('_git("log", "--full-history"', '_git("log", "-n", "1"'),
    "merges skipped": ('"--full-history", "-m",', ""),
    "root sentence dropped": ('    if entry.get("builder_root"):', "    if False:"),
    "own-file sentence dropped": ('    if entry.get("builder_dataset_file"):', "    if False:"),
    "note left out of README.md": ('note=(notes or {}).get(folder, "")', 'note=""'),
    "note left out of MANIFEST.json": ('"note": (notes or {}).get(f)', '"note": None'),
    "note folder unchecked": ("        if folder not in folders:", "        if False:"),
    "note prose lint skipped": (
        "problems += [*prose.check_text(where, note, None), *prose.asks(where, note)]",
        "problems += []",
    ),
}


def main() -> None:
    """Plant each mutation, count the failing tests, restore the file."""
    original = TARGET.read_bytes()
    survived = []
    try:
        for name, (old, new) in MUTATIONS.items():
            source = original.decode()
            assert source.count(old) == 1, f"{name}: {old!r} found {source.count(old)} times"
            TARGET.write_text(source.replace(old, new))
            done = subprocess.run(
                [
                    "bash",
                    str(W / "docs/plans/logs/wo13_k/K_pytest.sh"),
                    "-x",
                    "-k",
                    "rebuild_command or committed_version or note or gdc_count_matrix",
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            last = done.stdout.strip().splitlines()[-1]
            print(f"{name}: {'caught' if done.returncode else 'SURVIVED'} ({last})")
            if not done.returncode:
                survived.append(name)
    finally:
        TARGET.write_bytes(original)
    print("restored" if TARGET.read_bytes() == original else "NOT RESTORED")
    sys.exit(1 if survived else 0)


if __name__ == "__main__":
    main()
