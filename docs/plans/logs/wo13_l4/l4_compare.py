"""L.4: compare every table and README the before and after M7 builds wrote."""

import filecmp
import sys
from pathlib import Path

before, after = (Path(a) / "outputs" / "07_proteo" for a in sys.argv[1:3])
for pattern in ("*.csv", "*.parquet", "*_README.md", "*.svg"):
    names = sorted(p.name for p in before.glob(pattern))
    differ = [n for n in names if not filecmp.cmp(before / n, after / n, shallow=False)]
    print(f"{pattern}: {len(names)} files, {len(differ)} differ: {differ}")
