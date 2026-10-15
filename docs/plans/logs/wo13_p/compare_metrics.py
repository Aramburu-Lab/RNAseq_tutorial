"""Package P: Step-1's 92 recorded metrics against a rebuilt root's, value by value (D-046).

Run: python compare_metrics.py <recorded csv> <05_de metrics csv> <06_function metrics csv>
"""

import csv
import sys
from pathlib import Path


def read(path: str, module: str | None = None) -> dict[tuple[str, str], tuple[str, str]]:
    """(module, metric) to (value, detail); a module-less table takes ``module``."""
    lines = Path(path).read_text().splitlines()
    rows = csv.DictReader(line for line in lines if not line.startswith("#"))
    return {(r.get("module") or module, r["metric"]): (r["value"], r["detail"]) for r in rows}


recorded = read(sys.argv[1])
now = {**read(sys.argv[2], "de"), **read(sys.argv[3], "function")}
equal = [k for k in recorded if k in now and now[k][0] == recorded[k][0]]
differ = [k for k in recorded if k in now and now[k][0] != recorded[k][0]]
missing = [k for k in recorded if k not in now]
details = [k for k in recorded if k in now and now[k][1] == recorded[k][1]]
print(
    f"recorded {len(recorded)}; equal {len(equal)}; differ {len(differ)}; "
    f"missing {len(missing)}; details equal {len(details)}"
)
for k in differ + missing:
    print("  ", k, recorded[k], now.get(k))
