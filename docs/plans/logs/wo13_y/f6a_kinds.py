"""wo13-y item 5: the F.6a changes of g6_compare.py's output, counted by the problem added.

Usage: python f6a_kinds.py <g6_compare output>. Each changed variant's added problem line is
reduced to its key with indices and values dropped, then counted, string and bool apart.
"""

import re
import sys
from collections import Counter
from pathlib import Path

text = Path(sys.argv[1]).read_text()
found: Counter[str] = Counter()
for block in re.split(r"\n(?=CHANGED|ADDS)", text):
    if block.startswith("CHANGED"):
        added = [x.strip() for x in block.split("after:", 1)[1].splitlines()[1:] if x.strip()]
    elif block.startswith("ADDS"):
        added = [x.strip() for x in block.split("\t")[2].split(" | ")]
    else:
        continue
    for line in added:
        key, _, value = line.partition(": ")
        kind = "true" if value.endswith("True") or "True" in value else "string"
        key = re.sub(r"\[\d+\]", "[i]", key)
        key = re.sub(r"^(expect\.counts|order|proteo\.classes)\.[^.:]+", r"\1.<k>", key)
        key = re.sub(r"^(expect\.counts\.<k>)\.[^:]+", r"\1.<v>", key)
        found[f"{kind}\t{key}"] += 1
for name, n in sorted(found.items()):
    print(f"{n}\t{name}")
print(f"{sum(found.values())}\tproblems added")
