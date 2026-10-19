"""List every name the repository imports from ctx.dataset (F.5)."""

import ast
import collections
import re
import subprocess
from pathlib import Path

files = subprocess.run(
    ["git", "ls-files", "*.py"], capture_output=True, text=True, check=True
).stdout.split()
names = collections.defaultdict(set)
attrs = collections.defaultdict(set)
strings = []
for f in files:
    if f == "src/ctx/dataset.py" or f.startswith("src/ctx/dataset/"):
        continue
    src = Path(f).read_text()
    tree = ast.parse(src)
    pkg = f.removeprefix("src/").removesuffix(".py").split("/")
    if pkg[-1] == "__init__":
        pkg = pkg[:-1] + ["x"]
    mods = set()
    for n in ast.walk(tree):
        if isinstance(n, ast.ImportFrom):
            mod = n.module or ""
            full = ".".join(pkg[: -n.level] + ([mod] if mod else [])) if n.level else mod
            if full == "ctx.dataset":
                for a in n.names:
                    names[a.name].add(f)
            elif full == "ctx":
                mods |= {a.asname or "dataset" for a in n.names if a.name == "dataset"}
        elif isinstance(n, ast.Import):
            mods |= {a.asname or "ctx.dataset" for a in n.names if a.name == "ctx.dataset"}
    for m in mods:
        for x in re.findall(re.escape(m) + r"\.(\w+)", src):
            attrs[x].add(f)
    strings += [(f, x) for x in re.findall(r"[\"']ctx\.dataset\.?(\w*)", src)]
for k in sorted(names):
    print(k, len(names[k]))
print("count", len(names))
print("ATTR", {k: sorted(v) for k, v in attrs.items()})
print("STR", strings)
