"""F.6 and R2.5b's planted mutations: each undoes one fix, runs the tests that should catch it.

Run from the worktree root with a clean tree: ``python docs/plans/logs/g6_mutations.py``.
Each mutation is applied, its tests run, and the file restored with ``git checkout``; a
mutation is CAUGHT when the tests fail.
"""

import subprocess
import sys
from pathlib import Path

D = "src/ctx/dataset/"
F6 = ["tests/test_dataset_variants.py"]
R25B = ["tests/test_gather_counts_rsem.py::test_largest_rows_left_out_first",
        "tests/test_gather_counts_history.py::test_gather_reads_the_history_resource"]
MUTATIONS = [
    ("filter table unchecked", D + "checks.py", F6,
     '    table = _table(raw, where, problems)\n', "    table = raw\n"),
    ("filter values unchecked", D + "checks.py", F6,
     'tuple(str(x) for x in _list(v, f"{where}.{k}", problems))', "tuple(str(x) for x in v)"),
    ("column value map unchecked", D + "sources.py", F6,
     "if values is not None and not isinstance(values, Mapping):", "if False:"),
    ("split unchecked", D + "analyses.py", F6,
     '    table = _table(raw, f"{where}.split", problems)\n', "    table = raw\n"),
    ("pca factors unchecked", D + "analyses.py", F6,
     'factors = _sequence(p.get("factors"), f"{where}.factors", problems)',
     'factors = tuple(p["factors"])'),
    ("pca name read directly", D + "analyses.py", F6,
     '        if "name" not in p:\n', '        if p["name"] is None:\n'),
    ("scatters pairs unchecked", D + "analyses.py", F6,
     "if not all(isinstance(pair, list) and len(pair) == 2 for pair in scatters):", "if False:"),
    ("running type unchecked", D + "analyses.py", F6,
     "if not isinstance(running, list) or len(running) != 2:", "if len(running) != 2:"),
    ("cross_check pairs unchecked", D + "sources.py", F6,
     "    pairs = pairs if paired else []\n", ""),
    ("batch bins unchecked", D + "loader.py", F6,
     'bins = {k: _list(v, f"derived.batch.bins.{k}", problems) for k, v in bins.items()}',
     "bins = dict(bins)"),
    ("order values unchecked", D + "loader.py", F6,
     'tuple(_sequence(v, f"order.{k}", problems))', "tuple(v)"),
    ("expected counts unchecked", D + "loader.py", F6,
     "bad = {v: n for v, n in counts.items() if not isinstance(n, int)}", "bad = {}"),
    ("proteo lists unchecked", D + "proteo.py", F6,
     'lists = {k: _sequence(raw.get(k, []), f"proteo.{k}", problems)',
     "lists = {k: raw.get(k, [])"),
    ("proteo class sets unchecked", D + "proteo.py", F6,
     'sets = _sequence(entry.get("sets"), f"proteo.classes.{name}.sets", [])',
     'sets = entry.get("sets")'),
    ("problems raised one by one", D + "loader.py", F6,
     '        raise DatasetError(f"{path}:\\n  " + "\\n  ".join(problems))',
     '        raise DatasetError(f"{path}:\\n  " + problems[0])'),
    ("left-out detail dropped", "src/ctx/gather/counts_run.py", R25B,
     'rows = tuple((i, k, w, d, f"{t:,.0f}") for i, k, w, d, t in',
     'rows = tuple((i, k, w, "", f"{t:,.0f}") for i, k, w, d, t in'),
    ("left-out detail from another row", "src/ctx/gather/counts_run.py", R25B,
     'top.select("source_id", "kind", "reason", "detail", "total").rows())',
     'top.select("source_id", "kind", "reason", pl.col("detail").reverse(), "total").rows())'),
    ("left-out header key reverted", "src/ctx/report/content/gather_counts.toml", R25B,
     '"Reason", "Detail", "Counts held"]', '"Reason", "Note", "Counts held"]'),
]


def run(name: str, path: str, tests: list[str], old: str, new: str) -> str:
    """Apply one mutation, run its tests, restore the file; CAUGHT or SURVIVED."""
    file = Path(path)
    text = file.read_text()
    assert text.count(old) == 1, f"{name}: the text to mutate is not in {path} once"
    file.write_text(text.replace(old, new))
    try:
        done = subprocess.run([sys.executable, "-m", "pytest", "-q", "-x", "-p",
                               "no:cacheprovider", *tests], capture_output=True, text=True)  # fmt: skip
    finally:
        subprocess.run(["git", "checkout", "--", path], check=True)
    return "CAUGHT" if done.returncode else "SURVIVED"


def main() -> None:
    """Run every mutation and print one line each."""
    selected = sys.argv[1:]
    for name, path, tests, old, new in MUTATIONS:
        if not selected or any(s in name for s in selected):
            print(f"{run(name, path, tests, old, new)}\t{name}", flush=True)


if __name__ == "__main__":
    main()
