"""Plant one mutation per new rule of wo13-v, run the tests that should catch it, restore.

Each mutation is an exact replacement asserted to match once; a mutation survives when the
tests pass with it in place. Run through env.sh so the worktree's src is imported.
"""

import subprocess
import sys
from pathlib import Path

WT = Path(__file__).resolve().parents[4]
PAGE = WT / "src" / "ctx" / "clinical" / "page.py"
TEXT = WT / "src" / "ctx" / "report" / "content" / "clinical.toml"
TESTS = ["tests/test_clinical_page.py", "-k", "shared or note or sex or own_endpoint"]
MUTATIONS = [
    (
        "Question names breast again",
        TEXT,
        'prompt = "Why can [[overall survival]] be a weak endpoint where most patients outlive '
        "TCGA's follow-up?\"",
        'prompt = "Why does TCGA-CDR advise against [[overall survival]] as the endpoint for '
        'breast cancer?"',
    ),
    (
        "ER lead names breast again",
        TEXT,
        "ER-negative tumours progress sooner",
        "ER-negative breast cancer progresses sooner",
    ),
    (
        "counts of the dataset file's endpoint, not the fitted one",
        PAGE,
        'counts = clinical.published.get("table3", {}).get(ctx.endpoint)',
        'counts = clinical.published.get("table3", {}).get(clinical.endpoint)',
    ),
    (
        "--endpoint called the dataset file's",
        PAGE,
        "if ctx.endpoint == clinical.endpoint",
        "if ctx.endpoint != clinical.endpoint",
    ),
    ("refusals left out", PAGE, "    if against:\n", "    if not against:\n"),
    (
        "cautions left out",
        PAGE,
        "if caution else []  # fmt: skip",
        "if not caution else []  # fmt: skip",
    ),
    (
        "a misspelt mark passes",
        PAGE,
        "if marks[e] not in USES}",
        "if marks[e] not in (*USES, 'yes')}",
    ),
    (
        "an empty mark read as copied",
        PAGE,
        "if not all(marks.get(e) for e in ENDPOINTS)",
        "if not all(e in marks for e in ENDPOINTS)",
    ),
    (
        "sex sentence on a cohort that tests sex",
        PAGE,
        '    if "sex" in variables(ctx):\n        return ()\n',
        "",
    ),
    ("the frequent sex counted", PAGE, "fewer = min(SEXES", "fewer = max(SEXES"),
    ("a tie printed", PAGE, "or counts[fewer] == max(counts.values())", ""),
    (
        "the endpoint note built and dropped",
        PAGE,
        "(flow, endpoints, *endpoint_note(ctx), missing",
        "(flow, endpoints, missing",
    ),
]


def run(name: str, path: Path, old: str, new: str) -> bool:
    """True when the tests fail with the mutation in place."""
    body = path.read_text()
    assert body.count(old) == 1, (name, old)
    path.write_text(body.replace(old, new))
    try:
        command = [sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider", *TESTS]
        done = subprocess.run(command, cwd=WT, capture_output=True, text=True, check=False)
    finally:
        path.write_text(body)
    caught = done.returncode != 0
    print(f"{'caught  ' if caught else 'SURVIVED'} {name}: {done.stdout.strip().splitlines()[-1]}")
    return caught


results = [run(*m) for m in MUTATIONS]
print(f"{sum(results)} of {len(results)} caught")
sys.exit(0 if all(results) else 1)
