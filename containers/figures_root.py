# /// codebase-meta
# exempt: >-
#   build tooling; writes the fixture root the figures image's proof redraws from, outside outputs/
# ///
"""Write a fixture ctx root and run ``ctx de`` on it: the volcano the figures proof redraws.

    python containers/figures_root.py <empty folder>

The synthetic design is tests/synth.py's DE_GROUPS on tests/fixtures/datasets/de_synthetic.toml,
the same root the test suite's ``pipeline`` fixture builds, without its wheel.
"""

import subprocess
import sys
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DATASET = REPO / "tests" / "fixtures" / "datasets" / "de_synthetic.toml"


def main(root: Path) -> None:
    """Write the gather outputs and the gene sets into ``root``, then run ``ctx de``.

    Raises:
        SystemExit: If ``root`` is not empty, or ``ctx de`` fails.
    """
    if root.exists() and any(root.iterdir()):
        raise SystemExit(f"{root}: not empty")
    sys.path.insert(0, str(REPO / "tests"))
    from synth import (  # noqa: PLC0415 - tests/ is on the path only from here
        DE_GROUPS,
        design_counts,
        write_gather_outputs,
        write_gene_sets,
    )

    from ctx.dataset import load  # noqa: PLC0415 - after the path

    root.mkdir(parents=True, exist_ok=True)
    write_gather_outputs(root, "desynth", design_counts(DE_GROUPS, batch_lfc=2.0), date.today())
    write_gene_sets(root / load(DATASET).resource("sets").dest)
    command = ["ctx", "--root", str(root), "de", "--dataset", "desynth",
               "--dataset-file", str(DATASET)]  # fmt: skip
    if subprocess.run(command, cwd=root, check=False).returncode != 0:
        raise SystemExit("ctx de failed")


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
