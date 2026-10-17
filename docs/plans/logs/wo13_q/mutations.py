"""Package Q: planted mutations of the qc GDC door, each run against its tests.

Each mutation replaces one exact line (asserted to occur once), runs the named tests with
no bytecode written, and restores the file whatever happens. A mutation is caught when the
tests fail. Usage: py.sh python docs/plans/logs/wo13_q/mutations.py.
"""

import os
import subprocess
import sys
from pathlib import Path

WT = Path(__file__).resolve().parents[4]
RUN = WT / "src" / "ctx" / "qc" / "run.py"
TESTS = ["tests/test_qc_gdc_door.py"]
MUTATIONS = {
    "no door": (
        '    sheet = gdc.chosen(sheet) if "chosen" in sheet.columns else sheet\n',
        "",
    ),
    "counts not cut to the door": (
        '        gathered["counts"] = pl.read_parquet(inputs["counts"], '
        'columns=["gene_id", *samples])\n',
        '        gathered["counts"] = pl.read_parquet(inputs["counts"])\n',
    ),
    "door before the read, on the raw sheet's samples": (
        '    samples = sheet.unique("sample", keep="first", maintain_order=True)'
        '["sample"].to_list()\n',
        '    samples = pl.read_parquet(inputs["samplesheet"])["sample"].to_list()\n',
    ),
}


def main() -> int:
    """Run every mutation; print caught or SURVIVED; non-zero if any survived."""
    original = RUN.read_text()
    env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}
    survived = 0
    for name, (old, new) in MUTATIONS.items():
        assert original.count(old) == 1, name
        try:
            RUN.write_text(original.replace(old, new))
            command = [sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider",
                       *TESTS]  # fmt: skip
            done = subprocess.run(command, cwd=WT, env=env, capture_output=True, text=True,
                                  check=False)  # fmt: skip
        finally:
            RUN.write_text(original)
        caught = done.returncode != 0
        survived += not caught
        print(f"{'caught' if caught else 'SURVIVED'}: {name}")
    return survived


if __name__ == "__main__":
    sys.exit(main())
