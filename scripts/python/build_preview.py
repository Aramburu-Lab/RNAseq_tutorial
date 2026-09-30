#!/usr/bin/env python3
# /// codebase-meta
# exempt: instructor-side; builds the course preview, every module on made-up data until it is built
# ///
r"""Build the course preview: one page per module M0 to M10 and the hub, on one made-up dataset.

Built modules run as a student runs them, through ``ctx``, on the test suite's synthetic
design (tests/synth.py): counts with planted differential genes, a batch and three planted QC
problems. A module not built yet draws its page from ``scripts/python/preview_<module>.py``,
whose ``page(root, dataset, day)`` writes ``outputs/NN_<module>/<day>_<module>_<dataset>_
overview_v1.html`` from made-up data and the page text in ``content/<module>.toml``. Every
preview page says, in its first result, that its numbers are made up.

Usage:
    pixi run python scripts/python/build_preview.py --root ~/ctx_preview/synthetic
    pixi run python scripts/python/build_preview.py --root <dir> --only process
"""

from __future__ import annotations

import argparse
import importlib
import subprocess
import sys
import time
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(REPO / "tests"), str(Path(__file__).parent)]
CTX = str(Path(sys.executable).parent / "ctx")
DATASET = "desynth"
# Course order; each is a MODULES token with no command yet (D-044).
PREVIEWS = ("process", "quant", "proteo", "ncrna", "genecard", "clinical")
# The first result of every preview page: text.word("ui", "preview.note").
NOTE = ("ui", "preview.note")
# The course's batch names (D-043) in place of the test data's short and long.
BATCH = {"short": "Reads ≤ 76 nt (min 35 nt)", "long": "Reads > 76 nt (max 500 nt)"}


def ctx(root: Path, *args: str) -> None:
    """Run one ctx command in root, timed; stop the build on a failure."""
    start = time.perf_counter()
    done = subprocess.run([CTX, "--root", str(root), *args], cwd=root, capture_output=True,
                          text=True, check=False)  # fmt: skip
    print(f"{time.perf_counter() - start:6.1f} s  exit {done.returncode}  ctx {' '.join(args)}")
    if done.returncode:
        sys.exit(done.stderr[-3000:])


def synthetic(root: Path) -> Path:
    """The made-up gather outputs, gene types and gene sets; return the dataset file."""
    import polars as pl  # noqa: PLC0415
    from synth import write_gene_sets  # noqa: PLC0415
    from test_qc_run import plant, with_annotation  # noqa: PLC0415

    from ctx.dataset import load  # noqa: PLC0415

    plant(root)
    spec = with_annotation(root, root)
    sheet = next((root / "outputs" / "01_gather").glob("*_samplesheet_v1.parquet"))
    frame = pl.read_parquet(sheet)
    frame.with_columns(batch=pl.col("batch").replace(BATCH)).write_parquet(sheet)
    text = spec.read_text()
    for old, new in BATCH.items():
        text = text.replace(f'batch = ["{old}"]', f'batch = ["{new}"]')
    spec.write_text(text)
    write_gene_sets(root / load(spec).resource("sets").dest)
    return spec


def main() -> None:
    """Build every page, or only the preview modules named."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--only", nargs="*", choices=PREVIEWS, help="preview modules only")
    args = parser.parse_args()
    root, day = args.root.expanduser().resolve(), date.today()
    if args.only is None:
        if root.exists():
            sys.exit(f"{root} exists; move it aside first")
        root.mkdir(parents=True)
        spec = synthetic(root)
        ctx(root, "setup")
        for command in ("qc", "de", "function"):
            ctx(root, command, "--dataset", DATASET, "--dataset-file", str(spec))
    for module in args.only or PREVIEWS:
        start = time.perf_counter()
        page = importlib.import_module(f"preview_{module}").page(root, DATASET, day)
        print(f"{time.perf_counter() - start:6.1f} s  preview {module}: {page.name}")
    if args.only is None:
        ctx(root, "report", "--dataset", DATASET, "--dataset-file", str(spec))


if __name__ == "__main__":
    main()
