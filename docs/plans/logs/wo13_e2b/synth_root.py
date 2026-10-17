"""Write the pipeline fixture's gather outputs (tests/conftest.py pipeline) into a new root.

python synth_root.py <root>; then run each step as
python -m ctx --root <root> de --dataset desynth \
    --dataset-file tests/fixtures/datasets/de_synthetic.toml
"""

import sys
from datetime import date
from pathlib import Path

WT = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(WT / "tests"), str(WT / "src")]

from synth import DE_GROUPS, design_counts, write_gather_outputs, write_gene_sets  # noqa: E402

from ctx.dataset import load  # noqa: E402

SPEC = WT / "tests" / "fixtures" / "datasets" / "de_synthetic.toml"
root = Path(sys.argv[1])
write_gather_outputs(root, "desynth", design_counts(DE_GROUPS, batch_lfc=2.0), date.today())
write_gene_sets(root / load(SPEC).resource("sets").dest)
print("written", root)
