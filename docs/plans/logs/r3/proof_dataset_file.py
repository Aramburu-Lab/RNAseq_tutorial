"""R3 proof, the student's edits to the template copied as the guide copies it.

Its slug and label set, the [proteo] block pasted with its EDIT lines answered for a TMT
MaxQuant table, and the [[resources]] entry ``prepare fetch`` printed pasted after it.

Run through env.sh: python docs/plans/logs/r3/proof_dataset_file.py <dataset file> <fetch log>
"""

import sys
from pathlib import Path

from ctx import pipeline

path, fetched = Path(sys.argv[1]), Path(sys.argv[2])
text = path.read_text()
for old, new in (
    ('slug = "mine"', 'slug = "r3_pxd038437"'),
    ('label = "My study"', 'label = "Multiple myeloma, Ramberger 2024 (R3 proof)"'),
):
    assert text.count(old) == 1, old
    text = text.replace(old, new)
block = pipeline.block("[proteo]", slug="r3_pxd038437")
answers = (
    ('member = "proteinGroups.txt"\n', ""),
    ("tracks = []", 'tracks = ["status"]'),
    ("genes = []", 'genes = ["CD38", "IRF4"]'),
    ('quantity = "LFQ intensity"', 'quantity = "Reporter intensity corrected"'),
    ('acquisition = "label_free"', 'acquisition = "tmt"'),
    ('batch = "batch"', 'batch = "batch"\nreference = "reference"'),
)
for old, new in answers:
    assert block.count(old) == 1, old
    block = block.replace(old, new)
entry = "".join(line for line in fetched.read_text().splitlines(keepends=True)
                if not line.startswith("exit"))  # fmt: skip
path.write_text(text + "\n" + block + "\n" + entry)
print(path.read_text()[-1600:])
