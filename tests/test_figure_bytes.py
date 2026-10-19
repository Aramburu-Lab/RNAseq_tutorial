"""656: two renders of one figure are byte-equal, PDF and SVG, for every figure contract.

Each test names its mutation.
"""

from datetime import date

import pytest
from test_viz import KINDS, _contracts

from ctx.io.naming import stem
from ctx.viz.save import write_figure


@pytest.mark.parametrize("kind", KINDS)
def test_two_renders_are_byte_equal(kind, tmp_path):
    """The same contract written twice into two folders gives equal PDFs and SVGs, which carry
    no creation time. Catches the PDF's CreationDate and the SVG's dc:date kept, and the SVG's
    ids hashed with a random salt, each of which alone makes reruns differ."""
    figure_stem = stem("de", f"{kind}_bytes", 1, day=date(2026, 10, 19))
    data = _contracts()[kind]
    files = [write_figure(data, tmp_path / str(n), figure_stem, title="Demo", data_releases={})
             for n in (1, 2)]  # fmt: skip
    for first, second in zip(*files, strict=True):
        if first.suffix in {".pdf", ".svg"}:
            assert first.read_bytes() == second.read_bytes(), first.name
    for path in files[0]:
        if path.suffix == ".pdf":
            assert b"/CreationDate" not in path.read_bytes(), path.name
        if path.suffix == ".svg":
            assert b"<dc:date>" not in path.read_bytes(), path.name
