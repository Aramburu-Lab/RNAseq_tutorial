"""F.6: a broken dataset file parses or raises one DatasetError, never a traceback.

Before F.6, 407 of the 7,800 variants tests/dataset_variants.py makes raised AttributeError,
TypeError, ValueError or KeyError (docs/plans/logs/2026-10-19_g6_classification.log).
"""

from pathlib import Path

import pytest
from dataset_variants import sources, variants

from ctx.dataset import DatasetError, load, loader

# Every third variant, and every file's all-lines-retyped variants: about 2,600 of 7,800,
# under 2 s; the full run is docs/plans/logs/g6_classify.py.
EVERY = 3
SYNTHETIC = Path(__file__).parent / "fixtures" / "datasets" / "synthetic.toml"


@pytest.fixture
def in_memory(monkeypatch):
    """Load text as a file at a path never written: writing 7,800 files took a minute."""
    path, current, real = Path("never_written") / "variant.toml", {}, Path.read_text

    def read_text(self: Path, *args, **kwargs) -> str:
        return current["text"] if self == path else real(self, *args, **kwargs)

    monkeypatch.setattr(Path, "read_text", read_text)
    # load records the path it read per slug; the variants' slugs are real cohorts'.
    monkeypatch.setattr(loader, "_LOADED", {})

    def outcome(text: str) -> str:
        current["text"] = text
        try:
            load(path)
        except DatasetError as error:
            return "listed" if str(error).startswith(f"{path}:\n  ") else f"unnamed {error}"
        except Exception as error:  # noqa: BLE001 - the outcome under test
            return f"{type(error).__name__} {error}"
        return "parsed"

    return outcome


def test_every_sampled_variant_parses_or_lists(in_memory):
    """Each sampled variant is a Dataset or a DatasetError opening with the file. Catches any
    parser reading a value of the wrong type unchecked (columns.values = 7: AttributeError)."""
    found, failures = 0, []
    for source in sources():
        texts = variants(source.read_text())
        for j, text in enumerate(texts):
            if j % EVERY and j < len(texts) - 4:
                continue
            found += 1
            if (result := in_memory(text)) not in ("parsed", "listed"):
                failures.append(f"{source.name}#{j}: {result[:200]}")
    assert len(sources()) >= 22
    assert found > 2000
    assert not failures, f"{len(failures)} of {found}:\n" + "\n".join(failures[:20])


def test_crashing_kinds_listed_together(tmp_path):
    """Three values the old parsers crashed on, in one file, are three lines of one error.
    Catches a fix that raises on the first instead of listing every problem."""
    text = SYNTHETIC.read_text()
    text = text.replace('values = { "Control" = "control", "Case" = "case" }', "values = 7")
    text = text.replace('time_point = ["", "fu2", "fu10"]', "time_point = true")
    text = text.replace("control = 7", 'control = "x"')
    path = tmp_path / "broken.toml"
    path.write_text(text)
    with pytest.raises(DatasetError) as caught:
        load(path)
    assert str(caught.value).splitlines() == [
        f"{path}:",
        "  columns.condition.values: expected a table, found 7",
        "  order.time_point: expected a list, found True",
        "  expect.counts.condition.control: expected int, found 'x'",
    ]
