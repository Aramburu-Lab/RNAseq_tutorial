"""F.6: a broken dataset file parses or raises one DatasetError, never a traceback.

Before F.6, 407 of the 7,800 variants tests/dataset_variants.py makes raised AttributeError,
TypeError, ValueError or KeyError (docs/plans/logs/2026-10-19_g6_classification.log). F.6a
refuses a string for a list and true for a count: 72 variants that parsed are now refused, 27
for a string and 45 for true (docs/plans/logs/2026-10-20_y_f6a_classification.log).
"""

from pathlib import Path

import pytest
from dataset_variants import SWAPS, sources, variants

from ctx.dataset import DatasetError, load, loader

# Each key line retyped to an integer and to a string (SWAPS 1 and 3), and every file's
# all-lines-retyped variants: half of the 7,800, about 2 s; deleting a line and retyping it
# to a boolean are the full run's (docs/plans/logs/g6_classify.py). Every third variant
# missed the one scatters and one bins line, and their mutations survived.
KEPT = (1, 3)
SYNTHETIC = Path(__file__).parent / "fixtures" / "datasets" / "synthetic.toml"
CONTRAST = (
    '[[contrasts]]\nname = "c"\nvariable = "condition"\nreference = "control"\n'
    'tested = "case"\nn = [7, 7]\n'
)


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
            if j % len(SWAPS) not in KEPT and j < len(texts) - 4:
                continue
            found += 1
            if (result := in_memory(text)) not in ("parsed", "listed"):
                failures.append(f"{source.name}#{j}: {result[:200]}")
    assert len(sources()) >= 22
    assert found > 3500
    assert not failures, f"{len(failures)} of {found}:\n" + "\n".join(failures[:20])


def test_crashing_kinds_listed_together(tmp_path):
    """Four values the old parsers crashed on, in one file, are four lines of one error.
    Catches a fix that raises on the first instead of listing every problem."""
    text = SYNTHETIC.read_text()
    text = text.replace('values = { "Control" = "control", "Case" = "case" }', "values = 7")
    text = text.replace('time_point = ["", "fu2", "fu10"]', "time_point = true")
    text = text.replace("control = 7", 'control = "x"')
    text += CONTRAST + "filters = { condition = 7 }\n"
    path = tmp_path / "broken.toml"
    path.write_text(text)
    with pytest.raises(DatasetError) as caught:
        load(path)
    assert str(caught.value).splitlines() == [
        f"{path}:",
        "  columns.condition.values: expected a table, found 7",
        "  order.time_point: expected a list, found True",
        "  expect.counts.condition.control: expected int, found 'x'",
        "  contrasts[0].filters.condition: expected a list, found 7",
    ]


@pytest.mark.parametrize(("snippet", "problem"), [
    ('[de]\nfeatured = "c"\nscatters = [["c"]]\n',
     "de.scatters: expected [[contrast, contrast], ...], found [['c']]"),
    ('[de]\nfeatured = "c"\n[[de.pca]]\nfactors = ["condition"]\n',
     "de.pca[0].name: expected str, found None"),
    ("[derived.batch]\nbins = { short = 50 }\n",
     "derived.batch.bins.short: expected a list, found 50"),
])  # fmt: skip
def test_nested_values_listed(tmp_path, snippet, problem):
    """A value inside a table that the variants retype only whole is checked too. Catches a
    one-item scatter, a PCA with no name and a bin of one length crashing the parser."""
    path = tmp_path / "broken.toml"
    path.write_text(SYNTHETIC.read_text() + CONTRAST + snippet)
    with pytest.raises(DatasetError) as caught:
        load(path)
    assert str(caught.value).splitlines() == [f"{path}:", f"  {problem}"]


@pytest.mark.parametrize(("old", "new", "problem"), [
    ("", '[de]\nfeatured = "c"\n[[de.pca]]\nname = "p"\nfactors = "condition"\n',
     "de.pca[0].factors: expected a list, found 'condition'"),
    ("", '[de]\nfeatured = "c"\nbox_genes = "ACTB"\n',
     "de.box_genes: expected a list, found 'ACTB'"),
    ('time_point = ["", "fu2", "fu10"]', 'time_point = "fu2"',
     "order.time_point: expected a list, found 'fu2'"),
    ("control = 7", "control = true", "expect.counts.condition.control: expected int, found True"),
    ("n = [7, 7]", "n = [true, 7]",
     "contrasts[0].n: expected [reference n, tested n], found [True, 7]"),
])  # fmt: skip
def test_a_string_for_a_list_and_true_for_a_count_are_refused(tmp_path, old, new, problem):
    """F.6a (Rytis 2026-10-20): a string where a list belongs is refused, not read as its
    characters (``factors = "condition"`` was 9 one-letter factors), and true is no count,
    though isinstance takes a bool for an int. Catches either let through again."""
    text = SYNTHETIC.read_text() + CONTRAST
    if old:
        assert text.count(old) == 1, old
        text = text.replace(old, new)
    else:
        text += new
    path = tmp_path / "broken.toml"
    path.write_text(text)
    with pytest.raises(DatasetError) as caught:
        load(path)
    assert str(caught.value).splitlines() == [f"{path}:", f"  {problem}"]
