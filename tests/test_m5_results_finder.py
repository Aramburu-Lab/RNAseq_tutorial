"""Every reader of M5's results finds them through ``ctx.de.results`` alone (R1.3d).

The finder is replaced by a stub that answers from a folder no reader would build; a reader that
spells ``<contrast>_results`` itself looks in ``05_de`` and misses it.
"""

from pathlib import Path

import polars as pl
import pytest
from conftest import DE_SPEC, ROOT

from ctx.de import results as m5
from ctx.figures import kinds, record
from ctx.genecard import cards

DATASET_FILES = [DE_SPEC, *sorted((ROOT / "conf" / "datasets").glob("*.toml"))]


def _stub(asked: list[str], folder: Path):
    def results(_root: Path, _slug: str, contrast: str) -> Path:
        asked.append(contrast)
        return folder / f"{contrast}.parquet"

    return results


@pytest.mark.parametrize("path", DATASET_FILES, ids=[p.stem for p in DATASET_FILES])
def test_every_figure_kind_finds_m5_results_through_the_finder(path, tmp_path, monkeypatch):
    """Catches a kind that names its input ``("de", f"{contrast}_results")``: input_paths then
    finds it with ``newest`` in 05_de, and the stub is never asked."""
    from ctx.dataset import load  # noqa: PLC0415

    spec = load(path)
    asked: list[str] = []
    monkeypatch.setattr(m5, "results", _stub(asked, tmp_path / "elsewhere"))
    monkeypatch.setattr(record, "newest", lambda folder, _slug, what: folder / what)
    monkeypatch.setattr(record, "_fetched", lambda _root, _spec, name: tmp_path / name)
    monkeypatch.setattr(record.pl, "read_parquet", lambda _p: {"registry": ["registry.toml"]})
    checked = 0
    for k in kinds.KINDS:
        if not kinds.offered(k, spec) or k.parsed is not None:
            continue
        for params in k.made(spec):
            asked.clear()
            found = record.input_paths(k, params, spec, tmp_path)
            from_m5 = {n for n, p in found.items() if p.parent == tmp_path / "elsewhere"}
            spelled = {n for n, p in found.items() if p.name.endswith("_results")
                       and p.parent.name == "05_de"}  # fmt: skip
            assert not spelled, (k.module, k.figure, spelled)
            assert len(from_m5) == len(asked), (k.module, k.figure)
            checked += len(from_m5)
    assert checked > 0 or spec.de is None  # a file with [de] has M5 kinds that read results


def test_genecard_finds_m5_results_through_the_finder(tmp_path, monkeypatch):
    """Catches ``cards.results`` building ``<contrast>_results`` in 05_de itself."""
    folder = tmp_path / "elsewhere"
    folder.mkdir()
    pl.DataFrame({"gene_id": ["ENSG00000106236.4"], "log2FoldChange": [1.0],
                  "log2FoldChange_shrunk": [0.5], "padj": [0.01],
                  "call": ["up"]}).write_parquet(folder / "a_case.parquet")  # fmt: skip
    asked: list[str] = []
    monkeypatch.setattr(m5, "results", _stub(asked, folder))
    frame = cards.results(tmp_path, "any", "a_case")
    assert asked == ["a_case"]
    assert frame["number"].to_list() == [106236]
