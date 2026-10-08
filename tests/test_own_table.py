"""U.1: a study without an SDRF, curated from the student's own table beside the dataset file.

The table is written at test time from the metadata-only fetchngs fixture's own runs, saved as
a spreadsheet saves it: a byte-order mark and Windows line ends. Each test names its mutation.
"""

import csv
import hashlib
import io
from pathlib import Path

import polars as pl
import pytest
from conftest import ctx_run

from ctx.dataset import DatasetError, load
from ctx.gather.curation import CurationError, read_table
from ctx.gather.fetchngs import read
from ctx.gather.harmonise import curated_runs

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "fetchngs" / "metadata_only"
SPEC = """[dataset]
slug = "mine"

[guard]
ids = []
patterns = []

[sdrf]
file = "{table}"
join = ["run_accession", "run_accession"]

[columns.donor]
source = "donor"
meaning = "Individual."

[columns.cell_type]
source = "cell type"
meaning = "Cell type."

[columns.condition]
source = "condition"
meaning = "Condition."
values = {{ "Control" = "control", "Case" = "case" }}
"""


@pytest.fixture(scope="module")
def runs():
    return read(FIXTURE).runs


def table(runs: pl.DataFrame, delimiter: str = ",", drop: int = 0) -> bytes:
    """One row per run, but the last ``drop``; a byte-order mark and CRLF line ends."""
    buffer = io.StringIO()
    out = csv.writer(buffer, delimiter=delimiter, lineterminator="\r\n")
    out.writerow(["run_accession", "donor", "cell type", "condition"])
    names = runs["run_accession"].to_list()
    for i, run in enumerate(names[: len(names) - drop]):
        kind, condition = "type, A" if i % 2 else "type B", "Control" if i < 7 else "Case"
        out.writerow([run, f"d{i}", kind, condition])
    return chr(0xFEFF).encode() + buffer.getvalue().encode()


def dataset(folder: Path, name: str, data: bytes) -> Path:
    """The dataset file and its table, in a folder of their own."""
    folder.mkdir(parents=True, exist_ok=True)
    (folder / name).write_bytes(data)
    spec = folder / "mine.toml"
    spec.write_text(SPEC.format(table=name))
    return spec


@pytest.mark.parametrize(("name", "delimiter"), [("mine.csv", ","), ("mine.tsv", "\t")])
def test_own_table_curates(tmp_path, runs, name, delimiter, monkeypatch):
    """Curated from a CSV (a quoted comma in a cell) or a TSV, from another working folder.
    Catches the table read from the working folder, tabs forced on a CSV, the byte-order mark
    kept (the join header is then not found), and CRLF kept in the last column."""
    spec = load(dataset(tmp_path / "conf", name, table(runs, delimiter)))
    monkeypatch.chdir(tmp_path)
    (tmp_path / name).write_text("run_accession\nnot this file\n")
    releases: dict[str, str] = {}
    frame = curated_runs(tmp_path / "root", spec, read(FIXTURE), releases)
    assert frame.height == runs.height
    assert set(frame["cell_type"]) == {"type, A", "type B"}
    counted = frame["condition"].value_counts().sort("condition").rows()
    assert counted == [("case", 7), ("control", 7)]
    data = (tmp_path / "conf" / name).read_bytes()
    assert releases[f"sha256:{name}"] == hashlib.sha256(data).hexdigest()


def test_own_table_missing_run(tmp_path, runs):
    """A run the table does not list stops curation, naming it. Catches the reverse join
    dropped for own tables."""
    spec = load(dataset(tmp_path, "mine.csv", table(runs, drop=1)))
    with pytest.raises(CurationError) as caught:
        curated_runs(tmp_path / "root", spec, read(FIXTURE), {})
    last = runs["run_accession"][-1]
    assert [(c, k) for c, k, _ in caught.value.problems] == [("join", last)]


def test_own_table_absent(tmp_path, runs):
    """No table beside the dataset file: the message names the path and where to start it."""
    spec = load(dataset(tmp_path, "other.csv", b"x"))
    (tmp_path / "other.csv").unlink()
    with pytest.raises(DatasetError, match=r"other\.csv is missing; write it beside mine\.toml"):
        curated_runs(tmp_path / "root", spec, read(FIXTURE), {})


@pytest.mark.parametrize(
    ("sdrf", "said"),
    [('file = "t.csv"\nresource = "sdrf"', "not both"), ("", "not both"),
     ('file = "t.xlsx"', "read by its ending")],
)  # fmt: skip
def test_own_table_schema(tmp_path, sdrf, said):
    """The loader refuses both routes at once, neither, and an ending it cannot read. Catches
    a file and a resource both accepted, which curates from one silently."""
    path = tmp_path / "x.toml"
    path.write_text(f'[dataset]\nslug = "x"\n\n[sdrf]\n{sdrf}\njoin = ["a", "b"]\n')
    with pytest.raises(DatasetError, match=said):
        load(path)


def test_read_table_keeps_quoted_line_ends(tmp_path):
    """A quoted cell keeps its own line end, and the row ends go. Catches a CR stripped by a
    global replace, which also rewrote text inside cells."""
    path = tmp_path / "t.csv"
    path.write_bytes(chr(0xFEFF).encode() + b'a,b\r\n"x\r\ny",z\r\n')
    body, digest = read_table(path)
    assert list(csv.reader(io.StringIO(body, newline=""), delimiter="\t")) == [
        ["a", "b"],
        ["x\r\ny", "z"],
    ]
    assert digest == hashlib.sha256(path.read_bytes()).hexdigest()


def test_gather_reads_own_table(tmp_path, runs):
    """ctx gather end to end, the dataset file given from another folder: the sheet carries the
    curated values and its README the table's SHA-256 as saved. Catches the hash dropped from
    the README, and the table found by the working folder."""
    spec = dataset(tmp_path / "conf", "mine.csv", table(runs))
    root = tmp_path / "root"
    root.mkdir()
    done = ctx_run("--root", str(root), "gather", "--dataset", "mine", "--fetchngs-dir",
                   str(FIXTURE), "--dataset-file", str(spec), cwd=root)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    (sheet,) = (root / "outputs" / "01_gather").glob("*_gather_mine_samplesheet_v1.csv")
    assert set(pl.read_csv(sheet)["donor"]) == {f"d{i}" for i in range(runs.height)}
    readme = sheet.with_name(sheet.stem + "_README.md").read_text()
    assert hashlib.sha256((tmp_path / "conf" / "mine.csv").read_bytes()).hexdigest() in readme
