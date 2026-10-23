"""Output naming (one builder, one checker) and the README companion."""

import hashlib
import re
import subprocess
import sys
import textwrap
from datetime import date
from pathlib import Path

import polars as pl
import pytest
from test_packaging import build_wheel, run_from_wheel

from ctx.io.naming import MODULES, check, description, figure_files, lint, stem
from ctx.io.readme_companion import write_table

ROOT = Path(__file__).resolve().parents[1]
DAY = date(2026, 9, 26)
COLUMNS = {
    "gene_id": ("none", "Gene identifier."),
    "count": ("reads", "Reads assigned to the gene."),
    "score": ("log2 ratio", "A planted Float32 value."),
}


def _frame() -> pl.DataFrame:
    return pl.DataFrame(
        {"gene_id": ["g1", "g2", "g3"], "count": [5, 0, 12], "score": [0.5, -1.25, 2.0]},
        schema={"gene_id": pl.String, "count": pl.Int32, "score": pl.Float32},
    )


def test_stem_token_order():
    assert stem("de", "volcano", 1, day=DAY) == "2026-09-26_de_volcano_v1"
    names = [
        name
        for token, _ in MODULES
        for desc in ("volcano", "top_50_genes", "x")
        for version in (1, 12)
        for name in figure_files(stem(token, desc, version, day=DAY))
    ]
    assert len(names) == len(MODULES) * 3 * 2 * 7
    assert all(check(n) for n in names)


@pytest.mark.parametrize(
    "name",
    [
        "2026-09-26_de_volcano_white_v1.pdf",
        "2026-09-26_de_volcano_v1_light.pdf",
        "2026-09-26_de_volcano_v1_dark.pdf",
        "2026-09-26_de_vol cano_v1.pdf",
        # A dataset label as it reads, never its name (D-114): a space, parentheses, a colon.
        "2026-09-26_de_Breast tissue (BRCA)_overview_v1.html",
        "2026-09-26_de_Breast-tissue-(BRCA)_overview_v1.html",
        "2026-09-26_de_BRCA-immune-cells:-macrophages_overview_v1.html",
        "2026-09-26_de_Breast-tissue-BRCA,_overview_v1.html",
        # A hyphen joins two runs of letters or digits, once.
        "2026-09-26_de_-Breast_v1.pdf",
        "2026-09-26_de_Breast-_v1.pdf",
        "2026-09-26_de_Breast--tissue_v1.pdf",
        "2026-09-26_de_Breast_-tissue_v1.pdf",
        "2026-09-26_de_Breast-tissue-BRCA_White_v1.pdf",
        "2026-09-26_de_Breast-tissue-BRCA_overview_v1_White.pdf",
        "2026-02-30_de_volcano_v1.pdf",
        "2026-09-26_de_volcano.pdf",
        "2026-09-26_de_volcano_vx.pdf",
        "2026-09-26_de_volcano_v0.pdf",
        "2026-09-26_nomodule_volcano_v1.pdf",
        "2026-09-26_de_volcano_v1_white_README.md",
        "2026-09-26_de_volcano_v1.pdf.bak",
    ],
)
def test_rejects(name):
    assert not check(name)


@pytest.mark.parametrize(
    "name",
    [
        "2026-10-23_de_Breast-tissue-BRCA_overview_v1.html",
        "2026-10-23_de_BRCA-immune-cells-macrophages_volcano_v1_white.svg",
        "2026-10-23_gather_Non-cancerous-blood-teaching-subset_samples_v1_README.md",
        "2026-10-23_qc_Lung-tissue-LUAD_reads_v2.csv",
        "2026-10-23_de_tcga_brca_volcano_v1.pdf",
    ],
)
def test_accepts_a_dataset_name_from_its_label(name):
    """Catches: the D-114 names refused, or the older lower-case slugs refused beside them."""
    assert check(name)


def test_description_reads_the_dataset_name_back():
    name = "2026-10-23_de_Pancreas-tissue-PAAD_overview_v1.html"
    assert description(name) == "Pancreas-tissue-PAAD_overview"


def test_stem_raises_on_illegal_parts():
    with pytest.raises(ValueError, match="legal"):
        stem("de", "Breast tissue (BRCA)", 1, day=DAY)


def test_lint_counts(tmp_path):
    good = ["2026-09-26_de_volcano_v1.csv", "2026-09-26_de_volcano_v1_README.md"]
    bad = ["volcano.csv", "sub/2026-09-26_de_volcano_v1_light.pdf"]
    for name in good + bad:
        (tmp_path / name).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / name).write_text("x")
    (tmp_path / ".DS_Store").write_text("x")
    scanned, violations = lint(tmp_path)
    assert scanned == len(good) + len(bad)
    assert sorted(violations) == sorted(bad)
    empty = tmp_path / "empty"
    empty.mkdir()
    assert lint(empty) == (0, [])


@pytest.mark.parametrize(
    ("table_stem", "columns"),
    [
        ("2026-09-26_de_table_v1", {**COLUMNS, "count": ("", "Reads.")}),
        ("2026-09-26_de_table_v1", {**COLUMNS, "score": ("log2 ratio", "")}),
        ("2026-09-26_de_table_v1", {k: v for k, v in COLUMNS.items() if k != "score"}),
        ("2026-09-26_de_Ta ble_v1", COLUMNS),
    ],
)
def test_companion_validates_before_writing(tmp_path, table_stem, columns):
    out = tmp_path / "out"
    with pytest.raises(ValueError):
        write_table(_frame(), out, table_stem, columns=columns, data_releases={})
    assert not out.exists()


def test_companion_describes_file_as_written(tmp_path):
    table_stem = stem("de", "table", 1, day=DAY)
    paths = write_table(
        _frame(), tmp_path, table_stem, columns=COLUMNS, data_releases={"GENCODE": "test"}
    )
    readme = (tmp_path / f"{table_stem}_README.md").read_text()
    assert [p.name for p in paths] == [
        f"{table_stem}.parquet",
        f"{table_stem}.csv",
        f"{table_stem}_README.md",
    ]
    assert all(check(p.name) for p in paths)
    for data in paths[:2]:
        assert hashlib.sha256(data.read_bytes()).hexdigest() in readme
    assert "| `count` | reads | Reads assigned to the gene. | `Int32` | `Int64` |" in readme
    assert "| `score` | log2 ratio | A planted Float32 value. | `Float32` | `Float64` |" in readme
    assert "| GENCODE | test |" in readme
    back = pl.read_csv(paths[1], schema_overrides=_frame().schema)
    assert back.equals(pl.read_parquet(paths[0]))


def test_readme_written_last(tmp_path, monkeypatch):
    def fail(*args, **kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(pl.DataFrame, "write_csv", fail)
    with pytest.raises(OSError, match="disk full"):
        write_table(
            _frame(), tmp_path, stem("de", "t", 1, day=DAY), columns=COLUMNS, data_releases={}
        )
    assert [p.name for p in tmp_path.iterdir()] == ["2026-09-26_de_t_v1.parquet"]


def _sha(repo: Path) -> str:
    return subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "HEAD"], capture_output=True, text=True, check=True
    ).stdout.strip()


def _unrelated_repo(path: Path) -> str:
    """A one-commit repository whose identity is its own: a CI runner has no global one, and
    `git commit` there exits 128 (run 37807649888)."""
    path.mkdir(parents=True)
    own = (["config", "user.name", "t"], ["config", "user.email", "t@t"],
           ["config", "commit.gpgsign", "false"])  # fmt: skip
    for args in (["init", "-q"], *own, ["commit", "-q", "--allow-empty", "-m", "unrelated"]):
        subprocess.run(["git", "-C", str(path), *args], check=True, capture_output=True)
    return _sha(path)


WRITE = """
    from datetime import date
    from pathlib import Path
    import polars as pl
    from ctx.io.readme_companion import write_table
    write_table(pl.DataFrame({"a": [1]}), Path("out"), "2026-09-26_de_t_v1",
                columns={"a": ("none", "x")}, data_releases={})
    print(Path("out/2026-09-26_de_t_v1_README.md").read_text())
"""


def test_provenance_not_from_cwd(tmp_path):
    other = _unrelated_repo(tmp_path / "other")
    readme = subprocess.run(
        [sys.executable, "-c", textwrap.dedent(WRITE)],
        cwd=tmp_path / "other",
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    assert other not in readme
    assert re.search(rf"\| ctx \| [^|]*git {_sha(ROOT)}", readme)

    host = tmp_path / "host"
    host_sha = _unrelated_repo(host)
    site = build_wheel(host)
    readme = run_from_wheel(site, WRITE, cwd=host)
    assert host_sha not in readme
    assert "installed wheel" in readme
