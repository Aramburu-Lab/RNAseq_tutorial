"""S2 (R2.5): Ensembl gene IDs the annotation lacks, placed through Ensembl's ID history.

The annotation is test_gather_counts's nine made-up genes; the history is a made-up
``stable_id_event`` dump in the column order of Ensembl's core schema. Each test names the
mutation it catches.
"""

from pathlib import Path

import polars as pl
import pytest
from conftest import ctx_run
from test_gather_counts import GENES, matrix, raw, root_with, whole

from ctx.gather import counts
from ctx.report.builder import unpacked

NULL = "\\N"


def gene(n: int) -> str:
    """A made-up stable gene ID; 1 to 9 are the annotation's genes."""
    return f"ENSG{n:011d}"


# (old, new, type): Ensembl's stable_id_event rows, version 1 to 2, session 1, score 1.
EVENTS = [
    (gene(91), gene(3), "gene"),  # retired, one current successor
    (gene(92), gene(99), "gene"),  # through an ID itself retired, to a current one
    (gene(99), gene(4), "gene"),
    (gene(93), gene(7), "gene"),  # split into two current genes
    (gene(93), gene(8), "gene"),
    (gene(94), NULL, "gene"),  # retired with no successor
    (gene(95), gene(1), "gene"),  # onto a gene the matrix holds under its own ID
    (gene(96), gene(96), "gene"),  # a version step, no new ID
    (gene(97), gene(9), "transcript"),  # not a gene event
    (gene(98), gene(2), "gene"),  # onto a current gene that later split
    (gene(2), gene(6), "gene"),
    (gene(81), gene(5), "gene"),  # two retired IDs onto one gene
    (gene(82), gene(5), "gene"),
    (gene(83), gene(84), "gene"),  # a cycle, then a current gene
    (gene(84), gene(83), "gene"),
    (gene(84), gene(9), "gene"),
    (NULL, gene(9), "gene"),  # a gene created
]
IDS = [f"{GENES[0][0]}", *(f"{gene(n)}.3" for n in (91, 92, 93, 94, 95, 96, 97, 98, 81, 82, 83))]


def history(path: Path, events=EVENTS, columns: int = 7) -> Path:
    """The dump as Ensembl's FTP holds it: tab-separated, no header, NULL as backslash N."""
    fields = [[old, "1", new, "2", "1", kind, "1", "extra"][:columns] for old, new, kind in events]
    lines = ["\t".join(f) for f in fields]
    path.write_text("\n".join(lines) + "\n")
    return path


def annotation() -> pl.DataFrame:
    """The made-up annotation's genes and names."""
    return pl.DataFrame(
        {"gene_id": [g for g, _, _ in GENES], "gene_name": [n for _, n, _ in GENES]}
    )


def test_history_places_only_one_current_gene(tmp_path):
    """A retired ID is placed where its history, followed through every later event, reaches
    exactly one gene of the annotation, and never unseats a row placed by its own ID. Catches
    the walk stopped at the first event (92) or at the first current gene (98), a split given
    to one side (93), version steps or transcript events read as successors (96, 97), a cycle
    that never ends (83), a fallback row kept over the matrix's own (95), two fallback rows
    on one gene where none is the matrix's own both kept (81, 82), and the count of rows the
    annotation lacks."""
    found = counts.id_history(history(tmp_path / "h.txt"))
    assert found[gene(93)] == {gene(7), gene(8)} and gene(96) not in found
    assert gene(97) not in found and gene(94) not in found
    mapped = counts.place(raw(IDS, whole(len(IDS))), annotation(), None, None, found)
    assert mapped.genes.select("gene_name", "source_id").rows() == [
        ("ALPHA", IDS[0]), ("GAMMA", f"{gene(91)}.3"), ("DELTA", f"{gene(92)}.3"),
        ("IOTA", f"{gene(83)}.3"),
    ]  # fmt: skip
    left = {r[0].split(".")[0]: (r[2], r[3]) for r in mapped.left_out.rows()}
    assert left == {
        gene(93): ("not_in_annotation", f"ID history: 2 current genes, {gene(7)}, {gene(8)}"),
        gene(94): ("not_in_annotation", "ID history: 0 current genes"),
        gene(95): ("shared_gene", f"{GENES[0][0]}, with {IDS[0]}, {gene(95)}.3"),
        gene(96): ("not_in_annotation", "ID history: 0 current genes"),
        gene(97): ("not_in_annotation", "ID history: 0 current genes"),
        gene(98): ("not_in_annotation", f"ID history: 2 current genes, {gene(2)}, {gene(6)}"),
        gene(81): ("shared_gene", f"{GENES[4][0]}, with {gene(81)}.3, {gene(82)}.3"),
        gene(82): ("shared_gene", f"{GENES[4][0]}, with {gene(81)}.3, {gene(82)}.3"),
    }
    assert (mapped.through_history, mapped.lacked) == (3, 11)
    without = counts.place(raw(IDS, whole(len(IDS))), annotation(), None)
    assert without.genes["source_id"].to_list() == [IDS[0]]
    assert (without.through_history, without.lacked) == (0, 11)
    assert set(without.left_out["detail"].to_list()) == {""}


def test_history_with_another_shape_is_refused(tmp_path):
    """A file of six or eight columns is named, never read by guessed positions. Catches a
    dump of another table, or another schema, read as the ID history."""
    for columns in (6, 8):
        with pytest.raises(counts.CountsError, match="stable_id_event"):
            counts.id_history(history(tmp_path / "h.txt", columns=columns))


def test_gather_reads_the_history_resource(tmp_path):
    """With ensembl_id_history declared, gather places a retired ID through it and the page
    and README say so. Catches the resource declared and never read, and the page's count."""
    root, spec = root_with(tmp_path)
    history(root / "resources" / "h.txt")
    resource = ('\n[[resources]]\nname = "ensembl_id_history"\nurls = ["http://127.0.0.1/h"]\n'
                'bytes = 0\nsha256 = "' + "0" * 64 + '"\ndest = "resources/h.txt"\n')  # fmt: skip
    spec.write_text(spec.read_text().replace("[counts]", resource + "\n[counts]"))
    ids = [IDS[0], f"{gene(91)}.3", f"{gene(94)}.3"]
    path = matrix(tmp_path / "m.tsv", ids, whole(3))
    done = ctx_run("--root", str(root), "gather", "--dataset", "mat", "--counts", str(path),
                   "--dataset-file", str(spec), cwd=tmp_path)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    out = root / "outputs" / "01_gather"
    (genes,) = out.glob("*_gather_mat_genes_v1.parquet")
    assert pl.read_parquet(genes)["source_id"].to_list() == ids[:2]
    page = unpacked(next(out.glob("*_gather_mat_overview_v1.html")).read_text())
    assert "1 of 2 Ensembl rows whose ID the annotation lacks" in page
    # R2.5b: the largest rows left out show why, here a retired ID with no successor.
    assert "<td>ID history: 0 current genes</td>" in page
    readme = next(out.glob("*_counts_v1_README.md")).read_text()
    assert "Ensembl stable_id_event, h.txt" in readme
