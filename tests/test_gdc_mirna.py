"""ctx.gather.gdc_mirna on real miRNA and isoform rows in made-up GDC folders
(tests/gdc_mirna_fixture.py), each test with the mutation it catches and was seen to fail."""

import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from gdc_fixture import aliquot
from gdc_fixture import build as build_star
from gdc_mirna_fixture import FILES, build, mirna_aliquot, rows, spec_file

from ctx.dataset import load
from ctx.gather import gdc, gdc_mirna, gdc_run

CTX = Path(sys.executable).parent / "ctx"
# One tumour and one normal of two patients: a1 and a2 are tumours, a3 and a4 normals.
STAR = (aliquot(1), aliquot(2), aliquot(1, code="11"), aliquot(2, code="11"))


def folders(tmp: Path, edit=None, files=FILES, star=STAR) -> tuple[gdc.Folder, gdc.Folder]:
    """Both kinds' folders of the made-up project, read and checked."""
    spec = load(spec_file(tmp / "x.toml")).gdc
    assert spec is not None
    small = [mirna_aliquot(a) for a in star]
    out = []
    for kind in (gdc_mirna.MIRNA, gdc_mirna.ISOFORM):
        build(tmp / kind, kind, small, files, edit=edit if kind == gdc_mirna.ISOFORM else None)
        out.append(gdc_mirna.read_kind(tmp / kind, spec, kind))
    return out[0], out[1]


def by_hand(name: str) -> tuple[dict[str, int], dict[str, int], dict[str, int]]:
    """From one fixture's isoform rows, typed loops: mature reads per MIMAT, reads per region,
    and reads per hairpin."""
    mature: dict[str, int] = defaultdict(int)
    region: dict[str, int] = defaultdict(int)
    hairpin: dict[str, int] = defaultdict(int)
    header, *body = rows(name, gdc_mirna.ISOFORM)
    at = {h: i for i, h in enumerate(header)}
    for r in body:
        reads = int(r[at["read_count"]])
        kind, _, mimat = r[at["miRNA_region"]].partition(",")
        region[kind] += reads
        hairpin[r[at["miRNA_ID"]]] += reads
        if kind == "mature":
            mature[mimat] += reads
    return mature, region, hairpin


def test_a_mature_mirna_sums_every_hairpin_and_only_mature_rows(tmp_path):
    """let-7a-5p (MIMAT0000062) comes from three hairpins. Catches: one hairpin taken for the
    sum (its count falls to about a third), and precursor or stem-loop rows counted as mature
    (miR-21's mature count rises by its precursor reads)."""
    found = gdc_mirna.read_matrices(*folders(tmp_path))
    column = {s: j for j, s in enumerate(found.files["sample"])}
    for k, name in enumerate(FILES):
        j = column[mirna_aliquot(STAR[k]).sample]
        mature, region, _ = by_hand(name)
        got = {
            m: c for m, c in zip(found.mature, found.mature_counts[:, j].tolist(), strict=True) if c
        }
        assert got == dict(mature)
        stats = found.stats.row(j, named=True)
        assert {r: stats[f"{r}_reads"] for r in region} == dict(region)
    let7a = found.mature.index("MIMAT0000062")
    hairpins = [h for h in found.hairpins if h.startswith("hsa-let-7a-")]
    assert len(hairpins) == 3
    assert (found.mature_counts[let7a] > found.hairpin_counts[found.hairpins.index(hairpins[0])]
            ).all()  # fmt: skip


def test_each_hairpins_isoform_rows_must_sum_to_its_mirna_count(tmp_path):
    """Catches: the check between the two kinds left out, so an isoform file with a row lost
    passes; the aliquot is named."""

    def drop(name: str, text: str) -> str:
        """The second isoform row of a2 removed."""
        if name != "a2":
            return text
        lines = text.splitlines(keepends=True)
        return "".join(lines[:2] + lines[3:])

    with pytest.raises(gdc.GdcError, match="do not sum") as error:
        gdc_mirna.read_matrices(*folders(tmp_path, edit=drop))
    assert mirna_aliquot(STAR[1]).sample in str(error.value)


def test_the_two_kinds_must_list_the_same_aliquots(tmp_path):
    """Catches: aliquots joined by position, or an isoform folder of another set of files."""
    spec = load(spec_file(tmp_path / "x.toml")).gdc
    assert spec is not None
    small = [mirna_aliquot(a) for a in STAR]
    build(tmp_path / "m", gdc_mirna.MIRNA, small)
    build(tmp_path / "i", gdc_mirna.ISOFORM, [*small[:3], mirna_aliquot(aliquot(3))])
    with pytest.raises(gdc.GdcError, match="different aliquots"):
        gdc_mirna.read_matrices(
            gdc_mirna.read_kind(tmp_path / "m", spec, gdc_mirna.MIRNA),
            gdc_mirna.read_kind(tmp_path / "i", spec, gdc_mirna.ISOFORM),
        )


def test_a_query_of_another_workflow_is_refused(tmp_path):
    """Catches: a folder of another miRNA workflow read as this one."""
    spec = load(spec_file(tmp_path / "x.toml")).gdc
    assert spec is not None
    build(tmp_path / "m", gdc_mirna.MIRNA, [mirna_aliquot(a) for a in STAR],
          workflow="Other miRNA Profiling")  # fmt: skip
    with pytest.raises(gdc.GdcError, match="workflow_type"):
        gdc_mirna.read_kind(tmp_path / "m", spec, gdc_mirna.MIRNA)


def test_hairpin_loci_span_every_file_and_name_their_matures(tmp_path):
    """Catches: a locus from one file alone (its span is narrower than the union's)."""
    found = gdc_mirna.read_matrices(*folders(tmp_path))
    spans = defaultdict(list)
    for name in FILES:
        header, *body = rows(name, gdc_mirna.ISOFORM)
        at = {h: i for i, h in enumerate(header)}
        for r in body:
            _, chrom, span, strand = r[at["isoform_coords"]].split(":")
            start, end = (int(v) for v in span.split("-"))
            spans[r[at["miRNA_ID"]]].append((chrom, start, end, strand))
    for row in found.loci.iter_rows(named=True):
        seen = spans[row["miRNA_ID"]]
        assert (row["chrom"], row["start"], row["end"], row["strand"]) == (
            seen[0][0],
            min(s[1] for s in seen),
            max(s[2] for s in seen),
            seen[0][3],
        )
    let7a = found.loci.filter(pl.col("miRNA_ID") == "hsa-let-7a-1")["mimats"][0]
    assert "MIMAT0000062" in let7a.split(";")


def test_a_mirna_file_joins_the_star_file_of_its_sample_barcode_not_its_aliquot(tmp_path):
    """The small-RNA and the mRNA library of one sample are different aliquots. Catches: the
    join on the aliquot (no file would join), and an FFPE miRNA library taken."""
    star = (
        aliquot(1),
        aliquot(2),
        aliquot(1, code="11"),
        aliquot(2, code="11", preservation="FFPE"),
    )
    build_star(tmp_path / "star", star)
    spec = load(spec_file(tmp_path / "x.toml"))
    read = gdc_run.read_all(tmp_path / "star", spec)
    assert spec.gdc is not None
    small = [mirna_aliquot(a) for a in star]
    for kind in (gdc_mirna.MIRNA, gdc_mirna.ISOFORM):
        build(tmp_path / kind, kind, small)
    found = gdc_run.read_mirna(tmp_path / gdc_mirna.MIRNA, tmp_path / gdc_mirna.ISOFORM,
                               spec.gdc, read)  # fmt: skip
    taken = gdc.chosen(found.sheet)
    assert sorted(taken["sample"]) == sorted(s.sample for s in small[:3])
    assert dict(taken.select("sample", "rna_sample").rows()) == {
        small[i].sample: star[i].sample for i in range(3)
    }
    ffpe = found.sheet.filter(pl.col("sample") == small[3].sample)
    assert ffpe["drop_reason"][0] == "ffpe"


def test_gather_writes_the_mirna_tables_beside_the_star_ones(tmp_path):
    """The CLI path end to end. Catches: a count column out of the miRNA sheet's order, a table
    written without its README, and --mirna-dir accepted without --isoform-dir."""
    build_star(tmp_path / "star", STAR)
    small = [mirna_aliquot(a) for a in STAR]
    for kind in (gdc_mirna.MIRNA, gdc_mirna.ISOFORM):
        build(tmp_path / kind, kind, small)
    spec = spec_file(tmp_path / "fake_brca.toml")
    base = [str(CTX), "--root", str(tmp_path / "root"), "gather", "--dataset", "fake_brca",
            "--gdc-dir", str(tmp_path / "star"), "--dataset-file", str(spec)]  # fmt: skip
    half = subprocess.run([*base, "--mirna-dir", str(tmp_path / gdc_mirna.MIRNA)],
                          capture_output=True, text=True, check=False)  # fmt: skip
    assert half.returncode == 1 and "together" in half.stderr
    done = subprocess.run([*base, "--mirna-dir", str(tmp_path / gdc_mirna.MIRNA),
                           "--isoform-dir", str(tmp_path / gdc_mirna.ISOFORM)],
                          capture_output=True, text=True, check=False)  # fmt: skip
    assert done.returncode == 0, done.stderr
    out = tmp_path / "root" / "outputs" / "01_gather"
    sheet = pl.read_parquet(next(out.glob("*_fake_brca_mirna_samplesheet_v1.parquet")))
    for what, key in (("mirna_mature_counts", "mimat"), ("mirna_hairpin_counts", "hairpin")):
        table = pl.read_parquet(next(out.glob(f"*_fake_brca_{what}_v1.parquet")))
        assert table.columns == [key, *sheet["sample"]]
    hairpins = pl.read_parquet(next(out.glob("*_fake_brca_mirna_hairpin_counts_v1.parquet")))
    first = {r[0]: int(r[1]) for r in rows("a1", gdc_mirna.MIRNA)[1:]}
    column = sheet.filter(pl.col("sample") == mirna_aliquot(STAR[0]).sample)["sample"][0]
    assert dict(hairpins.select("hairpin", column).rows()) == first
    for table in out.glob("*.parquet"):
        assert table.with_name(table.name.replace(".parquet", "_README.md")).is_file()
    assert np.isclose(
        sheet["mirna_reads"].sum(),
        sum(sum(int(r[1]) for r in rows(n, gdc_mirna.MIRNA)[1:]) for n in FILES),
    )


# Measured 2026-10-08 by this reader on the folders on this Mac, and agreeing with the plan's
# prototype where it measured the same thing (docs/plans/2026-10-07_m8_real_data_plan_proposal.md,
# section 1): aliquots, files taken, MIMATs, hairpins with a read, miRNA-mapped reads,
# tumour-normal pairs taken, and pairs with both layers on one sample barcode.
FOUR = {
    "tcga_brca": ("ctx_breast", 1207, 1188, 2236, 1604, 4_217_403_486, 103, 103),
    "tcga_kirc": ("ctx_kidney", 616, 588, 2102, 1502, 1_852_722_023, 71, 71),
    "tcga_luad": ("ctx_lung", 567, 561, 2213, 1586, 2_730_063_500, 46, 19),
    "tcga_lusc": ("ctx_lung", 533, 533, 2241, 1603, 2_023_836_222, 45, 38),
}


def _paired(frame: pl.DataFrame) -> set[str]:
    """Patients with a tumour and a normal library in the frame."""
    return (set(frame.filter(pl.col("condition") == "tumour")["case"])
            & set(frame.filter(pl.col("condition") == "normal")["case"]))  # fmt: skip


@pytest.mark.realdata
@pytest.mark.parametrize("slug", sorted(FOUR))
def test_the_four_tcga_projects_give_the_measured_numbers(slug):
    """Opt-in (-m realdata): each project's three folders read as gather reads them. Catches
    rows dropped or doubled by the reader, the rule or the join, at real size."""
    home, *expected = FOUR[slug]
    raw = Path(os.environ.get("CTX_GDC_HOME", Path.home())) / home / "data/raw"
    if not all((raw / f"{slug}{k}").is_dir() for k in ("", "_mirna", "_isoform")):
        pytest.skip(f"{raw}/{slug} folders are not on this machine")
    spec = load(Path(__file__).resolve().parents[1] / "conf" / "datasets" / f"{slug}.toml")
    assert spec.gdc is not None
    star = gdc_run.read_all(raw / slug, spec)
    found = gdc_run.read_mirna(raw / f"{slug}_mirna", raw / f"{slug}_isoform", spec.gdc, star)
    taken = gdc.chosen(found.sheet)
    m = found.matrices
    measured = (found.sheet.height, taken.height, len(m.mature),
                int((m.hairpin_counts.sum(axis=1) > 0).sum()),
                int(found.sheet["mirna_reads"].sum()), len(_paired(taken)),
                len(_paired(taken.filter(pl.col("rna_sample").is_not_null()))))  # fmt: skip
    assert measured == tuple(expected)
