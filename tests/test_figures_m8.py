"""F.1: M8's figures, each redrawn from its spec's inputs alone in a clean folder.

Every figure with a spec writes M8's data table, PDF and SVG byte for byte, and its page
carries the same live controls, where the pairs not drawn first have their values alone.
F5, the library types, is drawn by no fixture: no made-up sheet carries GDC's read groups.
Each test names its mutation.
"""

import re
from datetime import date
from pathlib import Path

import polars as pl
import pytest
from conftest import ctx_run, once
from test_figures import copied, drawn_alike, redrawn, redrawn_each, specs

from ctx.figures.spec import parse

NCRNA = ("08_ncrna",)
STEP1 = {"biotype_shares", "filter_by_biotype", "volcano_by_biotype", "lncrna_neighbour"}
MIRNA = {"mirna_volcano", "target_shift", "mirna_target_pairs"}
CONTROLS = re.compile(r'<script type="application/json" id="ctx-controls-\d+">(.*?)</script>',
                      flags=re.S)  # fmt: skip


def _run(root: Path, spec: Path, dataset: str, *options: str) -> Path:
    """ctx ncrna on a root; its output folder."""
    done = ctx_run("--root", str(root), "ncrna", "--dataset", dataset, "--dataset-file",
                   str(spec), *options, cwd=root)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    return root / "outputs" / "08_ncrna"


@pytest.fixture(scope="session")
def tcga(ncrna_tcga_run) -> tuple[Path, Path]:
    """ctx ncrna on the made-up GDC root: the root and its dataset file."""
    return ncrna_tcga_run[:2]


@pytest.fixture(scope="session")
def optioned(tmp_path_factory) -> tuple[Path, Path]:
    """ctx ncrna on the made-up Step-1 root with three biotype groups and another count."""
    return once(tmp_path_factory, "ncrna_options", _optioned)


def _optioned(folder: Path) -> tuple[Path, Path]:
    from ncrna_root import DATASET, make  # noqa: PLC0415

    root, spec = make(folder, date.today())
    _run(root, spec, DATASET, "--biotypes", "protein_coding,lncRNA,pseudogene",
         "--min-count", "5")  # fmt: skip
    return root, spec


def _redraws_alike(root: Path, found, tmp_path: Path) -> set[str]:
    """Each spec redrawn in its own clean folder equals M8's files and live controls.

    Returns:
        The kinds whose page carries a live control.
    """
    live = set()
    for spec, csv in zip(found, redrawn_each(found, root, tmp_path), strict=True):
        module_csv = next(root.glob(f"outputs/08_ncrna/{spec.stem}.csv"))
        assert drawn_alike(csv, module_csv) == [], spec.stem
        shown = [CONTROLS.findall(p.with_suffix(".html").read_text()) for p in (module_csv, csv)]
        assert shown[0] == shown[1], spec.stem
        live |= {parse(spec.figure_id)[2]} if shown[0] else set()
    return live


def _kinds(found) -> set[str]:
    """The kinds of the specs found."""
    return {parse(s.figure_id)[2] for s in found}


def test_step1_figures_redraw_byte_for_byte(ncrna_run, tmp_path):
    """The four Step-1 figures, redrawn from gather's sheet and counts, M5's results and M8's
    gene records and options. Catches the gene records and options written after the first
    figure, which leaves it without a spec, and gene records out of the counts' row order."""
    root = ncrna_run[0]
    found = specs(root, NCRNA, root.parent / "ncsynth.toml")
    assert _kinds(found) == STEP1
    assert {"volcano_by_biotype", "lncrna_neighbour"} <= _redraws_alike(root, found, tmp_path)


def test_the_run_s_options_are_the_ones_redrawn(optioned, tmp_path):
    """A run with three biotype groups and a minimum count of 5: the filter curves and the
    volcano redraw byte for byte. Catches the drawer drawing the default groups or count
    instead of the options M8 wrote, and the groups read back in another order."""
    root, spec = optioned
    options = pl.read_parquet(next(root.glob("outputs/08_ncrna/*_figure_options_v1.parquet")))
    assert options.select("biotypes", "min_count").row(0) == ("protein_coding,lncRNA,pseudogene", 5)
    found = [s for s in specs(root, NCRNA, spec)
             if parse(s.figure_id)[2] in {"filter_by_biotype", "volcano_by_biotype"}]  # fmt: skip
    assert len(found) == 2
    _redraws_alike(root, found, tmp_path)


def test_gdc_figures_redraw_byte_for_byte(tcga, tmp_path):
    """The Step-1 figures and the miRNA half on a GDC root: the miRNA volcano from M8's fit
    unrefitted, the target shifts from the four target files, the pairs from gather's miRNA
    tables. Catches the target files left unread, and the target files' inputs left out of
    the target shift's spec."""
    root, spec = tcga
    found = specs(root, NCRNA, spec)
    assert _kinds(found) == STEP1 | MIRNA
    assert {"mirna_volcano", "mirna_target_pairs"} <= _redraws_alike(root, found, tmp_path)


def test_a_drawer_takes_the_sheet_through_the_gdc_door(tcga, tmp_path):
    """The redrawn biotype bars name the twelve patients' chosen libraries and never the
    unchosen extra aliquot of the fixture. Catches the drawer reading gather's sheet without
    gdc.chosen, which would draw S999's bar beside the others."""
    import ncrna_tcga_root  # noqa: PLC0415

    root, spec = tcga
    one = next(s for s in specs(root, NCRNA, spec) if parse(s.figure_id)[2] == "biotype_shares")
    folder = copied(one, root, tmp_path / "door")
    sheet = pl.read_parquet(folder / next(i.file for i in one.inputs if i.name == "samplesheet"))
    assert ncrna_tcga_root.EXTRA in set(sheet["sample"])
    drawn = pl.read_csv(redrawn(folder, one))
    assert ncrna_tcga_root.EXTRA not in set(drawn["sample"])
    assert drawn["sample"].n_unique() == 2 * ncrna_tcga_root.PAIRS


def test_the_mirna_sheet_passes_the_gdc_door(tcga, tmp_path):
    """The miRNA volcano and pairs, redrawn from a miRNA sheet that also lists an unchosen
    patient's two libraries, equal M8's. Catches either drawer reading the miRNA sheet without
    gdc.chosen: the volcano would count 13 patients a group, the pairs ask for a library the
    mature counts lack."""
    root, spec = tcga
    found = [s for s in specs(root, NCRNA, spec)
             if parse(s.figure_id)[2] in {"mirna_volcano", "mirna_target_pairs"}]  # fmt: skip
    assert len(found) == 2
    for one in found:
        folder = copied(one, root, tmp_path / one.stem)
        path = folder / next(i.file for i in one.inputs if i.name == "mirna_samplesheet")
        sheet = pl.read_parquet(path)
        extra = sheet.unique("condition", keep="first", maintain_order=True).with_columns(
            case=pl.lit("P99"), donor=pl.lit("P99"), sample_id="P99" + pl.col("condition"),
            sample="MP99" + pl.col("condition"), rna_sample=pl.lit("S999"),
            run_accession=pl.lit("AAA"), chosen=pl.lit(value=False),
            drop_reason=pl.lit("not_first"),
        )  # fmt: skip
        pl.concat([sheet, extra]).write_parquet(path)
        params = (folder / "params.toml").read_text().splitlines()
        kept = [x for x in params if not (x.startswith("mirna_samplesheet = ") and len(x) > 80)]
        (folder / "params.toml").write_text("\n".join(kept) + "\n")
        module_csv = next(root.glob(f"outputs/08_ncrna/{one.stem}.csv"))
        assert drawn_alike(redrawn(folder, one), module_csv) == [], one.stem


def test_the_gene_records_are_the_counts_rows(ncrna_run):
    """M8's gene records hold one row per count row, in count order, each with its group.
    Catches the records written in the annotation's order, which the pairs index by row."""
    root = ncrna_run[0]
    out = root / "outputs"
    genes = pl.read_parquet(next(out.glob("08_ncrna/*_gene_biotypes_v1.parquet")))
    counts = pl.read_parquet(next(out.glob("01_gather/*_counts_v1.parquet")))
    assert genes["gene_id"].to_list() == counts["gene_id"].to_list()
    assert genes["biotype"].null_count() == 0


def test_a_moved_target_file_is_refused(tcga, tmp_path):
    """A target file changed in the clean folder stops the redraw. Catches the target files
    left out of the spec's inputs, which a redraw would then find nowhere."""
    root, spec = tcga
    one = next(s for s in specs(root, NCRNA, spec) if parse(s.figure_id)[2] == "target_shift")
    names = {i.name for i in one.inputs}
    assert {"targetscan_predictions", "targetscan_families", "targetscan_genes",
            "mirtarbase_strong"} <= names  # fmt: skip
    folder = copied(one, root, tmp_path / "moved")
    target = folder / next(i.file for i in one.inputs if i.name == "mirtarbase_strong")
    target.write_text(target.read_text() + "\n")
    done = ctx_run("figure", one.figure_id, "--params", "params.toml", cwd=folder)
    assert done.returncode == 1 and "inputs.mirtarbase_strong" in done.stderr, done.stderr


def test_library_types_reads_every_file_of_the_sheet(tcga, tmp_path):
    """F5's drawer on tests/test_ncrna_tcga.py's two tumours sequenced both ways, its tables
    as a GDC root holds them: the same pairs as M8's own block on the whole sheet. Catches the
    drawer handing F5 the chosen files alone, where no tumour has its rRNA-depleted twin."""
    from test_ncrna_tcga import F5_COUNTS, F5_GENES, TUMOUR, f5_files  # noqa: PLC0415

    from ctx.dataset import load  # noqa: PLC0415
    from ctx.figures import draw_m8  # noqa: PLC0415
    from ctx.figures.kind import Plain  # noqa: PLC0415
    from ctx.figures.recorder import Recorder  # noqa: PLC0415
    from ctx.ncrna import libraries_page, run  # noqa: PLC0415

    spec = load(tcga[1])
    sheet = f5_files().with_columns(
        chosen=pl.col("drop_reason") == "",
        run_accession=pl.col("sample"),
        donor=pl.col("sample_id"),
    )
    tables = {"samplesheet": sheet, "counts": F5_COUNTS, "gene_biotypes": F5_GENES}
    inputs = {}
    for name, frame in tables.items():
        frame.write_parquet(inputs.setdefault(name, tmp_path / f"{name}.parquet"))
    data, _ = draw_m8.library_types(inputs, Plain(spec.slug), spec)
    recorder = Recorder(draw_m8.FOLDER, spec.slug)
    assert spec.ncrna is not None and spec.ncrna.mirna is not None
    whole = run.Tcga(spec.ncrna.mirna, sheet, sheet, *(pl.DataFrame(),) * 3, F5_COUNTS,
                     pl.DataFrame(), draw_m8.NO_TARGETS)  # fmt: skip
    libraries_page.group(recorder, whole, F5_GENES, TUMOUR)
    assert data.frame.equals(recorder.one("library_types")[0].frame)
    assert sorted(data.frame["sample_id"]) == ["A", "B"]
