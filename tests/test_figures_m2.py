"""F.1: M2's figures, each redrawn from its spec's inputs alone in a clean folder.

Every figure with a spec writes M2's data table, PDF and SVG byte for byte. Each test names
its mutation.
"""

import re
from datetime import date
from pathlib import Path

import polars as pl
import pytest
from conftest import DE_SPEC, once
from test_figures import copied, drawn_alike, redrawn, redrawn_each, specs
from test_qc_run import with_reads, with_reports  # noqa: F401 - session fixtures used below

from ctx.figures.spec import parse

QC = ("02_qc",)
# The linked PCA prints its variance shares on the page's axes alone, not in its PDF or SVG.
SHARES = re.compile(r"PC[12], \d+% of variance")


@pytest.fixture(scope="session")
def with_fastqc(tmp_path_factory):
    """ctx qc --rnaseq-dir on a copy of the rnaseq fixture holding FastQC's raw-read reports."""
    return once(tmp_path_factory, "qc_fastqc", _with_fastqc)


def _with_fastqc(tmp: Path) -> Path:
    from synth import Synthetic, design_counts, write_gather_outputs  # noqa: PLC0415
    from test_qc_fates import with_raw_fastqc  # noqa: PLC0415
    from test_qc_run import RNASEQ_SAMPLES, ctx  # noqa: PLC0415

    rnaseq, root = with_raw_fastqc(tmp), tmp / "root"
    synthetic = design_counts((("control", "short", 3), ("case", "short", 2)))
    names = dict(zip(synthetic.samples["sample"], RNASEQ_SAMPLES, strict=True))
    renamed = Synthetic(synthetic.counts.rename(names),
                        synthetic.samples.with_columns(pl.col("sample").replace_strict(names)),
                        synthetic.truth)  # fmt: skip
    write_gather_outputs(root, "desynth", renamed, date(2026, 9, 30))
    done = ctx("--root", str(root), "qc", "--dataset", "desynth", "--dataset-file", str(DE_SPEC),
               "--rnaseq-dir", str(rnaseq), cwd=root)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    return root


def _kinds(found) -> set[str]:
    """The kinds, a strip and a FastQC plot named by their parameter."""
    names = set()
    for spec in found:
        kind = parse(spec.figure_id)[2]
        extra = getattr(spec.params, "metric", None) or getattr(spec.params, "plot", None)
        names.add(f"{kind}:{extra}" if extra else kind)
    return names


def _redraws_alike(root: Path, found, tmp_path: Path) -> None:
    """Each spec redrawn in its own clean folder equals M2's files."""
    for spec, csv in zip(found, redrawn_each(found, root, tmp_path), strict=True):
        module_csv = next(root.glob(f"outputs/02_qc/{spec.stem}.csv"))
        assert drawn_alike(csv, module_csv) == [], spec.stem
        shown = (p.with_suffix(".html").read_text() for p in (module_csv, csv))
        assert len({tuple(SHARES.findall(page)) for page in shown}) == 1, spec.stem


def test_every_count_figure_redraws_byte_for_byte(qc_run, tmp_path):
    """The counts' figures, band map, PCA, size against genes, sex check and every strip,
    redrawn from M2's tables. Catches the PCA's shares read in the wrong order, the MAD
    multiplier taken from the registry file instead of the run, a sheet without its read-length
    groups, and a kind that no longer matches its stem."""
    root = qc_run[0]
    found = specs(root, QC, DE_SPEC)
    names = _kinds(found)
    assert {"linked_pca", "size_genes", "sex_check"} <= names, names
    assert {"strip:library_size", "strip:genes_detected", "strip:pca_distance"} <= names, names
    _redraws_alike(root, found, tmp_path)


def test_report_figures_redraw_byte_for_byte(with_reports, tmp_path):  # noqa: F811 - fixture
    """The pipeline reports' figures, read fates and flow, redrawn from M2's fates table.
    Catches the fates read with the values of another table."""
    root = with_reports[0]
    found = specs(root, QC, DE_SPEC)
    assert {"read_fates", "read_flow"} <= _kinds(found)
    _redraws_alike(root, found, tmp_path)


def test_fastqc_plots_redraw_from_their_tables(with_fastqc, tmp_path):
    """FastQC's five plots, redrawn from M2's copy of MultiQC's plot data, byte for byte, the
    GC plot over its theoretical curve. Catches a curve read back under another metric and the
    theoretical GC curve dropped."""
    found = specs(with_fastqc, QC, DE_SPEC)
    plots = {"base_quality", "fastqc_adapter", "fastqc_n_content", "gc_distribution",
             "overrepresented"}  # fmt: skip
    assert {f"fastqc:{n}" for n in plots} <= _kinds(found), _kinds(found)
    theoretical = next(with_fastqc.glob("outputs/02_qc/*_fastqc_theoretical_v1.parquet"))
    assert pl.read_parquet(theoretical).height > 0
    _redraws_alike(with_fastqc, [s for s in found if parse(s.figure_id)[2] == "fastqc"],
                   tmp_path)  # fmt: skip


def test_reads_figures_redraw_and_shared_stems_have_no_spec(with_reads, tmp_path):  # noqa: F811
    """With --reads, the grid and the runs' band map redraw byte for byte; no strip and no
    FastQC plot carries a spec, since the reads' figures take the samples' stems. Catches the
    unless rule left out, which would attach the samples' spec to the reads' strip."""
    out = with_reads[0]
    root = out.parent.parent
    found = specs(root, QC, DE_SPEC)
    names = _kinds(found)
    assert {"fastp_grid", "run_band_map"} <= names, names
    assert not [n for n in names if n.startswith(("strip:", "fastqc:"))], names
    _redraws_alike(root, found, tmp_path)


def test_the_registry_named_by_the_run_is_the_one_copied(qc_run, tmp_path):
    """A strip's panel copies the registry the run read, as its options table names it, and a
    multiplier moved there moves the redrawn lead. Catches k taken from the registry file
    instead of the run's options, which --mad overrides."""
    root = qc_run[0]
    spec = next(s for s in specs(root, QC, DE_SPEC)
                if getattr(s.params, "metric", "") == "library_size")  # fmt: skip
    folder = copied(spec, root, tmp_path / "k")
    options = folder / next(i.file for i in spec.inputs if i.name == "figure_options")
    row = pl.read_parquet(options)
    used = row["registry"][0]
    registry = next(i for i in spec.inputs if i.name == "registry")
    assert Path(used).name == registry.file
    row.with_columns(mad=pl.lit(7.0)).write_parquet(options)
    raw = (folder / "params.toml").read_text().splitlines()
    (folder / "params.toml").write_text(
        "\n".join(x for x in raw if not x.startswith("figure_options = ") or ".parquet" in x) + "\n"
    )
    page = redrawn(folder, spec).with_suffix(".html").read_text()
    assert "median ± 7 MAD" in page
    assert f"median ± {row['mad'][0]:g} MAD" not in page


def test_a_run_with_its_own_registry_names_that_registry(tmp_path):
    """ctx qc --registry with a copy of another name and another k: every spec copies that
    copy, and the strips redraw byte for byte. Catches the registry resolved afresh, the
    course's file, whose k draws other bands."""
    from test_qc_run import ROOT, ctx, plant  # noqa: PLC0415

    root = tmp_path / "root"
    plant(root)
    own = tmp_path / "own_registry.toml"
    raw = (ROOT / "conf" / "qc_metrics.toml").read_text()
    assert raw.count("mad_multiplier = 3.0") == 1
    own.write_text(raw.replace("mad_multiplier = 3.0", "mad_multiplier = 2.0"))
    done = ctx("--root", str(root), "qc", "--dataset", "desynth", "--dataset-file", str(DE_SPEC),
               "--registry", str(own), cwd=root)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    found = specs(root, QC, DE_SPEC)
    named = {i.file for s in found for i in s.inputs if i.name == "registry"}
    assert named == {own.name}, named
    _redraws_alike(root, [s for s in found if parse(s.figure_id)[2] == "strip"][:3], tmp_path)
