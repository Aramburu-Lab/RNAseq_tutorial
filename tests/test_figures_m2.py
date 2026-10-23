"""F.1: M2's figures, each redrawn from its spec's inputs alone in a clean folder.

Every figure with a spec writes M2's data table, PDF and SVG byte for byte. Each test names
its mutation.
"""

import re
from pathlib import Path

import polars as pl
from conftest import DE_SPEC
from test_figures import copied, drawn_alike, redrawn, redrawn_each, specs
from test_qc_run import with_reports  # noqa: F401 - session fixtures used below

from ctx.figures.spec import parse
from ctx.report.builder import unpacked

QC = ("02_qc",)
# The linked PCA prints its variance shares on the page's axes alone, not in its PDF or SVG.
SHARES = re.compile(r"PC[12], \d+% of variance")


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
        shown = (unpacked(p.with_suffix(".html").read_text()) for p in (module_csv, csv))
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


def test_fastqc_plots_redraw_from_their_tables(with_reports, tmp_path):  # noqa: F811 - fixture
    """FastQC's five plots, redrawn from M2's copy of nf-core's MultiQC plot data, byte for byte;
    nf-core's MultiQC draws no theoretical GC curve, so that table is empty (D-116). Catches a
    curve read back under another metric."""
    root = with_reports[0]
    found = specs(root, QC, DE_SPEC)
    plots = {"base_quality", "fastqc_adapter", "fastqc_n_content", "gc_distribution",
             "overrepresented"}  # fmt: skip
    assert {f"fastqc:{n}" for n in plots} <= _kinds(found), _kinds(found)
    theoretical = next(root.glob("outputs/02_qc/*_fastqc_theoretical_v1.parquet"))
    assert pl.read_parquet(theoretical).height == 0
    _redraws_alike(root, [s for s in found if parse(s.figure_id)[2] == "fastqc"], tmp_path)


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
    page = unpacked(redrawn(folder, spec).with_suffix(".html").read_text())
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
