"""F.1: M7's figures, each redrawn from its spec's inputs alone in a clean folder.

Every design's figures with a spec write M7's data table, PDF and SVG byte for byte, and the
words their page prints around them. Each test names its mutation.
"""

import re
import shutil
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from clinical_root import ctx
from conftest import proteo_study, writable_copy
from proteo_study_root import SLUG
from test_figures import copied, drawn_alike, redrawn, redrawn_each, specs

from ctx.dataset import load
from ctx.figures.kinds_m7 import TABLES
from ctx.figures.params import ParamsError, load_params
from ctx.figures.recorder import Recorder
from ctx.figures.spec import parse
from ctx.proteo import figure_tables, figures
from ctx.proteo.page import PAIRED_TOP
from ctx.proteo.run import Measured
from ctx.proteo.tables import Layer
from ctx.report.builder import unpacked

M7 = ("07_proteo",)
MISSING = {"quantified", "coverage"}
CORRELATION = {*MISSING, "rho_by_class", "worked:G9", "paired_zscores"}
# What each made-up root draws: its design's figures, F4 for each contrast with M5's results.
EXPECTED = {"across_root": {*CORRELATION, "contrast:basal_vs_luma"},
            "cross_root": {*MISSING, "fold_change", "shift_by_class"},
            "study_root": CORRELATION,
            "made": {*CORRELATION, "rho_by_class_normal", "contrast:tumour_normal"}}  # fmt: skip
# A figure page's words a PDF does not carry: its lead and its reading guide.
WORDS = re.compile(r'<p class="fig-lead">.*?</p>|<aside class="reading".*?</aside>', re.S)


@pytest.fixture(scope="session")
def across_root(proteo_across) -> tuple[Path, Path]:
    """The across-samples root, ctx proteo run, and its dataset file."""
    return proteo_across


@pytest.fixture(scope="session")
def cross_root(proteo_cross) -> tuple[Path, Path]:
    """The cross-study root, ctx proteo run, and its dataset file."""
    return proteo_cross


@pytest.fixture(scope="session")
def study_root(tmp_path_factory) -> tuple[Path, Path]:
    """A [proteo.table] study root, ctx proteo run, and its dataset file."""
    root = proteo_study(tmp_path_factory, "maxquant")
    return root, root.parent / f"{SLUG}.toml"


@pytest.fixture(scope="session")
def made(proteo_made) -> tuple[Path, Path]:
    """The tumour-against-normal root, ctx proteo run, and its dataset file."""
    return proteo_made


def _kinds(found) -> set[str]:
    """The kinds, a worked gene and a contrast named by their parameter."""
    names = set()
    for spec in found:
        kind = parse(spec.figure_id)[2]
        extra = getattr(spec.params, "gene", None) or getattr(spec.params, "contrast", None)
        names.add(f"{kind}:{extra}" if extra else kind)
    return names


def drawn(root: Path) -> set[str]:
    """The stem of every figure M7 wrote."""
    return {p.name.removesuffix("_white.pdf")
            for p in (root / "outputs" / "07_proteo").glob("*_v1_white.pdf")}  # fmt: skip


def _redraws_alike(root: Path, found, tmp_path: Path) -> None:
    """Each spec redrawn in its own clean folder equals M7's files and its page's words."""
    module_csvs = [next(root.glob(f"outputs/07_proteo/{spec.stem}.csv")) for spec in found]
    for spec, module_csv in zip(found, module_csvs, strict=True):
        # specs() runs after the module: an input written after its figure is found here too,
        # so the module's own page must carry the panel the spec was found for when drawn.
        assert 'id="ctx-repro-' in unpacked(module_csv.with_suffix(".html").read_text()), spec.stem
    for spec, module_csv, csv in zip(found, module_csvs, redrawn_each(found, root, tmp_path),
                                     strict=True):  # fmt: skip
        assert drawn_alike(csv, module_csv) == [], spec.stem
        said = [
            WORDS.findall(unpacked(p.with_suffix(".html").read_text())) for p in (module_csv, csv)
        ]
        assert said[0] and said[0] == said[1], spec.stem


@pytest.mark.parametrize("which", ["across_root", "cross_root", "study_root", "made"])
def test_every_figure_has_a_spec_and_redraws_byte_for_byte(which, request, tmp_path):
    """Every figure M7 draws, on each design, has a spec and redraws byte for byte from its
    tables, its lead and reading guide too. Catches a table written after its figure (no spec),
    the drawn genes' rows written for another gene set (a NaN row redrawn), the annotation's
    columns reordered, the coverage tertiles cut on another abundance, the class test run on
    the other method, a contrast's correlation read from another column, and a cross-study n
    counted from the wrong level."""
    root, dataset_file = request.getfixturevalue(which)
    found = specs(root, M7, dataset_file)
    assert {s.stem for s in found} == drawn(root), _kinds(found)
    assert _kinds(found) == EXPECTED[which]
    _redraws_alike(root, found, tmp_path)


def test_a_run_with_other_options_redraws_with_them(across_root, tmp_path):
    """ctx proteo with every option moved: the figures that print one change, and every figure
    redraws byte for byte from its tables. Catches an option read from Options' defaults
    instead of the run's table: the share in F1's lead and F2's line, the class F3 draws alone,
    and the method of F3's class test. F4 draws D-012's lfc whatever --lfc-min says."""
    root, dataset_file = across_root
    moved = tmp_path / "root"
    writable_copy(root, moved, ignore=shutil.ignore_patterns("07_proteo"))
    ctx("--root", str(moved), "proteo", "--dataset", "fake_brca", "--dataset-file",
        str(dataset_file), "--method", "pearson", "--min-quantified", "0.5", "--lfc-min", "0.5",
        "--gene-class", "complex")  # fmt: skip
    found = specs(moved, M7, dataset_file)
    assert _kinds(found) == EXPECTED["across_root"]
    changed = set()
    for spec in found:
        after = next(moved.glob(f"outputs/07_proteo/{spec.stem}.csv"))
        before = next(root.glob(f"outputs/07_proteo/*_{spec.stem[11:]}.csv"))
        if _read(after) != _read(before):
            changed.add(parse(spec.figure_id)[2])
    assert {"quantified", "coverage", "rho_by_class"} <= changed, changed
    _redraws_alike(moved, found, tmp_path / "redrawn")


def _read(csv: Path) -> tuple[bytes, list[str]]:
    """A figure's data table and its page's lead and reading guide."""
    return csv.read_bytes(), WORDS.findall(unpacked(csv.with_suffix(".html").read_text()))


def test_drawn_rows_hold_every_heatmap_gene_past_the_cap(across_root, tmp_path):
    """On 40 made-up genes, more than twice the heatmap's top and bottom count, the matrices
    read back from the drawn rows draw the same heatmap as the full ones. Catches the drawn
    genes cut at another count than the page's, which the 12-gene roots cannot show."""
    _, dataset_file = across_root
    spec = load(dataset_file)
    rng = np.random.default_rng(1)
    names = [f"G{i}" for i in range(1, 41)]
    genes = pl.DataFrame({"gene": names, "tested": [True] * 40, "rho": rng.uniform(-1, 1, 40)})
    rna, protein = rng.normal(0, 1, (40, 6)), rng.normal(0, 1, (40, 6))
    annotation = pl.DataFrame({"sample": [f"S{j}" for j in range(6)], "plex": ["P1"] * 6})
    full = Measured(Layer((), (), np.empty((0, 0))), [], genes, np.empty(0))
    full.extra |= {"rna": rna, "values": protein, "annotation": annotation}
    paths = {"proteo_genes": tmp_path / "genes.parquet"}
    genes.write_parquet(paths["proteo_genes"])
    names_of = {v: k for k, v in TABLES.items()}
    assert spec.proteo is not None
    for what, frame in figure_tables._correlation(full, spec.proteo, PAIRED_TOP):
        paths[names_of[what]] = tmp_path / f"{what}.parquet"
        frame.write_parquet(paths[names_of[what]])
    back = figure_tables.measured_from(paths)
    chosen = figures.heatmap_genes(genes, [], PAIRED_TOP)
    assert len(chosen) == 2 * PAIRED_TOP
    drawn = []
    for m in (full, back):
        recorder = Recorder("07_proteo", spec.slug)
        figures.paired(recorder, (genes, m.extra["rna"], m.extra["values"]),
                       m.extra["annotation"], chosen, (("plex", "plex"),))  # fmt: skip
        drawn.append(recorder.one("paired_zscores")[0].frame)
    assert drawn[0].equals(drawn[1])


def test_a_cross_study_n_counts_each_level_of_the_contrast(cross_root, tmp_path):
    """The fold-change figure's n is counted per level from the protein samples' groups: with
    one reference sample relabelled tested, the caption reads 5 and 3. Catches the two levels'
    n swapped, which the made-up root's 4 and 4 cannot show."""
    root, dataset_file = cross_root
    spec = next(s for s in specs(root, M7, dataset_file) if s.stem.endswith("_fold_change_v1"))
    folder = copied(spec, root, tmp_path / "n")
    samples = folder / next(i.file for i in spec.inputs if i.name == "protein_samples")
    frame = pl.read_parquet(samples)
    assert frame["group"].to_list().count("A") == frame["group"].to_list().count("B") == 4
    first_b = frame["group"].to_list().index("B")
    frame.with_columns(group=pl.Series([("A" if i == first_b else g)
                                        for i, g in enumerate(frame["group"])])
                       ).write_parquet(samples)  # fmt: skip
    raw = (folder / "params.toml").read_text().splitlines()
    (folder / "params.toml").write_text(
        "\n".join(x for x in raw if not x.startswith("protein_samples = ") or ".parquet" in x)
        + "\n"
    )
    page = unpacked(redrawn(folder, spec).with_suffix(".html").read_text())
    assert "5 and 3" in page and "4 and 4" not in page


def test_a_worked_gene_the_dataset_file_does_not_name_is_refused(across_root, tmp_path):
    """A params file naming another worked gene is refused by name. Catches the gene check
    left out, which would fail later, in the module's own code, on a gene no table holds."""
    root, dataset_file = across_root
    spec = next(s for s in specs(root, M7, dataset_file) if parse(s.figure_id)[2] == "worked")
    folder = copied(spec, root, tmp_path / "worked")
    raw = (folder / "params.toml").read_text()
    assert raw.count('gene = "G9"') == 1
    (folder / "params.toml").write_text(raw.replace('gene = "G9"', 'gene = "G1"'))
    with pytest.raises(ParamsError, match="gene: 'G1' is not a \\[proteo\\] gene"):
        load_params(folder / "params.toml", spec.figure_id)


def test_the_normals_class_figure_counts_the_normals_read_back(tmp_path):
    """The normals' n, which their class figure prints, is read back from the worked samples'
    tissues: 3 tumours and 2 normals give 3 and 2. Catches the normals counted as the tumours,
    which the made-up root's 4 cases with both tissues cannot show."""
    paths = {
        n: tmp_path / f"{n}.parquet"
        for n in ("proteo_genes", "paired_annotation", "drawn_rows", "worked_samples")
    }
    pl.DataFrame({"gene": ["G1"]}).write_parquet(paths["proteo_genes"])
    pl.DataFrame({"sample": ["t1", "t2", "t3"]}).write_parquet(paths["paired_annotation"])
    pl.DataFrame(schema=figure_tables.ROW_SCHEMA).write_parquet(paths["drawn_rows"])
    pl.DataFrame(
        {"sample": ["t1", "t2", "t3", "n1", "n2"], "tissue": ["tumour"] * 3 + ["normal"] * 2}
    ).write_parquet(paths["worked_samples"])
    back = figure_tables.measured_from(paths)
    assert (back.extra["n_samples"], back.extra["n_normal"]) == (3, 2)
