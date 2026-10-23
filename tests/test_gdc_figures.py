"""Which figures the registry offers a dataset by its route (fetchngs, count matrix, GDC), and
`ctx figure` on every offered kind on each route's made-up root (R1.3b). Each test names the
mutation it catches."""

import gzip
import shutil
import subprocess
import sys
from datetime import date
from pathlib import Path

import pytest
from clinical_root import CLINICAL
from conftest import DE_SPEC, ctx_run
from gdc_fixture import SPEC, aliquot, build

from ctx import pipeline
from ctx.dataset import DatasetError, load
from ctx.figures import command, kinds
from ctx.figures.kind import NoFigureError
from ctx.figures.params import ParamsError, load_params, params_file
from ctx.figures.record import input_paths
from ctx.figures.spec import FigureSpec, Input, figure_id

CTX = Path(sys.executable).parent / "ctx"
DE = ('\n[[contrasts]]\nname = "tumour_normal"\nvariable = "condition"\nreference = "normal"\n'
      'tested = "tumour"\ndesign = "~donor + group"\npaired = true\nn = [{n}, {n}]\n'
      '\n[de]\nfeatured = "tumour_normal"\n')  # fmt: skip
QC = ["band_map", "composition", "fastqc", "linked_pca", "read_fates", "read_flow", "sex_check",
      "size_genes", "strip"]  # fmt: skip
# M3's kinds, the laptop route's and nf-core's (D-116), offered on every fetchngs route.
PROCESS = ["agreement", "aligned_mapping", "biotypes", "counters", "duplication",
           "featured_settings", "gene_coverage", "junction_annotation", "library_complexity",
           "mapping", "read_distribution", "rseqc_share", "rseqc_strand", "strand_assigned",
           "strand_share"]  # fmt: skip
CARDS = ["box", "contrasts", "exons", "immune", "junctions", "tissues"]
# The kinds offered per module, measured 2026-10-21 from pipeline.route_modules and each kind's
# own ``made``: one dataset file per route, the GDC one with contrasts, [de] and [clinical].
OFFERED = {
    "fetchngs_curated": {"gather": ["composition", "library_sizes"], "qc": QC, "quant": ["fates"],
                         "process": PROCESS,
                         "genecard": CARDS},
    "fetchngs_de": {"qc": QC, "process": PROCESS,
                    "quant": ["estimates", "fates"],
                    "de": ["box", "calls", "dispersion", "lfc", "ma", "pca", "phist", "volcano"],
                    "function": ["background", "gsea", "nes", "ora", "running"],
                    "genecard": [*CARDS[:2], "evidence_view", *CARDS[2:]]},
    "counts": {"qc": QC, "de": ["dispersion", "ma", "phist", "volcano"],
               "genecard": CARDS},
    "gdc": {"de": ["dispersion", "ma", "phist", "volcano"], "genecard": CARDS,
            "clinical": ["competing", "confounding", "er", "expression", "km_score", "power",
                         "scan_score"]},
}  # fmt: skip


def gdc_spec(tmp: Path) -> Path:
    """A GDC dataset file with contrasts, [de] and [clinical], and no [function]."""
    path = tmp / "fake_brca.toml"
    path.write_text(
        SPEC.format(protocol="unstranded", cdr='attribution = "tcga"') + CLINICAL + DE.format(n=2)
    )
    return path


def curated_spec(tmp: Path) -> Path:
    """A fetchngs dataset file with its own table of the fetchngs fixture's runs."""
    from test_own_table import dataset, table  # noqa: PLC0415 - tests/ helper

    from ctx.gather.fetchngs import read  # noqa: PLC0415

    return dataset(tmp / "conf", "mine.csv", table(read(fetchngs_fixture()).runs))


def fetchngs_fixture() -> Path:
    """The fetchngs outdir fixture gather curates."""
    from test_own_table import FIXTURE  # noqa: PLC0415 - tests/ helper

    return FIXTURE


def counts_spec(tmp: Path) -> Path:
    """A count-matrix dataset file with a contrast and [de], beside its sample table."""
    from synth import counts as synthetic  # noqa: PLC0415 - tests/ helper
    from test_gather_counts import DE_SPEC as MATRIX_SPEC  # noqa: PLC0415

    (tmp / "conf").mkdir(parents=True, exist_ok=True)
    synthetic(n_genes=400).samples.select("sample", "donor", "cell_type", "condition").write_csv(
        tmp / "conf" / "matde.tsv", separator="\t"
    )
    (tmp / "conf" / "matde.toml").write_text(MATRIX_SPEC)
    return tmp / "conf" / "matde.toml"


def offered(spec_path: Path) -> dict[str, list[str]]:
    """The kinds the registry offers the dataset file, by module."""
    spec = load(spec_path)
    found: dict[str, list[str]] = {}
    for k in kinds.KINDS:
        if kinds.offered(k, spec):
            found.setdefault(k.module, []).append(k.figure)
    return {m: sorted(set(f)) for m, f in found.items()}


def test_the_registry_offers_each_route_only_the_modules_it_runs(tmp_path):
    """Each route's dataset file is offered exactly the listed kinds. Catches the route rule
    left out: M2 qc, M3 process and M4 quant offered to a GDC dataset whose route runs none of
    them, M1's bars offered to a GDC sheet, which `ctx figure` drew into a ColumnNotFoundError,
    and M6's offered from [de] alone, whose background lesson stopped on a KeyError."""
    files = {"fetchngs_curated": curated_spec(tmp_path / "f"), "fetchngs_de": DE_SPEC,
             "counts": counts_spec(tmp_path / "c"), "gdc": gdc_spec(tmp_path)}  # fmt: skip
    routes = {name: pipeline.route(load(path)) for name, path in files.items()}
    assert routes == {"fetchngs_curated": "fetchngs", "fetchngs_de": "fetchngs",
                      "counts": "counts", "gdc": "gdc"}  # fmt: skip
    assert {name: offered(path) for name, path in files.items()} == OFFERED


def redraw_every_offered(root: Path, spec_path: Path, out: Path) -> dict[str, str]:
    """Each offered kind the root holds the inputs of, drawn as `ctx figure` draws it.

    One figure per kind, the first its dataset file makes, from a folder holding its inputs
    alone. Returns ``drawn``, the module's own refusal, or ``no inputs`` per kind.
    """
    spec = load(spec_path)
    found = {}
    for k in kinds.KINDS:
        if not kinds.offered(k, spec) or k.parsed is not None:
            continue
        params = k.made(spec)[0]
        name = f"{k.module}/{k.figure}"
        try:
            paths = input_paths(k, params, spec, root)
        except DatasetError:
            found[name] = "no inputs"
            continue
        folder = out / f"{k.module}_{k.figure}"
        folder.mkdir(parents=True)
        for path in paths.values():
            shutil.copyfile(path, folder / path.name)
        fid = figure_id("breast", k.module, k.figure)
        listed = tuple(Input(n, p.name, "", "path", "") for n, p in paths.items())
        drawn = FigureSpec(fid, spec.slug, k.module, params, listed, k.plotting, *[""] * 7)
        (folder / "params.toml").write_text(params_file(drawn).split("[sha256]")[0])
        try:
            command.draw(load_params(folder / "params.toml", fid), folder / "drawn", date.today())
        except NoFigureError as error:
            found[name] = f"refused: {error}"
            continue
        found[name] = "drawn"
    return found


def _run(root: Path, spec: Path, *commands: tuple[str, ...]) -> None:
    """Each command on the root, for the dataset file; a failure stops with its stderr."""
    base = ("--dataset", load(spec).slug, "--dataset-file", str(spec))
    for one in commands:
        done = ctx_run("--root", str(root), one[0], *base, *one[1:], cwd=root.parent)
        assert done.returncode == 0, (one, done.stderr[-3000:])


def _fetchngs_curated(tmp: Path) -> tuple[Path, Path]:
    spec, root = curated_spec(tmp), tmp / "root"
    root.mkdir()
    _run(root, spec, ("gather", "--fetchngs-dir", str(fetchngs_fixture())))
    return root, spec


def _fetchngs_de(tmp: Path) -> tuple[Path, Path]:
    from synth import DE_GROUPS, design_counts, write_gene_sets  # noqa: PLC0415
    from test_qc_run import plant  # noqa: PLC0415

    root = tmp / "root"
    plant(root)
    write_gene_sets(root / load(DE_SPEC).resource("sets").dest,
                    synthetic=design_counts(DE_GROUPS))  # fmt: skip
    _run(root, DE_SPEC, ("qc",), ("de",), ("function",))
    return root, DE_SPEC


def _counts(tmp: Path) -> tuple[Path, Path]:
    from synth import counts as synthetic  # noqa: PLC0415
    from test_gather_counts import matrix  # noqa: PLC0415

    made, spec, root = synthetic(n_genes=400), counts_spec(tmp), tmp / "root"
    ids = [f"ENSG{i + 1:011d}.1" for i in range(made.counts.height)]
    (root / "resources").mkdir(parents=True)
    lines = [f'chr1\tt\tgene\t{i + 1}\t{i + 9}\t.\t+\t.\tgene_id "{g}"; gene_type "protein_coding";'
             f' gene_name "N{i}";' for i, g in enumerate(ids)]  # fmt: skip
    (root / "resources" / "a.gtf.gz").write_bytes(gzip.compress(("\n".join(lines) + "\n").encode()))
    path = matrix(tmp / "m.csv", ids, made.counts.drop("gene_id").to_numpy().astype(float),
                  sep=",", samples=made.counts.columns[1:])  # fmt: skip
    _run(root, spec, ("gather", "--counts", str(path)), ("qc",), ("de",))
    return root, spec


def _gdc_clinical(tmp: Path) -> tuple[Path, Path]:
    from clinical_root import make  # noqa: PLC0415

    root, spec = make(tmp)
    spec.write_text(spec.read_text() + DE.format(n=5))
    _run(root, spec, ("de",), ("clinical",))
    return root, spec


def _gdc_proteo(tmp: Path) -> tuple[Path, Path]:
    from proteo_root import across  # noqa: PLC0415

    return across(tmp), tmp / "fake_brca.toml"


def _gdc_ncrna(tmp: Path) -> tuple[Path, Path]:
    from ncrna_tcga_root import make  # noqa: PLC0415

    root, spec = make(tmp, date.today())
    _run(root, spec, ("ncrna",))
    return root, spec


NONE, DRAWN = "no inputs", "drawn"
# The module's own NoFigureError: it drew no such figure from the options it wrote.
REFUSED = "refused: {}: the module draws no such figure from these inputs"
READS_QC = {f"qc/{f}": NONE for f in ("read_fates", "read_flow", "fastqc")}
QC_RUN = {"qc/band_map": REFUSED.format("band_map"), "qc/linked_pca": DRAWN,
          "qc/size_genes": DRAWN, "qc/strip": REFUSED.format("strip_total_sequences"),
          **READS_QC}  # fmt: skip
PROCESS_NONE = {f"process/{f}": NONE for f in OFFERED["fetchngs_de"]["process"]}
# Measured 2026-10-21: every kind offered on each root, drawn, refused by its module, or with
# inputs the root does not hold (qc and process were not run on the gather root; process,
# quant, the reads half of qc and genecard's evidence view on no root here).
REDRAWN = {
    "fetchngs_curated": {"gather/composition": DRAWN, "gather/library_sizes": DRAWN,
                         **{f"qc/{f}": NONE for f in QC}, **PROCESS_NONE, "quant/fates": NONE},
    "fetchngs_de": {**QC_RUN, "qc/composition": REFUSED.format("composition"),
                    "qc/sex_check": DRAWN, **PROCESS_NONE, "quant/estimates": NONE,
                    "quant/fates": NONE,
                    **{f"de/{f}": DRAWN for f in OFFERED["fetchngs_de"]["de"]},
                    **{f"function/{f}": DRAWN for f in OFFERED["fetchngs_de"]["function"]},
                    "genecard/evidence_view": NONE},
    "counts": {**QC_RUN, "qc/composition": DRAWN, "qc/sex_check": NONE,
               **{f"de/{f}": DRAWN for f in OFFERED["counts"]["de"]}},
    "gdc_clinical": {**{f"de/{f}": DRAWN for f in OFFERED["gdc"]["de"]},
                     **{f"clinical/{f}": DRAWN for f in OFFERED["gdc"]["clinical"]}},
    "gdc_proteo": {f"proteo/{f}": DRAWN for f in ("contrast", "coverage", "paired_zscores",
                                                  "quantified", "rho_by_class", "worked")},
    "gdc_ncrna": {**{f"de/{f}": NONE for f in OFFERED["gdc"]["de"]},
                  **{f"ncrna/{f}": DRAWN for f in ("biotype_shares", "filter_by_biotype",
                                                   "lncrna_neighbour", "mirna_target_pairs",
                                                   "mirna_volcano", "target_shift",
                                                   "volcano_by_biotype")},
                  "ncrna/library_types": REFUSED.format("library_types")},
}  # fmt: skip
BUILDERS = {"fetchngs_curated": _fetchngs_curated, "fetchngs_de": _fetchngs_de,
            "counts": _counts, "gdc_clinical": _gdc_clinical, "gdc_proteo": _gdc_proteo,
            "gdc_ncrna": _gdc_ncrna}  # fmt: skip


@pytest.mark.parametrize("name", list(BUILDERS), ids=list(BUILDERS))
def test_every_offered_kind_draws_on_its_routes_root(tmp_path, name):
    """On each route's made-up root, every offered kind whose inputs the root holds draws, or
    its module refuses it by name; none stops on an error. Catches a kind offered to a route
    whose sheet lacks its columns, as M1's bars on a GDC sheet were."""
    root, spec = BUILDERS[name](tmp_path)
    assert redraw_every_offered(root, spec, tmp_path / "redrawn") == REDRAWN[name]


@pytest.fixture(scope="module")
def gathered(tmp_path_factory) -> tuple[Path, Path]:
    """A made-up GDC project gathered into a root: the root's gather folder and the dataset file."""
    tmp = tmp_path_factory.mktemp("gdc_figures")
    build(tmp / "gdc", [aliquot(1), aliquot(2), aliquot(1, code="11"), aliquot(2, code="11")])
    spec = gdc_spec(tmp)
    done = subprocess.run([str(CTX), "--root", str(tmp / "root"), "gather", "--dataset",
                           "fake_brca", "--gdc-dir", str(tmp / "gdc"), "--dataset-file", str(spec)],
                          capture_output=True, text=True, check=False)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    return tmp / "root" / "outputs" / "01_gather", spec


@pytest.mark.parametrize("figure", ["composition", "library_sizes"])
def test_ctx_figure_refuses_gathers_figures_for_a_gdc_sheet(gathered, tmp_path, figure):
    """`ctx figure breast_m1_<figure>` on a GDC sheet exits 1 naming the figure, with no
    traceback. Catches the refusal left out of the params check, which reaches the drawer and
    stops on a ColumnNotFoundError for cell_type or read_count."""
    out, spec = gathered
    sheet = next(out.glob("*_fake_brca_samplesheet_v1.parquet"))
    (tmp_path / sheet.name).write_bytes(sheet.read_bytes())
    (tmp_path / spec.name).write_bytes(spec.read_bytes())
    registry = "" if figure == "composition" else 'qc_registry = "qc_metrics.toml"\n'
    if registry:
        (tmp_path / "qc_metrics.toml").write_bytes(
            (Path(__file__).resolve().parents[1] / "conf" / "qc_metrics.toml").read_bytes()
        )
    params = tmp_path / "params.toml"
    params.write_text(
        f'figure = "breast_m1_{figure}"\ndataset = "fake_brca"\n\n[inputs]\n'
        f'samplesheet = "{sheet.name}"\ndataset_file = "{spec.name}"\n{registry}'
    )
    done = subprocess.run([str(CTX), "figure", f"breast_m1_{figure}", "--params", str(params),
                           "--out", str(tmp_path / "drawn")], capture_output=True, text=True,
                          check=False)  # fmt: skip
    assert done.returncode == 1, done.stdout
    assert f"{figure}: fake_brca's dataset file makes no such figure" in done.stderr, done.stderr
    assert "Traceback" not in done.stderr and "Error:" not in done.stderr, done.stderr[-2000:]
    assert not (tmp_path / "drawn").exists()


def test_m6s_background_without_function_is_refused_by_name(tmp_path):
    """A params file for M6's background lesson of a dataset file with [de] and no [function]
    is refused before drawing. Catches M6's kinds made from [de] alone, whose inputs
    then omit the gene sets the drawer reads."""
    spec = gdc_spec(tmp_path)
    (tmp_path / "results.parquet").write_bytes(b"")
    params = tmp_path / "params.toml"
    params.write_text('figure = "breast_m6_background"\ndataset = "fake_brca"\n'
                      'contrast = "tumour_normal"\n\n[inputs]\nresults = "results.parquet"\n'
                      f'dataset_file = "{spec.name}"\n')  # fmt: skip
    with pytest.raises(ParamsError, match="background: fake_brca's dataset file makes no such"):
        load_params(params, "breast_m6_background")
