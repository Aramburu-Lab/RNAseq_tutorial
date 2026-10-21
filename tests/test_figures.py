"""Phase 2D: figure IDs, FigureSpec and `ctx figure` (ctx.figures). Each test names its mutation.

A figure is redrawn in a clean folder from the files its spec names, copied out of the session
pipeline's root, and its data table must equal the module's byte for byte.
"""

import shutil
import tomllib
from pathlib import Path

import polars as pl
import pytest
from conftest import DE_SPEC, ctx_run
from test_own_table import FIXTURE, dataset, table

from ctx.dataset import load
from ctx.figures import kinds, record
from ctx.figures.params import ParamsError, load_params, params_file
from ctx.figures.spec import figure_id, parse
from ctx.gather.fetchngs import read
from ctx.report.builder import unpacked

ROOT = Path(__file__).resolve().parents[1]
# Every kind the synthetic dataset file makes: the featured contrast's volcano, MA plot,
# p histogram, dispersions, GSEA and ORA bars and background lesson; both declared scatters;
# the one lesson; both PCAs; the named genes' boxes; the featured set; the running score.
EXPECTED = {"volcano", "ma", "phist", "lfc", "calls", "gsea", "ora", "nes", "pca", "box",
            "dispersion", "background", "running"}  # fmt: skip


def test_ids_round_trip():
    """An ID reads back as its parts. Catches a module number taken from the output folder's
    NN (M9 is 11_genecard, M10 12_clinical) and a cohort that is not one word."""
    assert figure_id("breast", "de", "volcano") == "breast_m5_volcano"
    assert parse("breast_m5_volcano") == ("breast", "de", "volcano")
    assert parse("kidney_m10_km_os") == ("kidney", "clinical", "km_os")
    assert parse("myeloid_m9_tissues") == ("myeloid", "genecard", "tissues")
    for bad in ("Breast_m5_volcano", "breast_m11_x", "breast_5_volcano", "breast_m5_", "m5_x"):
        assert parse(bad) is None, bad


def specs(root: Path, folders=("05_de", "06_function"), dataset_file: Path = DE_SPEC) -> list:
    """The spec of every figure the modules of ``folders`` wrote, found as each module found it."""
    load(dataset_file)  # the dataset file this process read, as the drawing process had
    found = []
    for folder in folders:
        for pdf in sorted((root / "outputs" / folder).glob("*_v1_white.pdf")):
            csv = pdf.with_name(pdf.name.replace("_white.pdf", ".csv"))
            if (spec := record.describe(csv.name[:-4], csv.parent, csv)) is not None:
                found.append(spec)
    return found


def copied(spec, root: Path, folder: Path) -> Path:
    """A clean folder holding the spec's inputs and its params file."""
    folder.mkdir(parents=True)
    for item in spec.inputs:
        top = {"root": root, "repo": ROOT}.get(item.base)
        source = top / item.source if top else Path(item.source).expanduser()
        shutil.copyfile(source, folder / item.file)
    (folder / "params.toml").write_text(params_file(spec))
    return folder


def redrawn(folder: Path, spec) -> Path:
    done = ctx_run("figure", spec.figure_id, "--params", "params.toml", cwd=folder)
    assert done.returncode == 0, (spec.figure_id, done.stderr[-3000:])
    (csv,) = folder.glob(f"*_{spec.stem[11:]}.csv")
    return csv


# Redraws run this many at once: each ``ctx figure`` is a fresh process, most of it start-up,
# and one at a time the redraw tests were the suite's longest (M3.14).
REDRAWS = 4


def redrawn_each(found: list, root: Path, folder: Path) -> list[Path]:
    """``redrawn`` of each spec, its inputs copied from ``root`` into ``folder/<stem>``,
    ``REDRAWS`` at a time: the redrawn CSVs in ``found``'s order."""
    from concurrent.futures import ThreadPoolExecutor  # noqa: PLC0415

    with ThreadPoolExecutor(REDRAWS) as pool:
        return list(pool.map(lambda s: redrawn(copied(s, root, folder / s.stem), s), found))


def drawn_alike(redrawn_csv: Path, module_csv: Path) -> list[str]:
    """The data table, and the PDF and SVG of each theme, that differ from the module's; and
    ``panel`` where the module's own figure page carries no Reproduce panel.

    Both sets are named by the same stem; a redraw dated another day differs in name alone.
    A spec found after the run, as ``specs`` finds it, can see a table the module wrote after
    the figure; the panel on the module's page was made when the figure was written.
    """
    tail = (".csv", *(f"_{t}.{e}" for t in ("white", "black") for e in ("pdf", "svg")))
    new, old = redrawn_csv.name[: -len(".csv")], module_csv.name[: -len(".csv")]
    differ = [end for end in tail
              if (redrawn_csv.parent / f"{new}{end}").read_bytes()
              != (module_csv.parent / f"{old}{end}").read_bytes()]  # fmt: skip
    page = unpacked((module_csv.parent / f"{old}.html").read_text())
    return differ if 'id="ctx-repro-' in page else [*differ, "panel"]


def test_every_spec_redraws_its_data_byte_for_byte(pipeline, tmp_path):
    """Each figure with a spec, redrawn from its inputs alone in a clean folder, writes the
    module's data table byte for byte, and the module's README names its ID. Catches a fit read
    back with another row order, the samples taken from the contrasts table, the gene symbols
    joined twice, a figure whose README lost its ID, and a kind no module figure matched."""
    root = pipeline[0]
    found = specs(root)
    assert {parse(s.figure_id)[2] for s in found} == EXPECTED
    for spec, csv in zip(found, redrawn_each(found, root, tmp_path), strict=True):
        module_csv = next(root.glob(f"outputs/*/{spec.stem}.csv"))
        readme = module_csv.with_name(f"{spec.stem}_README.md").read_text()
        assert f"`{spec.figure_id}`" in readme, spec.stem
        assert drawn_alike(csv, module_csv) == [], spec.stem


def test_a_lesson_redrawn_at_another_setting_counts_as_the_grid(pipeline, tmp_path):
    """The lesson's bars at raw p < 0.01 with no fold-change floor equal M5's own grid table
    at that setting. Catches a setting read but not passed to the bars (the build's counts
    again) and a basis read as another column."""
    root = pipeline[0]
    spec = next(s for s in specs(root) if parse(s.figure_id)[2] == "calls")
    folder = copied(spec, root, tmp_path / "pvalue")
    params = (folder / "params.toml").read_text()
    for old, new in (('basis = "padj"', 'basis = "pvalue"'), ("alpha = 0.05", "alpha = 0.01"),
                     ("lfc_min = 1.0", "lfc_min = 0.0")):  # fmt: skip
        assert params.count(old) == 1, old
        params = params.replace(old, new)
    (folder / "params.toml").write_text(params)
    drawn = pl.read_csv(redrawn(folder, spec))
    (grid,) = (root / "outputs" / "05_de").glob("*_calls_grid_batch_v1.parquet")
    at = pl.read_parquet(grid).filter((pl.col("basis") == "pvalue") & (pl.col("alpha") == 0.01)
                                      & (pl.col("lfc_min") == 0.0))  # fmt: skip
    for row in at.iter_rows(named=True):
        for call in ("up", "down"):
            shown = drawn.filter((pl.col("contrast") == row["contrast"]) & (pl.col("call") == call))
            assert shown["genes"].to_list() == [row[call]], (row, call)
    build = pl.read_csv(next(root.glob(f"outputs/05_de/{spec.stem}.csv")))
    assert drawn["genes"].to_list() != build["genes"].to_list()


def test_a_bad_params_file_lists_every_problem(pipeline, tmp_path):
    """One run names each problem: a changed input, an unknown key, a wrong type, an undeclared
    contrast, a missing input, and the ID the command named against the file's. Catches a
    loader that stops at the first, and an input whose bytes changed drawn anyway."""
    root = pipeline[0]
    spec = next(s for s in specs(root) if parse(s.figure_id)[2] == "volcano")
    folder = copied(spec, root, tmp_path / "bad")
    raw = (folder / "params.toml").read_text()
    raw = raw.replace(f'contrast = "{spec.params.contrast}"', 'contrast = "nope"\ncolour = 3')
    raw = raw.replace('dataset = "desynth"', "dataset = 7")
    raw = raw.replace("size_factors = ", "# size_factors = ", 1)
    (folder / "params.toml").write_text(raw)
    results = next(i for i in spec.inputs if i.name == "results")
    (folder / results.file).write_bytes((folder / results.file).read_bytes() + b"\0")
    with pytest.raises(ParamsError) as error:
        load_params(folder / "params.toml", spec.figure_id.replace("volcano", "ma"))
    said = str(error.value)
    for part in (
        "figure: the file names",
        "colour: not a parameter",
        "dataset: expected str",
        "inputs.results:",
        "has SHA-256",
        "sha256.size_factors: names no input",
    ):
        assert part in said, (part, said)
    raw = raw.replace("dataset = 7", 'dataset = "desynth"')
    (folder / "params.toml").write_text(raw)
    (folder / results.file).write_bytes((folder / results.file).read_bytes()[:-1])
    with pytest.raises(ParamsError) as error:
        load_params(folder / "params.toml", spec.figure_id)
    said = str(error.value)
    for part in ("contrast: 'nope' is not a contrast", "inputs.size_factors: missing"):
        assert part in said, (part, said)


def test_figure_refuses_a_folder_inside_a_ctx_root(pipeline, tmp_path):
    """ctx figure writes nothing inside a ctx root, at its top or in a module folder. Catches the
    redraw written over a module's own figure files, whose names it shares on the same day."""
    root = pipeline[0]
    spec = next(s for s in specs(root) if parse(s.figure_id)[2] == "phist")
    folder = copied(spec, root, tmp_path / "ok")
    inside = shutil.copytree(root / "outputs" / "05_de", tmp_path / "copy" / "outputs" / "05_de")
    for out in (inside, inside.parent.parent):
        before = sorted(p.name for p in out.rglob("*"))
        done = ctx_run("figure", spec.figure_id, "--params", str(folder / "params.toml"),
                       "--out", str(out), cwd=folder)  # fmt: skip
        assert done.returncode == 1 and "outside any root" in done.stderr, done.stderr
        assert sorted(p.name for p in out.rglob("*")) == before


def test_gather_figures_redraw_from_the_sheet(tmp_path):
    """M1's composition and depth bars, redrawn from the sheet, the dataset file and the QC
    registry, equal gather's data byte for byte. Catches the depth standard read from the
    repository's registry instead of the copy named, and a sheet read back with other types."""
    root, conf = tmp_path / "root", tmp_path / "conf"
    root.mkdir()
    spec_path = dataset(conf, "mine.csv", table(read(FIXTURE).runs))
    done = ctx_run("--root", str(root), "gather", "--dataset", "mine", "--dataset-file",
                   str(spec_path), "--fetchngs-dir", str(FIXTURE), cwd=root)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    load(spec_path)
    found = [s for csv in sorted((root / "outputs" / "01_gather").glob("*_v1.csv"))
             if (s := record.describe(csv.name[:-4], csv.parent, csv))]  # fmt: skip
    assert sorted(parse(s.figure_id)[2] for s in found) == ["composition", "library_sizes"]
    for spec in found:
        folder = copied(spec, root, tmp_path / spec.stem)
        module_csv = root / "outputs" / "01_gather" / f"{spec.stem}.csv"
        assert redrawn(folder, spec).read_bytes() == module_csv.read_bytes(), spec.stem
    # The registry named is the one read: ENCODE's good depth moved in the copy moves the lead.
    spec = next(s for s in found if s.figure_id.endswith("library_sizes"))
    folder = copied(spec, root, tmp_path / "registry")
    registry = folder / next(i.file for i in spec.inputs if i.name == "qc_registry")
    raw = registry.read_text()
    assert tomllib.loads(raw)["metrics"]["total_sequences"]["plot_bands"][0] == 30_000_000
    registry.write_text(raw.replace("plot_bands = [30000000.0,", "plot_bands = [40000000.0,", 1))
    params = (folder / "params.toml").read_text().splitlines()
    # The registry's hash leaves the params file; its file name, which holds a dot, stays.
    kept = [p for p in params if not p.startswith("qc_registry = ") or "." in p]
    (folder / "params.toml").write_text("\n".join(kept) + "\n")
    page = unpacked(redrawn(folder, spec).with_suffix(".html").read_text())
    assert "fewer than 40 million" in page and "fewer than 30 million" not in page


def test_every_kind_names_inputs_its_drawer_reads():
    """Each kind's drawer and plotting function exist, and its parameters are a frozen dataclass.
    Catches a kind pointing at a renamed function, which only a redraw would find."""
    import importlib  # noqa: PLC0415

    for kind in kinds.KINDS:
        for target in (kind.draw, kind.plotting):
            module, function = target.split(":")
            assert callable(getattr(importlib.import_module(module), function)), target
        assert kind.params.__dataclass_params__.frozen, kind.figure


def test_figure_options_have_cards():
    """Every option of ctx figure, and the global --cohort, has a card. Catches an option added
    without one, which no module page would show missing, since ctx figure writes none."""
    import typer  # noqa: PLC0415

    from ctx.cli import app  # noqa: PLC0415
    from ctx.report.cards import parameter_cards  # noqa: PLC0415

    group = typer.main.get_command(app)
    figure = group.commands["figure"]  # type: ignore[attr-defined]
    names = [c.name for c in parameter_cards(group, figure)]
    assert {"--cohort", "figure_id", "--params", "--out"} <= set(names), names


def test_a_pca_redraw_reads_the_sheet_through_the_gdc_door(pipeline, tmp_path):
    """A GDC sheet's PCA, redrawn, holds the chosen files alone: an aliquot gather did not take,
    its counts present, stays out as M5 keeps it out. Catches the drawer reading the sheet
    without ``gdc.chosen``, which would draw that aliquot as one more sample."""
    from ctx.figures import draw_m5  # noqa: PLC0415

    root = pipeline[0]
    spec = next(s for s in specs(root) if s.figure_id.endswith("_pca") and s.params.pca == "all")
    folder = copied(spec, root, tmp_path / "gdc")
    files = {i.name: folder / i.file for i in spec.inputs}
    sheet = pl.read_parquet(files["samplesheet"])
    extra = sheet.head(1).with_columns(sample=pl.lit("S_EXTRA"), chosen=pl.lit(False))
    pl.concat([sheet.with_columns(chosen=pl.lit(True)), extra]).write_parquet(files["samplesheet"])
    counts = pl.read_parquet(files["counts"])
    counts.with_columns(S_EXTRA=pl.col(sheet["sample"][0])).write_parquet(files["counts"])
    drawn, _ = draw_m5.pca(files, spec.params, load(files["dataset_file"]))
    module = pl.read_csv(next(root.glob(f"outputs/05_de/{spec.stem}.csv")), comment_prefix="#")
    assert sorted(drawn.frame["sample"]) == sorted(module["sample"])
    assert "S_EXTRA" not in drawn.frame["sample"].to_list()


def test_a_pca_the_dataset_file_does_not_declare_is_refused(pipeline, tmp_path):
    """A params file naming a PCA the dataset file lacks is refused with the name. Catches the
    kind's check left out, which would stop with a StopIteration inside the drawer."""
    root = pipeline[0]
    spec = next(s for s in specs(root) if s.figure_id.endswith("_pca"))
    folder = copied(spec, root, tmp_path / "pca")
    raw = (folder / "params.toml").read_text()
    assert raw.count(f'pca = "{spec.params.pca}"') == 1
    (folder / "params.toml").write_text(raw.replace(f'pca = "{spec.params.pca}"', 'pca = "nope"'))
    with pytest.raises(ParamsError, match=r"pca: 'nope' is not a \[\[de.pca\]\] name"):
        load_params(folder / "params.toml", spec.figure_id)


def test_a_kind_parsed_from_its_description_must_read_back(monkeypatch):
    """A kind the dataset file cannot list matches a description its parser reads back to the
    same description, and no other. Catches the read-back check left out, which would give a
    parser that strips a prefix every stem it was handed."""
    from dataclasses import dataclass  # noqa: PLC0415

    from ctx.figures.kind import Kind  # noqa: PLC0415

    @dataclass(frozen=True)
    class Card:
        dataset: str
        gene: str

    def parsed(what: str, spec) -> Card:
        return Card(spec.slug, what.removeprefix("card_"))

    card = Kind("genecard", "card", Card, "m:f", "m:f", lambda _: [], lambda p: f"card_{p.gene}",
                lambda p, s: {}, parsed=parsed)  # fmt: skip
    monkeypatch.setattr(kinds, "KINDS", (card,))
    spec = load(DE_SPEC)
    assert kinds.matched("genecard", "card_nptx2", spec) == (card, Card("desynth", "nptx2"))
    assert kinds.matched("genecard", "tissues_nptx2", spec) is None
