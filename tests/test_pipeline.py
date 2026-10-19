"""U.6: `ctx run` and `ctx run --check` (ctx.pipeline).

Datasets are written at test time over the metadata-only fetchngs fixture with an own table
(tests/test_own_table.py), so no accession is typed here. Each test names its mutation.
"""

import gzip
import html
import re
import subprocess
from dataclasses import replace
from pathlib import Path

import polars as pl
import pytest
import typer
from conftest import ctx_run
from gdc_fixture import PROJECT, aliquot, build
from test_gdc_projects import OTHER, folders, two_projects
from test_own_table import FIXTURE, SPEC, dataset, table

from ctx import pipeline
from ctx.cli import app
from ctx.dataset import load
from ctx.gather.fetchngs import read
from ctx.runtime import launcher

ROOT = Path(__file__).resolve().parents[1]
# Commands that are not modules a dataset runs: the machine, the read-level Track A, and these.
# genecard fetches from HPA and Ensembl, which one command must not do unasked; reads lays out
# a student's FASTQ files before nf-core/rnaseq, outside any ctx root; figure redraws one figure
# in the reader's own folder, from a params file a page wrote.
NOT_RUN = {"setup", "process", "genecard", "run", "guide", "reads", "figure"}
CONTRAST = """
[[contrasts]]
name = "case_control"
variable = "condition"
reference = "control"
tested = "case"
n = [{reference}, {tested}]
{extra}"""


@pytest.fixture
def runs():
    return read(FIXTURE).runs


def _spec(folder: Path, runs: pl.DataFrame, extra: str = "") -> Path:
    spec = dataset(folder, "mine.csv", table(runs))
    spec.write_text(spec.read_text() + extra)
    return spec


def test_every_command_has_a_place():
    """Every subcommand is a module ctx run places in ORDER, or one it never runs. Catches a new
    module command left out of ctx run (M8 and M9 join ORDER in their place)."""
    commands = set(typer.main.get_command(app).commands)
    assert commands == set(pipeline.ORDER) | NOT_RUN, commands ^ (set(pipeline.ORDER) | NOT_RUN)


def test_steps_in_course_order(tmp_path, runs):
    """The modules a dataset declares run in course order, each with its own options. Catches
    two modules swapped (the test reads the order, not the set), quant without an annotation,
    and function without contrasts."""
    featured = '\n[de]\nfeatured = "case_control"\n'
    spec_path = _spec(tmp_path, runs, CONTRAST.format(reference=7, tested=7, extra=featured))
    spec = load(spec_path)
    rnaseq = tmp_path / "rnaseq"
    options = pipeline.Options("mine", spec_path, FIXTURE, rnaseq)
    assert [m for m, _ in pipeline.steps(spec, options)] == ["gather", "qc", "de", "report"]
    gtf = replace(spec, resources=(*spec.resources, replace_resource()))
    full = replace(gtf, function=object())  # type: ignore[arg-type]
    planned = pipeline.steps(full, options)
    assert [m for m, _ in planned] == ["gather", "qc", "quant", "de", "function", "report"]
    gather = dict(planned)["gather"]
    assert gather[gather.index("--rnaseq-dir") + 1] == str(rnaseq)
    no_counts = pipeline.steps(full, replace(options, rnaseq_dir=None))
    assert [m for m, _ in no_counts] == ["gather", "qc", "report"]
    gdc = replace(full, clinical=object())  # type: ignore[arg-type]
    planned = pipeline.steps(gdc, pipeline.Options("mine", spec_path, gdc_dirs=(tmp_path,)))
    assert [m for m, _ in planned] == ["gather", "de", "function", "clinical", "report"]
    both = replace(gdc, ncrna=object())  # type: ignore[arg-type]
    planned = pipeline.steps(both, pipeline.Options("mine", spec_path, gdc_dirs=(tmp_path,)))
    assert [m for m, _ in planned] == ["gather", "de", "function", "ncrna", "clinical", "report"]
    # M7 reads gather's counts and, across studies, M5's results: after de, before ncrna.
    across = replace(both, proteo=replace_proteo("across_samples"))
    planned = pipeline.steps(across, pipeline.Options("mine", spec_path, gdc_dirs=(tmp_path,)))
    assert [m for m, _ in planned] == ["gather", "de", "function", "proteo", "ncrna", "clinical",
                                       "report"]  # fmt: skip
    no_de = replace(across, de=None, proteo=replace_proteo("cross_dataset"))
    gdc_options = pipeline.Options("mine", spec_path, gdc_dirs=(tmp_path,))
    assert "proteo" not in dict(pipeline.steps(no_de, gdc_options))


def replace_proteo(design):
    """A stand-in [proteo]: steps reads only its design."""
    from ctx.dataset import ProteoSpec  # noqa: PLC0415

    return ProteoSpec(design, "protein")


def replace_resource():
    """A stand-in annotation resource: steps reads only its name."""
    from ctx.dataset import Resource  # noqa: PLC0415

    return Resource("annotation_gtf", ("http://127.0.0.1/a.gtf.gz",), 1, "0" * 64, "a.gtf.gz")


def test_execute_stops_at_the_first_failure(tmp_path, monkeypatch, capsys):
    """Each module runs as `python -m ctx` with the root and profile first; the first failure
    stops the run and names it and what was not run. Catches running on past a failure."""
    calls = []

    def fake(argv, check):
        calls.append(argv)
        return subprocess.CompletedProcess(argv, 3 if "de" in argv else 0)

    monkeypatch.setattr(pipeline.subprocess, "run", fake)
    planned = [("gather", ("--dataset", "x")), ("de", ("--dataset", "x")), ("report", ())]
    assert pipeline.execute(tmp_path, "laptop", planned) == 3
    assert [c[3:8] for c in calls] == [
        ["--root", str(tmp_path), "--profile", "laptop", "gather"],
        ["--root", str(tmp_path), "--profile", "laptop", "de"],
    ]
    assert calls[0][1:3] == ["-m", "ctx"]
    assert "de failed with exit 3; not run: report" in capsys.readouterr().err


def test_check_lists_every_problem(tmp_path, runs):
    """An rnaseq outdir that is not there, unfetched resources, a wrong n and a missing
    [function] are listed together. Catches the check stopping at its first problem."""
    extra = CONTRAST.format(reference=1, tested=1, extra='\n[de]\nfeatured = "case_control"\n')
    spec_path = _spec(tmp_path, runs, extra)
    gene_sets = pipeline._resource(pipeline.GENE_SETS)
    spec_path.write_text(spec_path.read_text() + "\n" + gene_sets)
    root = tmp_path / "root"
    missing = tmp_path / "no_rnaseq"
    options = pipeline.Options("mine", spec_path, FIXTURE, missing)
    problems, planned = pipeline.check(root, options)
    texts = [p.text for p in problems]
    assert texts[0] == f"--rnaseq-dir {missing} is not a folder", texts
    assert any(t.startswith("resource msigdb_hallmark is not fetched") for t in texts), texts
    assert any(t.startswith("case_control: n = [1, 1] declared, [7, 7] found") for t in texts)
    assert any(t.startswith("no [function]") for t in texts) and len(texts) == 4, texts
    assert planned == ["gather", "qc", "de", "report"]


def test_suggested_lines_pass_once_pasted(tmp_path, runs):
    """Every line --check suggests, pasted into the file as printed, clears its problem, until
    none is left. Catches a suggested contrast whose n fails, a block that is not TOML, and a
    [function] naming a resource the file lacks."""
    spec_path = _spec(tmp_path, runs)
    root = tmp_path / "root"
    for _ in range(3):
        # Stand-ins where ctx setup would fetch: the check reads only that they are there.
        for resource in load(spec_path).resources:
            (root / resource.dest).parent.mkdir(parents=True, exist_ok=True)
            (root / resource.dest).write_text("stand-in")
        problems, _ = pipeline.check(root, pipeline.Options("mine", spec_path, FIXTURE))
        lines = [p.lines for p in problems if p.lines]
        assert len(lines) == len(problems) and not any(p.blocks for p in problems), problems
        if not problems:
            break
        spec_path.write_text(spec_path.read_text() + "\n" + "\n".join(lines))
    assert problems == []
    spec = load(spec_path)
    assert [c.expected_n for c in spec.contrasts] == [(7, 7)]
    assert spec.de is not None and spec.function is not None


def test_n_after_filters_and_replicates(tmp_path, runs):
    """A contrast's n counts what its filters and replicate ranking keep, as ctx de does. Catches
    n counted before the filters, and replicates counted twice."""
    body = table(runs).decode("utf-8-sig").replace("\r\n", "\n").splitlines()
    # Run 1 joins run 0's donor, cell type and condition: a replicate, ranked out.
    first = body[1].split(",", 1)[1]
    body[2] = body[2].split(",", 1)[0] + "," + first
    (tmp_path / "mine.csv").write_text("\n".join(body) + "\n")
    spec_path = tmp_path / "mine.toml"
    filters = 'filters = { cell_type = ["type B"] }'
    spec_path.write_text(SPEC.format(table="mine.csv")
                         + CONTRAST.format(reference=0, tested=0, extra=filters))  # fmt: skip
    problems, _ = pipeline.check(tmp_path / "root", pipeline.Options("mine", spec_path, FIXTURE))
    (found,) = [p.text for p in problems if p.text.startswith("case_control")]
    # Type B is every even run: runs 0, 2, 4, 6 control (run 1 now repeats run 0, ranked out) and
    # 8, 10, 12 case.
    assert "found; set n = [4, 3]" in found, found


def test_run_matches_the_module_commands(tmp_path, runs):
    """`ctx run` writes the tables the module commands write one by one, and each page keeps
    the command ctx run gave it. Catches a module run with other options than its own command
    would take, and a page whose Reproduce names ctx run instead."""
    spec_path = _spec(tmp_path / "conf", runs)
    one, many = tmp_path / "one", tmp_path / "many"
    for root in (one, many):
        root.mkdir()
    done = ctx_run("--root", str(one), "run", "--dataset", "mine", "--dataset-file",
                   str(spec_path), "--fetchngs-dir", str(FIXTURE), cwd=one)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    printed = re.findall(r"^ctx run: \d+ of 3: (.+)$", done.stdout, re.M)
    # Through the environment that ran it, so a plain shell can paste it (runtime.launcher).
    assert all(c.startswith(" ".join(launcher()) + " --root ") for c in printed), printed
    assert [c.split("--root ")[1].split()[1] for c in printed] == ["gather", "qc", "report"]
    # Each module as its own page tells a reader to run it, typed here, not taken from ctx run.
    named = ("--dataset", "mine", "--dataset-file", str(spec_path))
    for module, args in (("gather", (*named, "--fetchngs-dir", str(FIXTURE))), ("qc", named),
                         ("report", named)):  # fmt: skip
        done = ctx_run("--root", str(many), module, *args, cwd=many)
        assert done.returncode == 0, (module, done.stderr[-3000:])
    tables = {p.relative_to(one): p.read_bytes() for p in one.rglob("*.csv")}
    assert len(tables) > 3
    assert tables == {p.relative_to(many): p.read_bytes() for p in many.rglob("*.csv")}
    pages = [p.read_text() for p in sorted((one / "outputs").glob("*/*_overview_v1.html"))]
    shown = [html.unescape(re.search(r'id="ctx-command">(.*?)</code>', p, re.S).group(1))
             for p in pages]  # fmt: skip
    # gather's and qc's pages; report writes a hub, whose page shows no command.
    assert shown == printed[:2]


def test_check_finds_no_matrix_and_fixes_paste(tmp_path, runs):
    """An rnaseq outdir without its count matrix is a problem before setup, and the setup line
    a fix gives runs from a plain shell, through the environment. Catches the check that listed
    only the two resources while Part 2 had died (U.10), and a bare `ctx` fix line."""
    spec_path = _spec(tmp_path, runs)
    spec_path.write_text(spec_path.read_text() + "\n" + pipeline._resource(pipeline.GENE_SETS))
    rnaseq = tmp_path / "rnaseq"
    (rnaseq / "salmon").mkdir(parents=True)
    options = pipeline.Options("mine", spec_path, FIXTURE, rnaseq)
    texts = [p.text for p in pipeline.check(tmp_path / "root", options)[0]]
    assert any(t.startswith("--rnaseq-dir holds no count matrix") for t in texts), texts
    (fetch,) = [t for t in texts if t.startswith("resource ")]
    assert f"; run {' '.join(launcher())} --root " in fetch, fetch


def test_folders_carried_to_the_modules_that_take_them(tmp_path, runs):
    """--salmon-dir reaches quant and --mirna-dir with --isoform-dir reach gather on the GDC
    route, each parsed by that module's own command line. Catches quant planned without
    --salmon-dir, gather without the miRNA folders, and a flag the module does not take."""
    spec_path = _spec(tmp_path, runs)
    spec = replace(load(spec_path), resources=(*load(spec_path).resources, replace_resource()))
    salmon, mirna, isoform = tmp_path / "salmon", tmp_path / "mirna", tmp_path / "isoform"
    reads = pipeline.Options("mine", spec_path, FIXTURE, tmp_path / "rnaseq", salmon_dir=salmon)
    gdc = pipeline.Options("mine", spec_path, gdc_dirs=(tmp_path,), mirna_dir=mirna,
                           isoform_dir=isoform)  # fmt: skip
    quant = dict(pipeline.steps(spec, reads))["quant"]
    gather = dict(pipeline.steps(spec, gdc))["gather"]
    assert quant[quant.index("--salmon-dir") + 1] == str(salmon)
    assert gather[gather.index("--mirna-dir") + 1] == str(mirna)
    assert gather[gather.index("--isoform-dir") + 1] == str(isoform)
    commands = typer.main.get_command(app).commands  # type: ignore[attr-defined]
    for module, args in (*pipeline.steps(spec, reads), *pipeline.steps(spec, gdc)):
        parsed = commands[module].make_context(module, list(args)).params
        for flag, value in zip(args[::2], args[1::2], strict=True):
            found = parsed[flag.lstrip("-").replace("-", "_")]
            # A repeatable option (gather's --gdc-dir) parses to a tuple of each value given.
            found = [str(f) for f in found] if isinstance(found, tuple) else [str(found)]
            assert found == [value], (module, flag)


def test_check_refuses_folders_no_module_reads(tmp_path, runs):
    """--check lists, together, a miRNA folder without its isoform folder or --gdc-dir, a salmon
    folder without --rnaseq-dir, and a folder that is not there, each stopping the run. Catches
    ``_pairing`` returning nothing, and a missing folder among the three left unchecked."""
    spec_path = _spec(tmp_path, runs)
    mirna = tmp_path / "mirna"
    options = pipeline.Options("mine", spec_path, FIXTURE, mirna_dir=mirna,
                               salmon_dir=tmp_path)  # fmt: skip
    problems, _ = pipeline.check(tmp_path / "root", options)
    texts = {p.text for p in problems if p.blocks}
    assert {f"--mirna-dir {mirna} is not a folder",
            "--mirna-dir and --isoform-dir go together: gather reads both",
            "--mirna-dir and --isoform-dir take --gdc-dir",
            "--salmon-dir takes --rnaseq-dir: quant reads both"} <= texts, texts  # fmt: skip
    # With --rnaseq-dir but no annotation_gtf, quant is not planned and nothing reads it.
    no_gtf = replace(options, rnaseq_dir=tmp_path, mirna_dir=None)
    texts = {p.text for p in pipeline.check(tmp_path / "root", no_gtf)[0]}
    assert any(t.startswith("--salmon-dir is read by quant alone") for t in texts), texts


def test_run_check_takes_the_folders(tmp_path, runs):
    """`ctx run --check` takes --salmon-dir, --mirna-dir and --isoform-dir and hands them to the
    check. Catches the command dropping one of them before ``pipeline.Options``."""
    spec_path = _spec(tmp_path, runs)
    root = tmp_path / "root"
    root.mkdir()
    done = ctx_run("--root", str(root), "run", "--dataset", "mine", "--dataset-file",
                   str(spec_path), "--fetchngs-dir", str(FIXTURE), "--salmon-dir", str(tmp_path),
                   "--mirna-dir", str(tmp_path / "m"), "--isoform-dir", str(tmp_path / "i"),
                   "--check", cwd=root)  # fmt: skip
    assert done.returncode == 1
    for text in ("--salmon-dir takes --rnaseq-dir", f"--mirna-dir {tmp_path / 'm'} is not",
                 f"--isoform-dir {tmp_path / 'i'} is not"):  # fmt: skip
        assert text in done.stderr, done.stderr


def test_two_projects_reach_gather_each_folder(tmp_path):
    """S.3: a dataset of two GDC projects runs through `ctx run` as `ctx gather` runs it, one
    --gdc-dir per project, and the sheet holds both. Catches every folder but one dropped on
    the way to gather, and the run command taking one --gdc-dir."""
    a, b = folders(tmp_path)
    spec_path = two_projects(tmp_path)
    options = pipeline.Options("fake_brca", spec_path, gdc_dirs=(b, a))
    gather = dict(pipeline.steps(load(spec_path), options))["gather"]
    given = [v for f, v in zip(gather, gather[1:], strict=False) if f == "--gdc-dir"]
    assert given == [str(b), str(a)], gather
    commands = typer.main.get_command(app).commands  # type: ignore[attr-defined]
    parsed = commands["gather"].make_context("gather", list(gather)).params["gdc_dir"]
    assert [str(p) for p in parsed] == [str(b), str(a)], parsed
    root = tmp_path / "root"
    root.mkdir()
    done = ctx_run("--root", str(root), "run", "--dataset", "fake_brca", "--dataset-file",
                   str(spec_path), "--gdc-dir", str(a), "--gdc-dir", str(b), cwd=root)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    assert "ctx run runs: gather, report" in done.stdout, done.stdout
    sheets = (root / "outputs" / "01_gather").glob("*_samplesheet_v1.parquet")
    assert sorted(set(pl.read_parquet(next(sheets))["project"])) == [PROJECT, OTHER]


def test_check_names_every_projects_folder_problem_together(tmp_path):
    """S.3: --check lists, at once, a folder that is not there, one of a project the file does
    not declare, two of one project, and each declared project with no folder, each stopping
    the run. Catches the check stopping at the first folder, a doubled project passing, and a
    project with no folder left unnamed."""
    a, _ = folders(tmp_path)
    twin = build(tmp_path / "twin", [aliquot(4)]).path
    stray = build(tmp_path / "stray", [aliquot(5)], project="FAKE-KIRC").path
    gone = tmp_path / "gone"
    spec_path = two_projects(tmp_path)
    options = pipeline.Options("fake_brca", spec_path, gdc_dirs=(a, twin, stray, gone))
    problems, planned = pipeline.check(tmp_path / "root", options)
    texts = [p.text for p in problems]
    assert all(p.blocks for p in problems) and planned == ["gather", "report"], texts
    assert texts == [
        f"--gdc-dir {gone} is not a folder",
        f"{a.resolve()} and {twin.resolve()} are both {PROJECT!r}",
        f"{stray.resolve()}: project 'FAKE-KIRC', not one of {(PROJECT, OTHER)}",
        f"no --gdc-dir for {OTHER!r}: give one folder per declared project",
    ], texts
    root = tmp_path / "root"
    root.mkdir()
    done = ctx_run("--root", str(root), "run", "--dataset", "fake_brca", "--dataset-file",
                   str(spec_path), "--gdc-dir", str(stray), "--gdc-dir", str(a), "--check",
                   cwd=root)  # fmt: skip
    assert done.returncode == 1
    assert "not one of" in done.stderr and f"no --gdc-dir for {OTHER!r}" in done.stderr


def test_check_refuses_gdc_folders_gather_would_refuse(tmp_path, runs):
    """--gdc-dir with a dataset file that has no [gdc], and a miRNA folder for two projects,
    are problems before any run. Catches either passing --check while gather refuses it."""
    a, b = folders(tmp_path)
    own = _spec(tmp_path, runs)
    found = pipeline.check(tmp_path, pipeline.Options("mine", own, gdc_dirs=(a,)))[0]
    texts = [p.text for p in found]
    assert f"--gdc-dir takes a dataset file with a [gdc] table; {own} has none" in texts, texts
    two = pipeline.Options("fake_brca", two_projects(tmp_path), gdc_dirs=(a, b),
                           mirna_dir=tmp_path, isoform_dir=tmp_path)  # fmt: skip
    texts = [p.text for p in pipeline.check(tmp_path, two)[0]]
    assert texts == ["--mirna-dir reads one project; this dataset file names several"], texts


def _matrix_root(tmp_path: Path) -> tuple[Path, Path, list[Path]]:
    """A count matrix dataset with [de] and [function], and two roots holding its resources."""
    from synth import counts as synthetic  # noqa: PLC0415 - tests/ helper
    from test_gather_counts import DE_SPEC, matrix  # noqa: PLC0415 - its matrix dataset

    made = synthetic(n_genes=400)
    ids = [f"ENSG{i + 1:011d}.1" for i in range(made.counts.height)]
    up = made.truth.with_row_index().filter(pl.col("lfc_condition_type_a") > 0)["index"]
    gmt = [("UP_SET", [f"N{i}" for i in up][:100])]
    gmt += [(f"RANDOM_{j}", [f"N{i}" for i in range(200 + 10 * j, 230 + 10 * j)]) for j in range(4)]
    conf = tmp_path / "conf"
    conf.mkdir()
    made.samples.select("sample", "donor", "cell_type", "condition").write_csv(
        conf / "matde.tsv", separator="\t"
    )
    sets = ('\n[[resources]]\nname = "sets"\nurls = ["http://127.0.0.1/s.gmt"]\nbytes = 0\n'
            'sha256 = "' + "0" * 64 + '"\ndest = "resources/sets.gmt"\n')  # fmt: skip
    function = ('\n[function]\ngene_sets = "sets"\nlicence = "Synthetic; none."\n'
                'featured_set = "UP_SET"\nrunning = ["a_case", "UP_SET"]\n')  # fmt: skip
    head, rest = DE_SPEC.split("[counts]", 1)
    (conf / "matde.toml").write_text(head + sets + "\n[counts]" + rest + function)
    path = matrix(tmp_path / "m.csv", ids, made.counts.drop("gene_id").to_numpy().astype(float),
                  sep=",", samples=made.counts.columns[1:])  # fmt: skip
    gtf = [f'chr1\tt\tgene\t{i + 1}\t{i + 9}\t.\t+\t.\tgene_id "{g}"; gene_type "protein_coding";'
           f' gene_name "N{i}";' for i, g in enumerate(ids)]  # fmt: skip
    roots = [tmp_path / "one", tmp_path / "many"]
    for root in roots:
        (root / "resources").mkdir(parents=True)
        (root / "resources" / "a.gtf.gz").write_bytes(
            gzip.compress(("\n".join(gtf) + "\n").encode())
        )
        (root / "resources" / "sets.gmt").write_text(
            "".join(f"{n}\tsynthetic\t" + "\t".join(g) + "\n" for n, g in gmt)
        )
    return conf / "matde.toml", path, roots


def test_run_matches_the_module_commands_through_de_and_function(tmp_path):
    """650, 738: on a count matrix with [de] and [function], `ctx run` writes every table and
    figure data file the module commands write one by one, gather, qc, de, function and report,
    byte for byte. Catches de or function run by ctx run with options their own command would
    not take, which the gather, qc and report proof above cannot see."""
    spec, path, (one, many) = _matrix_root(tmp_path)
    base = ("--dataset", "matde", "--dataset-file", str(spec))
    done = ctx_run("--root", str(one), "run", *base, "--counts", str(path), cwd=tmp_path)
    assert done.returncode == 0, done.stderr[-3000:]
    printed = re.findall(r"^ctx run: \d+ of \d+: (.+)$", done.stdout, re.M)
    modules = [c.split("--root ")[1].split()[1] for c in printed]
    assert modules == ["gather", "qc", "de", "function", "report"], modules
    for command in (("gather", "--counts", str(path)), ("qc",), ("de",), ("function",),
                    ("report",)):  # fmt: skip
        done = ctx_run("--root", str(many), command[0], *base, *command[1:], cwd=tmp_path)
        assert done.returncode == 0, (command, done.stderr[-3000:])
    tables = {p.relative_to(one): p.read_bytes() for p in one.rglob("*.csv")}
    assert any("06_function" in str(p) for p in tables) and any("05_de" in str(p) for p in tables)
    assert tables == {p.relative_to(many): p.read_bytes() for p in many.rglob("*.csv")}
