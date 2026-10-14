"""`ctx reads`: a student's own FASTQ files laid out as fetchngs.sbatch lays out a download.

The local files are the chr22 fixtures, gzipped at test time. The proof runs the worked example
both ways: the same reads as a fetchngs outdir, written here under the real fixture's headers
and pipeline_info, and as own files through `ctx reads`, give gather tables equal value by
value. No accession is typed here: four_runs.toml and its table name them. Each test names the
mutation it catches.
"""

import csv
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest
import typer
from conftest import ctx_run

from ctx.cli import app
from ctx.dataset import load
from ctx.gather import fetchngs, harmonise, own_fastq
from ctx.gather.curation import read_table
from ctx.report.cards import parameter_cards

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"
REAL = FIXTURES / "fetchngs" / "download"
FOUR = ROOT / "conf" / "datasets" / "four_runs.toml"
TEMPLATE = ROOT / "conf" / "dataset_template.toml"
# What the fixture's study holds by ENA's vocabulary: human RNA-seq. The fetchngs outdir below
# carries them as ENA would, the own table as a student would; the rest of ENA's descriptions
# stay empty on both sides.
DESCRIBED = {"scientific_name": "Homo sapiens", "tax_id": "9606",
             "library_source": "TRANSCRIPTOMIC", "library_strategy": "RNA-Seq"}  # fmt: skip


def _gz(text: str, path: Path) -> Path:
    """The text gzipped with no time stamp, so two writes are the same bytes."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(gzip.compress(text.encode(), mtime=0))
    return path


def _records(name: str) -> list[str]:
    """A fixture FASTQ's records, four lines each."""
    lines = (FIXTURES / name).read_text().splitlines()
    return ["\n".join(lines[i : i + 4]) + "\n" for i in range(0, len(lines), 4)]


def _counted(text: str) -> tuple[int, int]:
    """Reads and bases counted from the text, apart from ctx reads' own counting."""
    lines = text.splitlines()
    return len(lines) // 4, sum(len(s) for s in lines[1::4])


def four_runs_reads() -> dict[str, str]:
    """Each worked-example run's reads: halves of chr22_se's two 101 nt samples, se_c and
    se_d, so the four runs are 101 nt as four_runs.toml expects."""
    halves = []
    for name in ("chr22_se/se_c.fq", "chr22_se/se_d.fq"):
        records = _records(name)
        halves += [records[: len(records) // 2], records[len(records) // 2 :]]
    return {run: "".join(part) for run, part in zip(load(FOUR).accessions, halves, strict=True)}


def fetchngs_outdir(out: Path, reads: dict[str, str], download: bool) -> Path:
    """The reads as fetchngs 1.12.0 lays out a download, under the real fixture's headers.

    The real fixture's versions file and params stand; fetchngs's experiment is the run here,
    since these reads have no ENA experiment, and its FASTQ names follow fetchngs's
    ``<experiment>_<run>``.
    """

    def header(path: Path, delimiter: str) -> list[str]:
        with path.open(newline="") as fh:
            return next(csv.reader(fh, delimiter=delimiter))

    sheet_cols = header(REAL / "samplesheet" / "samplesheet.csv", ",")
    meta_cols = header(next((REAL / "metadata").glob("*.tsv")), "\t")
    map_cols = header(REAL / "samplesheet" / "id_mappings.csv", ",")
    (out / "pipeline_info").mkdir(parents=True)
    shutil.copy(REAL / "pipeline_info" / fetchngs.VERSIONS_FILE, out / "pipeline_info")
    (params_file,) = (REAL / "pipeline_info").glob("params_*.json")
    params = json.loads(params_file.read_text()) | {"skip_fastq_download": not download}
    (out / "pipeline_info" / params_file.name).write_text(json.dumps(params))
    (out / "ids.csv").write_text("".join(f"{run}\n" for run in reads))
    sheet = []
    for run, text in reads.items():
        name = f"{run}_{run}.fastq.gz"
        data = gzip.compress(text.encode(), mtime=0)
        md5 = hashlib.md5(data).hexdigest()
        if download:
            _gz(text, out / "fastq" / name)
            (out / "fastq" / "md5").mkdir(exist_ok=True)
            (out / "fastq" / "md5" / f"{name}.md5").write_text(f"{md5}  {name}\n")
        n, bases = _counted(text)
        row = dict.fromkeys(meta_cols, "") | DESCRIBED | {
            "id": f"{run}_{run}", "run_accession": run, "experiment_accession": run,
            "study_accession": load(FOUR).info["ena_study"], "library_layout": "SINGLE",
            "base_count": str(bases), "read_count": str(n), "fastq_md5": md5,
            "fastq_bytes": str(len(data)), "fastq_ftp": name, "fastq_1": name, "md5_1": md5,
            "single_end": "true",
        }  # fmt: skip
        (out / "metadata").mkdir(exist_ok=True)
        with (out / "metadata" / f"{run}.runinfo_ftp.tsv").open("w", newline="") as fh:
            writer = csv.DictWriter(fh, meta_cols, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerow(row)
        sheet.append(row | {"sample": run, "fastq_1": f"out/fastq/{name}", "fastq_2": "",
                            "strandedness": "auto"})  # fmt: skip
    for name, cols in (("samplesheet.csv", sheet_cols), ("id_mappings.csv", map_cols)):
        (out / "samplesheet").mkdir(exist_ok=True)
        with (out / "samplesheet" / name).open("w", newline="") as fh:
            writer = csv.DictWriter(fh, cols, quoting=csv.QUOTE_ALL, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(sheet)
    return out


def own_files(folder: Path, reads: dict[str, str]) -> Path:
    """The reads as a student holds them: one gzipped file per run, named by the run."""
    for run, text in reads.items():
        _gz(text, folder / f"{run}.fastq.gz")
    return folder


def own_table(path: Path) -> Path:
    """four_runs' own table with ENA's descriptions a student would record, as a CSV."""
    rows, _ = own_fastq.table_rows(load(FOUR).table())
    described = {"study_accession": load(FOUR).info["ena_study"], **DESCRIBED}
    with path.open("w", newline="") as fh:
        csv.writer(fh).writerows([rows[0] + list(described),
                                  *(r + list(described.values()) for r in rows[1:])])  # fmt: skip
    return path


def reads_cli(tmp_path: Path, run_folder: Path, folder: Path, table: Path, column: str):
    """`ctx reads` as a student types it."""
    return ctx_run("reads", "--run-folder", str(run_folder), "--fastq-dir", str(folder),
                   "--table", str(table), "--sample-column", column, cwd=tmp_path)  # fmt: skip


@pytest.fixture(scope="module")
def both_ways(tmp_path_factory):
    """The worked example's reads both ways, and gather run on each with four_runs.toml."""
    tmp = tmp_path_factory.mktemp("both")
    reads = four_runs_reads()
    done = reads_cli(tmp, tmp / "run", own_files(tmp / "mine", reads),
                     own_table(tmp / "mine.csv"), "run_accession")  # fmt: skip
    assert done.returncode == 0, done.stderr
    outdirs = {"fetchngs": fetchngs_outdir(tmp / "fetchngs" / "results", reads, True),
               "own": tmp / "run" / own_fastq.DOWNLOAD / "results"}  # fmt: skip
    roots = {}
    for route, outdir in outdirs.items():
        roots[route] = tmp / f"root_{route}"
        done = ctx_run("--root", str(roots[route]), "gather", "--dataset", "four_runs",
                       "--dataset-file", str(FOUR), "--fetchngs-dir", str(outdir),
                       cwd=tmp)  # fmt: skip
        assert done.returncode == 0, (route, done.stderr[-3000:])
    return reads, outdirs, roots, tmp


def _tables(root: Path) -> dict[str, pl.DataFrame]:
    return {p.name: pl.read_csv(p, infer_schema=False)
            for p in sorted((root / "outputs").rglob("*.csv"))}  # fmt: skip


def _fastq_md5(outdir: Path, frame: pl.DataFrame) -> pl.DataFrame:
    """Each FASTQ path column replaced by the MD5 of the file it names."""

    def digest(value: str | None) -> str:
        return hashlib.md5((outdir / value).read_bytes()).hexdigest() if value else ""

    return frame.with_columns(
        [pl.col(c).map_elements(digest, return_dtype=pl.String) for c in ("fastq_1", "fastq_2")]
    )


def test_both_ways_give_equal_gather_tables(both_ways):
    """The proof: the worked example's reads as a fetchngs download and as own files give gather
    tables equal value by value: the sheet, the composition and the library sizes. The FASTQ
    columns name files, each route's own name, so they are compared by the bytes they name.
    Catches reads or bases counted wrongly (a newline counted, a header line as sequence), the
    MD5 taken of the unzipped text, the layout or a table description lost, run_index left
    empty, and curation joined on another column."""
    _, outdirs, roots, _ = both_ways
    fetched, own = _tables(roots["fetchngs"]), _tables(roots["own"])
    assert len(fetched) >= 3 and set(fetched) == set(own), (set(fetched), set(own))
    sheet = next(n for n in fetched if "_samplesheet_" in n)
    for name in fetched:
        a, b = fetched[name], own[name]
        if name == sheet:
            a, b = _fastq_md5(outdirs["fetchngs"], a), _fastq_md5(outdirs["own"], b)
        assert a.equals(b), (name, a, b)
    frame = own[sheet]
    reads = four_runs_reads()
    assert frame["read_count"].to_list() == [str(_counted(reads[r])[0])
                                             for r in frame["run_accession"]]  # fmt: skip
    assert set(frame["status_1"]) == {"verified"} and set(frame["read_length"]) == {"101.0"}
    assert set(frame["library_layout"]) == {"SINGLE"} and frame["run_index"].to_list() == ["1"] * 4


def test_both_ways_metadata_outdirs_curate_alike(both_ways):
    """Part 3 reads the metadata run's outdir: ctx reads writes one, and it curates as the
    fetchngs one does. Catches the metadata outdir left out, or written unlike the download's."""
    reads, _, _, tmp = both_ways
    meta = fetchngs_outdir(tmp / "fetchngs_meta" / "results", reads, False)
    spec = load(FOUR)
    sheets = []
    for outdir in (meta, tmp / "run" / own_fastq.METADATA / "results"):
        checked = fetchngs.read(outdir)
        sheets.append(harmonise.sheet(checked, harmonise.curated_runs(tmp, spec, checked, {})))
    assert sheets[0].equals(sheets[1])
    assert set(sheets[0]["status_1"]) == {"not_downloaded"}


def test_headers_are_fetchngs():
    """The samplesheet, metadata and mapping headers are fetchngs 1.12.0's, from its real
    output. Catches a column dropped, added or moved in ctx reads' copy."""
    with (REAL / "samplesheet" / "samplesheet.csv").open(newline="") as fh:
        assert tuple(next(csv.reader(fh))) == own_fastq.SHEET_COLUMNS
    with (REAL / "samplesheet" / "id_mappings.csv").open(newline="") as fh:
        assert tuple(next(csv.reader(fh))) == own_fastq.MAPPING_COLUMNS
    for path in (REAL / "metadata").glob("*.tsv"):
        assert tuple(path.read_text().splitlines()[0].split("\t")) == own_fastq.METADATA_COLUMNS


@pytest.fixture
def lanes(tmp_path):
    """chr22's pairs over two lanes of P1, Illumina names, the second lane in a subfolder, and
    chr22_se's se_a single-end; a table one row per sample."""
    fq = tmp_path / "fq"
    for lane, sub in (("L001", ""), ("L002", "sub")):
        for mate in ("1", "2"):
            text = (FIXTURES / "chr22" / f"reads_R{mate}.fq").read_text()
            _gz(text, fq / sub / f"P1_S1_{lane}_R{mate}_001.fastq.gz")
    _gz((FIXTURES / "chr22_se" / "se_a.fq").read_text(), fq / "se_a.fq.gz")
    table = tmp_path / "mine.csv"
    table.write_text("sample,donor,cell_type,condition\nP1,d1,monocyte,control\n"
                     "se_a,d2,monocyte,case\n")  # fmt: skip
    return fq, table


def test_lanes_pairs_and_single_end(tmp_path, lanes):
    """Two lanes of a paired sample are two runs of it, mates matched by name, beside a
    single-end sample; counts are the fixtures' own: chr22's 2,500 pairs of 100 nt, chr22_se's
    se_a 1,000 reads of 76 nt (its samples.tsv). The filled template, joined on sample, curates
    all three runs. Catches mates paired across lanes, the lane taken for the sample, a pair
    counted as two reads, run_index empty for a name with no number, and curation that drops
    run_accession when it joins on sample."""
    fq, table = lanes
    run, sheet = own_fastq.write(tmp_path / "run", fq, table, "sample")
    runs = fetchngs.read(run).runs.sort("run_accession")
    assert runs["run_accession"].to_list() == ["P1_S1_L001_001", "P1_S1_L002_001", "se_a"]
    assert runs["sample"].to_list() == ["P1", "P1", "se_a"]
    assert runs["n_runs"].to_list() == [2, 2, 1] and runs["run_index"].to_list() == [1, 2, 1]
    pairs, bases = _counted((FIXTURES / "chr22" / "reads_R1.fq").read_text())
    _, bases_2 = _counted((FIXTURES / "chr22" / "reads_R2.fq").read_text())
    assert runs["read_count"].to_list() == [str(pairs), str(pairs), "1000"]
    assert runs["base_count"].to_list() == [str(bases + bases_2)] * 2 + [str(1000 * 76)]
    assert runs["single_end"].to_list() == [False, False, True]
    assert runs["fastq_2"].to_list() == ["fastq/P1_S1_L001_R2_001.fastq.gz",
                                          "fastq/P1_S1_L002_R2_001.fastq.gz", ""]  # fmt: skip
    spec_path = tmp_path / "mine.toml"
    spec_path.write_text(TEMPLATE.read_text().replace(
        'join = ["run_accession", "run_accession"]', 'join = ["sample", "sample"]'))  # fmt: skip
    curated = harmonise.curated_runs(tmp_path / "root", load(spec_path), fetchngs.read(run), {})
    assert sorted(curated["run_accession"]) == sorted(runs["run_accession"])
    assert curated.sort("run_accession")["condition"].to_list() == ["control", "control", "case"]
    assert sheet == tmp_path / "run" / own_fastq.METADATA / "results"


def test_files_stay_where_they_are(tmp_path, lanes):
    """The FASTQ files are linked, never moved or copied, and the samplesheet nf-core/rnaseq
    reads names each where it is, by its full path. Catches a copy into fastq/, a move, and a
    samplesheet naming the link."""
    fq, table = lanes
    before = {p: (p.stat().st_ino, p.read_bytes()) for p in fq.rglob("*.gz")}
    run, _ = own_fastq.write(tmp_path / "run", fq, table, "sample")
    assert {p: (p.stat().st_ino, p.read_bytes()) for p in fq.rglob("*.gz")} == before
    links = list((run / "fastq").glob("*.gz"))
    assert len(links) == len(before) and all(p.is_symlink() for p in links)
    assert {p.resolve() for p in links} == {p.resolve() for p in before}
    with (run / "samplesheet" / "samplesheet.csv").open(newline="") as fh:
        named = {r[k] for r in csv.DictReader(fh) for k in ("fastq_1", "fastq_2") if r[k]}
    assert named == {str(p.resolve()) for p in before}


def test_every_problem_listed_at_once(tmp_path):
    """One call lists every problem of the names and the table, and writes nothing. Catches a
    check that stops at its first problem, and a run folder written before the checks pass."""
    fq = tmp_path / "fq"
    text = (FIXTURES / "chr22_se" / "se_a.fq").read_text()
    for name in (
        "A_R1.fastq.gz",
        "B_R1.fastq.gz",
        "B_1.fastq.gz",
        "B_R2.fastq.gz",
        "stray.fastq.gz",
        "C_x.fq.gz",
        "D.fastq.gz",
        "D_R1.fastq.gz",
    ):
        _gz(text, fq / name)
    (fq / "E.fastq").write_text(text)
    table = tmp_path / "mine.csv"
    table.write_text("sample,read_count\nA,1\nB,1\nC,1\nC_x,1\nD,1\nF,1\nA,1\n")
    with pytest.raises(own_fastq.ReadsError) as caught:
        own_fastq.write(tmp_path / "run", fq, table, "sample")
    found = "\n".join(caught.value.problems)
    for said in ("A_R1.fastq.gz: mate 1 of run A has no mate 2",
                 "B_R1.fastq.gz and B_1.fastq.gz: both are run B mate 1",
                 "stray.fastq.gz: names no sample", "C_x.fq.gz: names C and C_x",
                 "run D is both single-end and paired", "E.fastq: not gzipped",
                 "sample F: no FASTQ file starts with its name",
                 "line 8: sample A is in an earlier row too",
                 "column 'read_count' is measured from the files"):  # fmt: skip
        assert said in found, (said, found)
    assert not (tmp_path / "run").exists()


def test_unreadable_files_refused_before_writing(tmp_path):
    """A truncated record and mates of unequal length are refused together, after every file is
    read, and nothing is written. Catches a record cut short accepted, unequal mates accepted,
    and outdirs left behind by a failed read."""
    fq = tmp_path / "fq"
    records = _records("chr22/reads_R1.fq")
    _gz("".join(records[:10]), fq / "P_R1.fastq.gz")
    _gz("".join(_records("chr22/reads_R2.fq")[:9]), fq / "P_R2.fastq.gz")
    _gz("".join(records[:3]) + records[3].rsplit("\n", 2)[0] + "\n", fq / "Q.fastq.gz")
    table = tmp_path / "mine.tsv"
    table.write_text("sample\nP\nQ\n")
    (tmp_path / "run").mkdir()
    with pytest.raises(own_fastq.ReadsError) as caught:
        own_fastq.write(tmp_path / "run", fq, table, "sample")
    found = "\n".join(caught.value.problems)
    assert "Q.fastq.gz record 4: not a whole FASTQ record" in found
    assert "run P: mates hold P_R1.fastq.gz 10 and P_R2.fastq.gz 9 reads" in found
    assert list((tmp_path / "run").iterdir()) == []


def test_written_outdir_is_not_overwritten(tmp_path, lanes):
    """A second run into the same folder is refused and the first outdirs stand. Catches an
    outdir overwritten in place."""
    fq, table = lanes
    run, _ = own_fastq.write(tmp_path / "run", fq, table, "sample")
    first = (run / "samplesheet" / "samplesheet.csv").read_bytes()
    with pytest.raises(own_fastq.ReadsError, match="exists; move it aside"):
        own_fastq.write(tmp_path / "run", fq, table, "sample")
    assert (run / "samplesheet" / "samplesheet.csv").read_bytes() == first
    assert sorted(p.name for p in (tmp_path / "run").iterdir()) == [own_fastq.DOWNLOAD,
                                                                   own_fastq.METADATA]  # fmt: skip


def test_script_runs_on_the_standard_library(tmp_path, lanes):
    """On Kebnekaise the file runs on the fetch environment's Python, where ctx is not installed:
    here with no site-packages and no PYTHONPATH, it writes what `ctx reads` writes. Catches a
    third-party or ctx import in own_fastq."""
    fq, table = lanes
    args = ["--fastq-dir", str(fq), "--table", str(table), "--sample-column", "sample"]
    script = ROOT / "src" / "ctx" / "gather" / "own_fastq.py"
    done = subprocess.run([sys.executable, "-I", "-S", str(script), "--run-folder",
                           str(tmp_path / "a"), *args], capture_output=True, text=True,
                          check=False, cwd=tmp_path)  # fmt: skip
    assert done.returncode == 0, done.stderr
    pairs, _ = _counted((FIXTURES / "chr22" / "reads_R1.fq").read_text())
    said = f"2 samples; 3 runs, 2 paired; {2 * pairs + 1000:,} reads or pairs"
    assert done.stdout.splitlines()[-1] == said
    cli = ctx_run("reads", "--run-folder", str(tmp_path / "b"), *args, cwd=tmp_path)
    assert cli.returncode == 0, cli.stderr
    sheets = [(tmp_path / f / own_fastq.DOWNLOAD / "results" / "samplesheet" / "samplesheet.csv")
              .read_bytes() for f in ("a", "b")]  # fmt: skip
    assert sheets[0] == sheets[1]


def test_both_entry_points_take_the_same_options():
    """`ctx reads` and the file run as a script take the same options with the same defaults,
    and each has a parameter card. Catches an option renamed or a default changed in one."""
    group = typer.main.get_command(app)
    command = group.commands["reads"]
    typed = {p.opts[0]: p.default for p in command.params}
    parsed = {a.option_strings[0]: a.default for a in own_fastq.parser()._actions
              if a.option_strings and a.option_strings[0] != "-h"}  # fmt: skip
    assert typed == parsed
    cards = {c.name: c.default for c in parameter_cards(command)}
    assert cards == {k: "required" if v is None else str(v) for k, v in typed.items()}


def test_gather_refuses_another_layout(tmp_path, lanes):
    """gather reads ctx reads' layout version alone, as it reads fetchngs 1.12.0 alone. Catches
    the layout version left unchecked."""
    fq, table = lanes
    run, _ = own_fastq.write(tmp_path / "run", fq, table, "sample")
    versions = run / "pipeline_info" / own_fastq.VERSIONS_FILE
    versions.write_text(versions.read_text().replace(f"layout {own_fastq.LAYOUT}", "layout 0"))
    with pytest.raises(fetchngs.FetchngsError, match="ctx reads layout 0"):
        fetchngs.read(run)


def test_own_table_read_one_way():
    """curation reads an own table through ctx reads' reader, so both read it alike. Catches a
    second copy of the reader drifting."""
    path = load(FOUR).table()
    body, digest = read_table(path)
    rows, same = own_fastq.table_rows(path)
    assert digest == same and body.splitlines() == ["\t".join(r) for r in rows]


def test_run_takes_the_route(tmp_path, lanes):
    """`ctx run --check`, then `ctx run`, take ctx reads' outdir as --fetchngs-dir and run gather,
    qc and report on it. Catches the route refused by run's own checks."""
    fq, table = lanes
    run, _ = own_fastq.write(tmp_path / "run", fq, table, "sample")
    spec_path = tmp_path / "mine.toml"
    shutil.copy(table, tmp_path / "own.csv")
    text = TEMPLATE.read_text().replace('join = ["run_accession", "run_accession"]',
                                        'join = ["sample", "sample"]')  # fmt: skip
    # Without resources, so nothing needs fetching and nothing stops the run.
    text = re.sub(r"^\[\[resources\]\]\n(?:.+\n)+", "", text, flags=re.M)
    spec_path.write_text(text.replace('file = "mine.csv"', 'file = "own.csv"'))
    root = tmp_path / "root"
    args = [
        "--root",
        str(root),
        "run",
        "--dataset",
        "mine",
        "--dataset-file",
        str(spec_path),
        "--fetchngs-dir",
        str(run),
    ]
    check = ctx_run(*args, "--check", cwd=tmp_path)
    # The one problem left is the contrast it suggests, which leaves de out and stops nothing.
    assert "1 problems, 0 of them stopping the run; ctx run runs: gather, qc, report" in (
        check.stdout
    ), check.stdout + check.stderr
    assert "no [[contrasts]]" in check.stderr and "curation" not in check.stderr
    done = ctx_run(*args, cwd=tmp_path)
    assert done.returncode == 0, done.stderr[-3000:]
    (sheet,) = (root / "outputs" / "01_gather").glob("*_samplesheet_v1.csv")
    assert pl.read_csv(sheet)["donor"].to_list() == ["d1", "d1", "d2"]
