"""ED.6e-1 to -3: M2's read stage from a cluster's report folder, against the local run.

The folder is built here as ED.6's block C writes it and block D copies it, from the chr22_se
reads and fetchngs's paired run: the tools run per sample, then everything block D leaves on
the cluster is deleted, and the FASTQ files go too, so only the reports can be read.
"""

import gzip
import html
import shutil
from dataclasses import replace
from datetime import date
from pathlib import Path

import polars as pl
import pytest
from conftest import once, writable_copy
from test_qc_run import GRID_SAMPLE, SPEC, ctx, with_reads  # noqa: F401 - session fixture

from ctx import tools
from ctx.de.run import Writer
from ctx.qc import Fastp, QcError, cluster, multiqc_io, reads, registry, run
from ctx.qc.blocks import Reads
from ctx.report.builder import unpacked

ROOT = Path(__file__).resolve().parents[1]
SE = ROOT / "tests" / "fixtures" / "chr22_se"
PAIRED = ROOT / "tests" / "fixtures" / "fetchngs" / "download" / "fastq"
PAIR = "SRX9626017_SRR13191702"
REG = registry.load(ROOT / "conf" / "qc_metrics.toml")
# Task samples as fetchngs names them (experiment accessions), unlike the sheet's; one holds
# two runs, as BRCA-immune-cells-monocytes's two-run samples do.
TASKS = (("EX_ab", (("R_se_a", ["se_a.fastq.gz"]), ("R_se_b", ["se_b.fastq.gz"]))),
         ("EX_c", (("R_se_c", ["se_c.fastq.gz"]),)), ("EX_d", (("R_se_d", ["se_d.fastq.gz"]),)),
         ("EX_e", (("R_se_e", ["se_e.fastq.gz"]),)),
         ("EX_pair", (("SRR13191702", [f"{PAIR}_1.fastq.gz", f"{PAIR}_2.fastq.gz"]),)))  # fmt: skip
# What block D copies of a run folder, for the read stage.
KEPT = ("tasks.tsv", "*.fastp.json", "*_fastqc.zip", "*.seqkit_stats.tsv")


def build_align(fastq: Path, tasks, dest: Path, grid_run: str) -> Path:
    """An ED.6 report folder: block C's calls per sample, the grid on one run, then pruned."""
    lines = ["sample\tmates\truns\tr1\tr2\tbytes"]
    for sample, runs in tasks:
        folder = dest / "runs" / sample
        for sub in ("fastp", "fastqc", "seqkit"):
            (folder / sub).mkdir(parents=True)
        mates = len(runs[0][1])
        raw = [str(fastq / names[m]) for m in range(mates) for _, names in runs]
        for run_name, names in runs:
            reads.run_fastp([str(fastq / n) for n in names], folder / "fastp", run_name, Fastp(),
                            2, trimmed=folder / "trimmed")  # fmt: skip
        tools.run("fastqc", ["--quiet", "--threads", "2", "--outdir", str(folder / "fastqc"), *raw],
                  cwd=folder / "fastqc", log=folder / "fastqc" / "fastqc.log")  # fmt: skip
        tools.run("seqkit", ["stats", "-a", "-T", "-j", "2", "-o",
                             str(folder / "seqkit" / f"{sample}.seqkit_stats.tsv"), *raw],
                  cwd=folder / "seqkit", log=folder / "seqkit" / "seqkit.log")  # fmt: skip
        r1 = ",".join(n[0] for _, n in runs)
        r2 = ",".join(n[1] for _, n in runs) if mates == 2 else "-"
        lines.append(f"{sample}\t{mates}\t{','.join(r for r, _ in runs)}\t{r1}\t{r2}\t1")
    names = next(n for _, runs in tasks for r, n in runs if r == grid_run)
    (dest / "grid").mkdir()
    for title, (q, n) in reads.grid_cells((reads.GRID_Q, reads.GRID_L)).items():
        reads.run_fastp([str(fastq / f) for f in names], dest / "grid", title,
                        replace(Fastp(), q=q, l=n), 2)  # fmt: skip
    (dest / "tasks.tsv").write_text("\n".join(lines) + "\n")
    for path in sorted(dest.rglob("*"), reverse=True):
        if path.is_file() and not any(path.match(k) for k in KEPT):
            path.unlink()
        elif path.is_dir() and not any(path.iterdir()):
            path.rmdir()
    return dest


def sheet() -> pl.DataFrame:
    """The gather sheet of a metadata-only outdir: no FASTQ named, every run not downloaded."""
    runs = [(r, f"S_{r}", "PAIRED" if len(n) == 2 else "SINGLE") for _, rs in TASKS for r, n in rs]
    return pl.DataFrame(runs, orient="row", schema=["run_accession", "sample", "library_layout"]
                        ).with_columns(fastq_1=pl.lit(""), status_1=pl.lit("not_downloaded"),
                                       fastq_2=pl.lit(""), status_2=pl.lit(""))  # fmt: skip


def _built(folder: Path):
    fastq = folder / "fastq"
    fastq.mkdir()
    for name in "abcde":
        with (SE / f"se_{name}.fq").open("rb") as src, gzip.open(fastq / f"se_{name}.fastq.gz",
                                                                   "wb") as dst:  # fmt: skip
            shutil.copyfileobj(src, dst)
    for mate in (1, 2):
        shutil.copy(PAIRED / f"{PAIR}_{mate}.fastq.gz", fastq)
    align = build_align(fastq, TASKS, folder / "align", "R_se_a")
    # The local read stage on the same files: the values the cluster's reports must give.
    files = pl.DataFrame([(r, f"S_{r}", m + 1, str(fastq / n)) for _, rs in TASKS
                          for r, ns in rs for m, n in enumerate(ns)],
                         schema=reads.FILE_COLUMNS, orient="row")  # fmt: skip
    folders = reads.run_tools(files, folder / "local" / "logs", threads=2, fastp=Fastp())
    config = reads.raw_config(folder / "local" / "config.json")
    mqc = multiqc_io.run(folders, folder / "local", config=config, modules=reads.MODULES)
    values = reads.run_values(mqc, files, REG)
    grid = reads.fastp_grid([str(fastq / "se_a.fastq.gz")], folder / "local" / "grid_logs",
                            folder / "local" / "grid", grid=(reads.GRID_Q, reads.GRID_L),
                            fastp=Fastp(), threads=2)  # fmt: skip
    shutil.rmtree(fastq)  # the original made unreachable: only the reports remain
    return align, values, grid


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """(the report folder, the local run's values per run, the local grid on R_se_a)."""
    return once(tmp_path_factory, "qc_cluster_align", _built)


def cluster_reads(align: Path, out: Path, the_sheet: pl.DataFrame, **given) -> Reads:
    writer = Writer(out / "02_qc", "test", date(2026, 10, 22), {})
    options = run.Options(SPEC, reads=True, align_dir=align, **given)
    return run.cluster_reads(writer, the_sheet, REG, options, align)


def test_the_reports_give_the_local_values_run_by_run(built, tmp_path):
    """ED.6e-1, -3: every run mapped by run and stem, seqkit read from one table per sample.

    Mutation caught: the reads stage reading the first sample's seqkit table alone."""
    align, local, _ = built
    seen = cluster_reads(align, tmp_path, sheet())
    order = local["run_accession"].to_list()
    got = seen.values.with_columns(sample=pl.lit("")).sort("run_accession")
    want = local.with_columns(sample=pl.lit("")).sort("run_accession")
    assert got.equals(want), (got, want)
    assert seen.values["run_accession"].to_list() == order
    assert seen.values["seqkit_seqkit-Q30_pct"].null_count() == 0
    assert len(list(align.glob("runs/*/seqkit/*.seqkit_stats.tsv"))) == len(TASKS)
    assert seen.skipped.height == 0 and seen.source


def test_the_grid_is_the_clusters_and_equals_the_local_one(built, tmp_path):
    """ED.6e-2: the grid's MultiQC half reads grid/*.fastp.json; the default run is the
    sheet's first, so a grid run on any other is refused.

    Mutations caught: the default grid run taken from the sheet's second run; the grid's
    input files left unchecked."""
    align, _, local_grid = built
    seen = cluster_reads(align, tmp_path / "a", sheet())
    assert seen.grid_run == "R_se_a" and seen.grid.equals(local_grid)
    with pytest.raises(QcError, match=r"grid_q10_l15: ran \['se_a.fastq.gz'\].*'se_c.fastq.gz'"):
        cluster_reads(align, tmp_path / "b", sheet(), grid_run="R_se_c")


def test_a_task_list_unlike_the_sheet_is_refused_naming_the_runs(built, tmp_path):
    align, _, _ = built
    other = sheet().filter(pl.col("run_accession") != "R_se_e").vstack(
        sheet().head(1).with_columns(run_accession=pl.lit("R_new")))  # fmt: skip
    with pytest.raises(QcError, match=r"only in the sheet \['R_new'\]; only in the task list "
                                      r"\['R_se_e'\]"):  # fmt: skip
        cluster_reads(align, tmp_path, other)
    swapped = sheet().with_columns(library_layout=pl.lit("PAIRED"))
    with pytest.raises(QcError, match=r"not the sheet's library_layout: \['R_se_a', 'R_se_b'"):
        cluster_reads(align, tmp_path, swapped)


def test_a_missing_report_or_another_filter_is_refused(built, tmp_path):
    align = writable_copy(built[0], tmp_path / "align")
    with pytest.raises(QcError, match=r"fastp reports not run as this ctx qc asks: R_se_a: "
                                      r"ran \['se_a.fastq.gz'\] with \['-q', '15'"):  # fmt: skip
        cluster_reads(align, tmp_path / "a", sheet(), fastp=Fastp(q=20))
    (align / "runs" / "EX_ab" / "fastqc" / "se_b_fastqc.zip").unlink()
    (align / "runs" / "EX_pair" / "fastp" / "SRR13191702.fastp.json").unlink()
    with pytest.raises(QcError, match=r"2 reports absent .*EX_ab/fastqc/se_b_fastqc.zip.*"
                                      r"EX_pair/fastp/SRR13191702.fastp.json"):  # fmt: skip
        cluster_reads(align, tmp_path / "b", sheet())


def test_a_corrupt_report_multiqc_cannot_read_is_refused(built, tmp_path):
    """A zip present but unreadable: MultiQC skips it, and ctx must stop rather than draw."""
    align = writable_copy(built[0], tmp_path / "align")
    (align / "runs" / "EX_c" / "fastqc" / "se_c_fastqc.zip").write_bytes(b"not a zip")
    with pytest.raises(QcError, match=r"no fastqc report for \['se_c'\]"):
        cluster_reads(align, tmp_path / "a", sheet())


def test_the_command_split_keeps_the_filter_alone():
    command = ("/x/fastp -i /fq/a_1.fastq.gz -I /fq/a_2.fastq.gz -R a -j a.fastp.json -h a.html "
               "-w 8 -q 15 -u 40 -l 15 -o t_1.fastq.gz -O t_2.fastq.gz ")  # fmt: skip
    named, rest = cluster.split_command(command)
    assert rest == Fastp().args()
    assert [Path(named[f]).name for f in ("-i", "-I")] == ["a_1.fastq.gz", "a_2.fastq.gz"]


@pytest.fixture(scope="session")
def from_cluster(with_reads, tmp_path_factory):  # noqa: F811 - the imported fixture
    """``ctx qc --reads --align-dir`` on with_reads's gather outputs, its FASTQs gone."""
    return once(tmp_path_factory, "qc_cluster_root", lambda f: _from_cluster(with_reads, f))


def _from_cluster(local, folder: Path):
    out, runs = local
    root = writable_copy(out.parents[1], folder / "root", ignore=shutil.ignore_patterns("02_qc"))
    raw = root / "data" / "raw"
    tasks = [(f"EX_{s}", ((r, [f"{s}.fastq.gz"]),)) for s, r in runs.items()]
    build_align(raw / "desynth" / "fastq", tasks, raw / "desynth_align", runs[GRID_SAMPLE])
    shutil.rmtree(raw / "desynth")
    sheet_path = next((root / "outputs" / "01_gather").glob("*_samplesheet_v1.parquet"))
    gone = {"fastq_1": pl.lit(""), "status_1": pl.lit("not_downloaded")}
    pl.read_parquet(sheet_path).with_columns(**gone).write_parquet(sheet_path)
    done = ctx("--root", str(root), "qc", "--dataset", "desynth", "--dataset-file", str(SPEC),
               "--reads", "--grid-run", runs[GRID_SAMPLE], "--align-dir",
               str(raw / "desynth_align"), cwd=root)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    refused = ctx("--root", str(root), "qc", "--dataset", "desynth", "--dataset-file", str(SPEC),
                  "--align-dir", str(raw / "desynth_align"), cwd=root)  # fmt: skip
    return root / "outputs" / "02_qc", refused.returncode, refused.stderr


def test_ctx_qc_from_the_cluster_writes_the_local_tables(with_reads, from_cluster):  # noqa: F811
    """End to end: the read tables equal the local run's, and the page says where they ran."""
    out, code, stderr = from_cluster
    for what in ("read_metrics", "read_flags", "grid", "fastp_grid", "status", "skipped"):
        got = pl.read_csv(next(out.glob(f"*_desynth_{what}_v1.csv")))
        want = pl.read_csv(next(with_reads[0].glob(f"*_desynth_{what}_v1.csv")))
        assert got.equals(want), what
    page = html.unescape(unpacked(next(out.glob("*_overview_v1.html")).read_text()))
    note = ("FastQC, fastp and seqkit ran on the cluster; MultiQC read their reports for 5 runs "
            "of 5 samples from desynth_align.")  # fmt: skip
    assert note in page and "tasks.tsv" in page
    assert code == 1 and "--align-dir is read with --reads" in stderr
