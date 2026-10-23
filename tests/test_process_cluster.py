"""`ctx process --align-dir`: M3 from a cluster's alignment (TODO.md ED.6e-4 to ED.6e-8).

The cluster folder is tests/process_cluster.py's, laid out as ED.6 block D copies it from the
local HISAT2 run of ``conftest.process_root``, so the two routes read the same alignments and
their tables must agree. The reads are deleted from the cluster root before the run, so a
route that opened a FASTQ would fail. Each test names the mutation it catches.
"""

import gzip
import hashlib
import shutil
from pathlib import Path

import polars as pl
import pytest
from conftest import ctx_run, once, writable_copy
from process_root import DATASET, plant

from ctx import pipeline
from ctx.dataset import load
from ctx.process import ProcessError, cluster, count, fc_compact
from ctx.process import run as process_run
from ctx.process.stages import Options
from ctx.report.builder import unpacked

ROOT = Path(__file__).resolve().parents[1]
OUT = ("outputs", "03_process")


def dataset_file(root: Path) -> Path:
    """A dataset file whose annotation is the slice's own GTF, as the cluster counted with."""
    import process_cluster  # noqa: PLC0415 - the fixture module, after conftest's checks

    gtf = root / "resources" / "gencode" / process_cluster.GTF.name
    gtf.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(process_cluster.GTF, gtf)
    digest = hashlib.sha256(gtf.read_bytes()).hexdigest()
    spec = root / f"{DATASET}.toml"
    spec.write_text(f'[dataset]\nslug = "{DATASET}"\n\n[guard]\nids = []\npatterns = []\n\n'
                    '[[resources]]\nname = "annotation_gtf"\n'
                    f'urls = ["http://127.0.0.1/{gtf.name}"]\nbytes = {gtf.stat().st_size}\n'
                    f'sha256 = "{digest}"\ndest = "{gtf.relative_to(root)}"\n')  # fmt: skip
    return spec


@pytest.fixture(scope="session")
def cluster_run(process_root, tmp_path_factory):
    """(root, dataset file, align folder, completed run) of `ctx process --align-dir`."""

    def build(folder: Path):
        import process_cluster  # noqa: PLC0415

        root = folder / "root"
        root.mkdir()
        plant(root)
        # The reads must be unreachable: this route reads the cluster's files, never a FASTQ.
        shutil.rmtree(root / "data" / "raw" / DATASET / "fastq")
        spec = dataset_file(root)
        align = process_cluster.build(process_root[0], root, DATASET)
        done = ctx_run("--root", str(root), "process", "--dataset", DATASET, "--dataset-file",
                       str(spec), "--align-dir", str(align), "--settings", "cut,untrimmed",
                       cwd=root)  # fmt: skip
        return root, spec, align, done

    return once(tmp_path_factory, "process_cluster", build)


def newest(root: Path, what: str) -> pl.DataFrame:
    """A table ctx process wrote, read from its Parquet file."""
    return pl.read_parquet(next(root.joinpath(*OUT).glob(f"*_{DATASET}_{what}_v1.parquet")))


def test_the_cluster_route_writes_the_page_from_no_fastq(cluster_run):
    """ED.6e-4: with every FASTQ deleted, --align-dir writes the page and its tables from the
    copied files, no index built and the cluster's index named. Catches a route that opens a
    FASTQ, builds an index here, or names another index."""
    root, _, _, done = cluster_run
    assert done.returncode == 0, done.stderr[-3000:]
    folders = dict(newest(root, "folders").select("step", "folder").rows())
    assert folders["index_main"] == "star_2.7.11b_gencode50_sjdb100"
    assert folders["index_cut"] == "star_2.7.11b_gencode50_half"
    assert (
        not newest(root, "folders").filter(pl.col("step").str.starts_with("index"))["built"].any()
    )
    assert {"align_main", "align_cut", "align_untrimmed", "count_plain", "count_variants",
            "count_cut"} <= set(folders)  # fmt: skip
    assert not list(root.joinpath(*OUT, "tools", "index").glob("star_*"))


def test_both_routes_give_the_same_tables(cluster_run, process_root):
    """ED.6e-4, ED.6e-5: the cluster's per-sample files of the local BAMs, joined, give the local
    run's gene counts at the called strand, its strand calls, every counter's Assigned and the
    featured sample's half-reference genes. Catches a join that drops or reorders a sample, a
    variant taken at another strand, and the cut counted on the whole annotation."""
    root, local = cluster_run[0], process_root[0]
    for what, key in (("gene_counts", "gene_id"), ("strandedness", "sample"),
                      ("counting", ["sample", "counter"])):  # fmt: skip
        mine, theirs = newest(root, what), newest(local, what)
        columns = [c for c in mine.columns if c in theirs.columns and not c.startswith(
            ("compare_", "agrees", "star_"))]  # fmt: skip
        assert mine.select(columns).sort(key).equals(theirs.select(columns).sort(key)), what
    incomplete = ("genes_kept", "genes_removed", "kept_reads_main", "kept_reads_cut")
    assert newest(root, "incomplete").select(incomplete).equals(
        newest(local, "incomplete").select(incomplete))  # fmt: skip


def test_the_variants_are_read_at_each_samples_called_strand(cluster_run, tmp_path):
    """ED.6e-5: each joined -O column is that sample's O_s<k> file at its own strand k, and a
    file whose first line names another -s is refused (I5). Catches the strand index ignored,
    every sample read at one strand, and the per-file option check dropped."""
    _, _, align, _ = cluster_run
    source, problems = cluster.read(align)
    assert source is not None and not problems
    samples = source.tasks["sample"].to_list()
    strands = tuple((0, 1, 2, 1, 0)[: len(samples)])
    alignment = count.Alignment(tmp_path, False, count.Aligner.STAR, tuple(samples), False)
    dirs = {s: source.sample_dir(s) for s in samples}
    settings = {"O": count.FcOptions(strands, overlap=True)}
    kept = align / fc_compact.ANNOTATION
    joined = cluster.Counter(dirs, alignment, tmp_path / "keyed", source.gtf, kept)(settings)
    table = count.fc_genes(joined.matrix("O"), fraction=False)
    for sample, k in zip(samples, strands, strict=True):
        _, rows = fc_compact.read(
            str(dirs[sample] / "featurecounts" / f"O_s{k}{fc_compact.SUFFIX}")
        )
        assert table[f"{sample}.bam"].to_list() == [int(c) for _, c in rows[1:]], (sample, k)
    tampered = tmp_path / "tampered"
    writable_copy(align, tampered)
    path = tampered / "runs" / samples[0] / "featurecounts" / f"O_s0{fc_compact.SUFFIX}"
    with gzip.open(path, "rt") as handle:
        text = handle.read()
    with gzip.open(path, "wt") as handle:
        handle.write(text.replace('"-s" "0"', '"-s" "2"', 1))
    bad = {s: tampered / "runs" / s for s in samples}
    with pytest.raises(ProcessError, match="I5: .*O_s0"):
        cluster.Counter(bad, alignment, tmp_path / "keyed2", source.gtf, kept)(settings)


def test_the_compact_form_reads_back_as_featurecounts_wrote_it(tmp_path):
    """The cluster keeps Geneid and counts, gzipped, and the annotation columns once; ctx
    rebuilds the matrix from both. Catches a compact file that drops or shifts a row, an
    original deleted before the check, and a second annotation passed as the first."""
    first = tmp_path / "s0.featureCounts.txt"
    head = '# Program:featureCounts v2.1.1; Command:"featureCounts" "-s" "0" "a.bam" \n'
    columns = "Geneid\tChr\tStart\tEnd\tStrand\tLength\ta.bam\n"
    rows = "g1\tc1\t1\t9\t+\t9\t5\ng2\tc1;c1\t20;40\t30;50\t-;-\t22\t0\n"
    first.write_text(head + columns + rows)
    kept = tmp_path / fc_compact.ANNOTATION
    out = Path(fc_compact.compact(str(first), str(kept)))
    assert not first.exists()
    command, back = fc_compact.read(str(out))
    assert command == head.rstrip("\n") and back == [("Geneid", "a.bam"), ("g1", "5"), ("g2", "0")]
    assert fc_compact.annotation_rows(str(kept))[2] == [
        "g2",
        "c1;c1",
        "20;40",
        "30;50",
        "-;-",
        "22",
    ]
    other = tmp_path / "s1.featureCounts.txt"
    other.write_text(head + columns + rows.replace("\t9\t+", "\t8\t+"))
    with pytest.raises(ValueError, match="annotation columns differ"):
        fc_compact.compact(str(other), str(kept))
    assert other.exists()


def test_the_door_refuses_what_the_cluster_did_not_run(tmp_path):
    """ED.6e-4: --genome is needed without --align-dir; with it, the options for reads here or
    for another STAR setting are refused, each named. Catches --genome required always, and a
    --multimap-max silently ignored."""
    spec = tmp_path / "x.toml"
    assert any("--genome is needed" in p for p in process_run.option_problems(
        Options(spec, None)))  # fmt: skip
    assert not [p for p in process_run.option_problems(Options(spec, None, align_dir=tmp_path))
                if "--genome" in p]  # fmt: skip
    said = process_run.option_problems(
        Options(spec, None, multimap=5, reads_per_run=10, fetchngs_dir=tmp_path, align_dir=tmp_path)
    )
    for flag in ("--multimap-max", "--reads-per-run", "--fetchngs-dir"):
        assert any(p.startswith(f"{flag} is not read with --align-dir") for p in said), flag


def test_a_sheet_run_the_cluster_did_not_align_is_named(cluster_run):
    """ED.6e-4: every run of the sheet must be in tasks.tsv under the same sample, and every task
    run in the sheet. Catches a missing run passed over and a sample regrouped by the cluster."""
    root, _, align, _ = cluster_run
    source, _ = cluster.read(align)
    assert source is not None
    sheet = pl.read_parquet(
        next(root.joinpath("outputs", "01_gather").glob("*_samplesheet_v1.parquet"))
    )
    files, problems = cluster.files_of(sheet, source)
    assert not problems and files.height == sheet.height
    moved = sheet.with_columns(sample=pl.when(pl.col("sample") == sheet["sample"][0])
                               .then(pl.lit("other")).otherwise("sample"))  # fmt: skip
    _, problems = cluster.files_of(moved, source)
    assert any("on the cluster" in p for p in problems)
    _, problems = cluster.files_of(sheet.head(sheet.height - 1), source)
    assert any("not in the sheet" in p for p in problems)


def test_the_page_says_where_it_ran_and_with_what(cluster_run, process_root):
    """ED.6e-7: the cluster page's Method and Reproduce say Kebnekaise, carry versions.txt's tool
    versions and the cost from the time files; the laptop page says neither. Catches the notes
    dropped, versions typed instead of read, and the cluster words on the laptop page."""
    root, _, align, _ = cluster_run
    (page,) = root.joinpath(*OUT).glob(f"*_{DATASET}_overview_v1.html")
    text = unpacked(page.read_text())
    stated = (align / "versions.txt").read_text()
    assert "On Kebnekaise" in text and "versions.txt" in text
    for tool in ("STAR 2.7.11b", "fastp 1.1.0", "featureCounts v2.1.1"):
        assert tool in stated
    assert "STAR 2.7.11b, fastp 1.1.0" in text
    steps = newest(root, "cluster_steps")
    assert steps.filter(pl.col("step") == "fc")["files"].item() == 5 * 12 + 3
    (local,) = process_root[0].joinpath(*OUT).glob(f"*_{DATASET}_overview_v1.html")
    assert "Kebnekaise" not in unpacked(local.read_text())


def test_ctx_run_plans_process_from_an_align_dir(tmp_path):
    """ED.6e-8: ctx run plans process before qc, with --align-dir, only when it is given; the
    core bundles list process for the four projects. Catches process planned without the folder
    and the bundles left without it."""
    spec_path = dataset_file(tmp_path)
    spec = load(spec_path)
    base = pipeline.Options(DATASET, spec_path, tmp_path)
    assert "process" not in dict(pipeline.steps(tmp_path, spec, base))
    planned = pipeline.steps(tmp_path, spec, pipeline.Options(DATASET, spec_path, tmp_path,
                                                              align_dir=tmp_path))  # fmt: skip
    names = [m for m, _ in planned]
    assert names.index("process") < names.index("qc")
    args = dict(planned)["process"]
    assert args[args.index("--align-dir") + 1] == str(tmp_path)
    import tomllib  # noqa: PLC0415

    for bundle, slugs in (("myeloid", ["Non-cancerous-blood"]), ("lung", ["Lung-tissue-LUAD"]),
                          ("breast", ["Breast-tissue-BRCA", "BRCA-immune-cells-macrophages",
                                      "BRCA-immune-cells-monocytes"])):  # fmt: skip
        entries = tomllib.loads((ROOT / "conf" / "cohorts" / f"{bundle}.toml").read_text())
        found = {e["slug"]: e["modules"] for e in entries["datasets"]}
        for slug in slugs:
            assert "process" in found[slug], (bundle, slug)
