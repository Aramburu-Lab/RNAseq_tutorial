"""ctx.qc.multiqc_io on the rnaseq fixture: MultiQC's tables against the reports themselves."""

import json
import os
import re
import shutil
from pathlib import Path

import polars as pl
import pytest

from ctx import tools
from ctx.io.naming import lint
from ctx.qc import QcError, multiqc_io

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "rnaseq"
SAMPLES = ("RAP1_IAA_30M_REP1", "RAP1_UNINDUCED_REP1", "RAP1_UNINDUCED_REP2", "WT_REP1", "WT_REP2")
FATES = ("reads_written", "reads_too_short", "reads_too_long", "reads_too_many_n",
         "reads_discarded_untrimmed")  # fmt: skip


def aliases() -> dict[str, str]:
    """Every name the fixture's samples appear under in MultiQC's tables."""
    names = {}
    for sample in SAMPLES:
        for suffix in ("", "_trimmed", "_trimmed_1", "_trimmed_2"):
            names[sample + suffix] = sample
    return names


@pytest.fixture(scope="module")
def ran(tmp_path_factory):
    inputs = tmp_path_factory.mktemp("inputs")
    for part in ("trimgalore", "salmon"):
        shutil.copytree(FIXTURE / part, inputs / part)
    out = tmp_path_factory.mktemp("outputs") / "02_qc" / "multiqc"
    return multiqc_io.run([inputs], out, modules=("trim_galore", "salmon")), out


def unit(name: str) -> tuple[str, str]:
    """(sample, mate) of a report or table name: MultiQC drops `_trimmed` from single-end ones."""
    return aliases()[name], "2" if name.endswith("_trimmed_2") else "1"


def reports() -> dict[tuple[str, str], dict]:
    """The Trim Galore JSON reports, read here without MultiQC."""
    found = {}
    for path in (FIXTURE / "trimgalore").glob("*_trimming_report.json"):
        stem = path.name.removesuffix(".fastq.gz_trimming_report.json")
        found[unit(stem)] = json.loads(path.read_text())["read_processing"]
    return found


def test_trim_galore_table_equals_the_reports(ran):
    run, _ = ran
    table = multiqc_io.numeric(multiqc_io.table(run, "trim_galore"), ("total_reads", *FATES))
    by_name = {unit(r["Sample"]): r for r in table.iter_rows(named=True)}
    source = reports()
    assert set(by_name) == set(source)
    for name, values in source.items():
        for column in ("total_reads", *FATES):
            assert by_name[name][column] == values[column], (name, column)


def test_every_read_has_a_fate(ran):
    """total = written + too short + too long + too many N + discarded, exactly, per report."""
    run, _ = ran
    table = multiqc_io.numeric(multiqc_io.table(run, "trim_galore"), ("total_reads", *FATES))
    gap = table.select(pl.col("total_reads") - pl.sum_horizontal(FATES))
    assert gap.to_series().to_list() == [0.0] * table.height


def test_reads_written_reach_salmon(ran):
    """What Trim Galore wrote is what Salmon processed: read 1 for paired samples (D-022)."""
    run, _ = ran
    trim = multiqc_io.table(run, "trim_galore")
    names = multiqc_io.map_names(trim["Sample"], aliases())
    written = {}
    for row in trim.iter_rows(named=True):
        if not row["Sample"].endswith("_trimmed_2"):
            written[names[row["Sample"]]] = float(row["reads_written"])
    for sample in SAMPLES:
        meta = json.loads((FIXTURE / "salmon" / sample / "aux_info" / "meta_info.json").read_text())
        assert written[sample] == meta["num_processed"], sample
    salmon = multiqc_io.numeric(multiqc_io.table(run, "salmon"), ["num_processed"])
    assert dict(salmon.select("Sample", "num_processed").rows()) == {
        s: float(json.loads((FIXTURE / "salmon" / s / "aux_info" / "meta_info.json").read_text())
                 ["num_processed"]) for s in SAMPLES
    }  # fmt: skip


def test_unmapped_names_are_listed(ran):
    with pytest.raises(QcError, match=r"\['STRANGER', 'ZED_trimmed_1'\]"):
        multiqc_io.map_names(["WT_REP1", "STRANGER", "ZED_trimmed_1"], aliases())


def test_a_missing_report_stops_the_run(ran, tmp_path):
    run, _ = ran
    multiqc_io.require_sources(run, "trim_galore", SAMPLES, aliases())
    copy = tmp_path / "multiqc_data"
    shutil.copytree(run.data, copy)
    sources = (copy / "multiqc_sources.txt").read_text().splitlines()
    kept = [line for line in sources if "WT_REP2_trimmed" not in line]
    (copy / "multiqc_sources.txt").write_text("\n".join(kept) + "\n")
    with pytest.raises(QcError, match=r"no trim_galore report for \['WT_REP2'\]"):
        multiqc_io.require_sources(multiqc_io.MultiqcRun(copy, run.version), "trim_galore",
                                   SAMPLES, aliases())  # fmt: skip


def test_multiqc_folder_is_exempt_from_the_naming_lint(ran):
    _, out = ran
    scanned, violations = lint(out.parent)
    assert violations == []
    assert any(out.rglob("multiqc_trim_galore.txt"))
    (out.parent / "badly named.csv").write_text("x")
    assert lint(out.parent)[1] == ["badly named.csv"]
    assert scanned == 0


def test_relative_inputs_are_read_from_the_callers_directory(tmp_path, monkeypatch):
    """MultiQC runs in its output folder; a relative input once pointed nowhere from there."""
    shutil.copytree(FIXTURE / "salmon", tmp_path / "salmon")
    monkeypatch.chdir(tmp_path)
    run = multiqc_io.run([Path("salmon")], tmp_path / "out", modules=("salmon",))
    assert multiqc_io.table(run, "salmon").height == len(SAMPLES)


def test_native_folders_are_exempt_at_the_first_level_only(tmp_path):
    """`outputs/03_process/{multiqc,tools}/` keep the tools' names; the same names deeper are
    linted (D-037).

    Mutation caught: a native folder matched at any depth (the rule before M3)."""
    module = tmp_path / "03_process"
    names = ["tools/index/star_0123456789ab/SAindex", "tools/README.md",
             "multiqc/fc_s0/multiqc_data/multiqc_featurecounts.txt",
             "figures/tools/Aligned.out.bam", "figures/multiqc/multiqc_report.html"]  # fmt: skip
    for name in names:
        (module / name).parent.mkdir(parents=True, exist_ok=True)
        (module / name).write_text("x")
    assert lint(module) == (2, ["figures/multiqc/multiqc_report.html",
                                "figures/tools/Aligned.out.bam"])  # fmt: skip


SE = ROOT / "tests" / "fixtures" / "chr22_se"
SLICE = ROOT / "tests" / "fixtures" / "chr22" / "slice.fa"
# STAR's ReadsPerGene.out.tab: the N_ rows, then genes; columns unstranded, read 1 sense,
# read 2 sense. Every category is nonzero in some sample, so MultiQC keeps every column.
GENE_COUNTS = {
    "se_a": "N_unmapped\t3\t3\t3\nN_multimapping\t1\t1\t1\nN_noFeature\t4\t950\t2\n"
            "N_ambiguous\t437\t0\t388\nG1\t500\t40\t600\nG2\t62\t9\t11\n",
    "se_c": "N_unmapped\t30\t30\t30\nN_multimapping\t49\t49\t49\nN_noFeature\t0\t846\t0\n"
            "N_ambiguous\t437\t5\t372\nG1\t300\t50\t250\nG2\t184\t25\t299\n",
}  # fmt: skip


def test_star_gene_count_totals_are_plot_data(tmp_path):
    """STAR's per-column totals exist only as MultiQC plot data, as numbers (TODO.md M3.2).

    MultiQC leaves out a category that is zero in every sample (the next test), so a column
    is read only after checking it is there.

    Mutations caught: plot_table without its numeric cast; reading `multiqc_<name>.txt`."""
    for sample, text in GENE_COUNTS.items():
        (tmp_path / "in" / f"{sample}.ReadsPerGene.out.tab").parent.mkdir(exist_ok=True)
        (tmp_path / "in" / f"{sample}.ReadsPerGene.out.tab").write_text(text)
    run = multiqc_io.run([tmp_path / "in"], tmp_path / "out")
    reverse = multiqc_io.plot_table(run, "star_gene_counts_Reverse_Stranded")
    assert reverse.columns == ["Sample", "Overlapping Genes", "No Feature", "Ambiguous Features",
                               "Multimapping", "Unmapped"]  # fmt: skip
    assert reverse.sort("Sample").rows() == [("se_a", 611.0, 2.0, 388.0, 1.0, 3.0),
                                             ("se_c", 549.0, 0.0, 372.0, 49.0, 30.0)]  # fmt: skip
    same = multiqc_io.plot_table(run, "star_gene_counts_Same_Stranded").sort("Sample")
    assert same["Overlapping Genes"].to_list() == [49.0, 75.0]
    unstranded = multiqc_io.plot_table(run, "star_gene_counts_Unstranded")
    assert "No Feature" in unstranded.columns  # se_a's 4
    with pytest.raises(QcError, match="no star_gene_counts_Absent.txt"):
        multiqc_io.plot_table(run, "star_gene_counts_Absent")


def test_a_plot_value_that_is_not_a_number_stops_the_run(tmp_path):
    """Mutation caught: plot_table without its numeric cast (the text passes through)."""
    (tmp_path / "star_gene_counts_Unstranded.txt").write_text("Sample\tUnmapped\nse_a\t3x\n")
    with pytest.raises(QcError, match="not a number"):
        multiqc_io.plot_table(
            multiqc_io.MultiqcRun(tmp_path, "1.35"), "star_gene_counts_Unstranded"
        )


def test_a_category_zero_everywhere_is_left_out(tmp_path):
    """Pins MultiQC 1.35's behaviour, which M3.5's readers must allow for: no code of ours is
    under test, so no mutation of it applies; a MultiQC that kept the column fails here."""
    (tmp_path / "in").mkdir()
    (tmp_path / "in" / "se_c.ReadsPerGene.out.tab").write_text(GENE_COUNTS["se_c"])
    run = multiqc_io.run([tmp_path / "in"], tmp_path / "out")
    assert "No Feature" not in multiqc_io.plot_table(run, "star_gene_counts_Unstranded").columns


def test_hisat2_summaries_are_read(tmp_path):
    """HISAT2's `--new-summary --summary-file` is what MultiQC's hisat2 module reads; the
    sample is the summary's name before `.hisat2.txt`.

    Mutation caught: "hisat2" missing from MODULES (no table is written)."""
    index, aligned, logs = tmp_path / "index", tmp_path / "aligned", tmp_path / "logs"
    for folder in (index, aligned, logs):
        folder.mkdir()
    tools.run("hisat2-build", [str(SLICE), str(index / "slice")], cwd=index,
              log=logs / "build.log")  # fmt: skip
    for sample in ("se_a", "se_c"):
        # The log stays out of MultiQC's input: hisat2 prints the same summary to stderr.
        tools.run("hisat2", ["-x", str(index / "slice"), "-U", str(SE / f"{sample}.fq"),
                             "--new-summary", "--summary-file",
                             str(aligned / f"{sample}.hisat2.txt"), "-S", os.devnull],
                  cwd=aligned, log=logs / f"{sample}.log")  # fmt: skip
    run = multiqc_io.run([aligned], tmp_path / "multiqc")
    columns = ["unpaired_total", "unpaired_aligned_none", "unpaired_aligned_one",
               "unpaired_aligned_multi"]  # fmt: skip
    table = multiqc_io.numeric(multiqc_io.table(run, "hisat2"), columns)
    for row in table.iter_rows(named=True):
        summary = (aligned / f"{row['Sample']}.hisat2.txt").read_text()
        said = [int(n) for n in re.findall(r": (\d+)(?:\n| \()", summary)]
        assert [row[c] for c in columns] == said, row["Sample"]
    assert sorted(table["Sample"]) == ["se_a", "se_c"]
    assert table.filter(pl.col("Sample") == "se_a")["unpaired_total"].item() == 1000


def test_relative_dest_and_config_are_resolved(tmp_path, monkeypatch):
    """MultiQC runs in dest: with a relative -c it exits 2 (M3 review 1), with a relative -o
    it exits 1.

    Mutations caught: the config passed as given; dest passed as given."""
    (tmp_path / "in").mkdir()
    (tmp_path / "in" / "se_a.ReadsPerGene.out.tab").write_text(GENE_COUNTS["se_a"])
    (tmp_path / "config.json").write_text(json.dumps({"report_comment": "relative config"}))
    monkeypatch.chdir(tmp_path)
    run = multiqc_io.run([Path("in")], Path("out") / "multiqc", config=Path("config.json"))
    assert run.data == tmp_path / "out" / "multiqc" / "multiqc_data"
    assert "relative config" in (tmp_path / "out" / "multiqc" / "multiqc_report.html").read_text()


def test_at_once_runs_the_calls_together_and_keeps_their_order():
    """Three calls that each wait for the other two finish only if they run together; the
    one that ends last comes first, as given; an error comes out after every call has ended.

    Mutations caught: the calls run one at a time (the barrier breaks); the results in the
    order they ended; an error swallowed, or raised before the other calls end."""
    import threading  # noqa: PLC0415
    import time  # noqa: PLC0415

    together, ended = threading.Barrier(3, timeout=10), []

    def call(name: str, wait: float) -> str:
        together.wait()
        time.sleep(wait)
        ended.append(name)
        return name.upper()

    work = {"c": lambda: call("c", 0.3), "a": lambda: call("a", 0.0), "b": lambda: call("b", 0.1)}
    assert list(multiqc_io.at_once(work, 3).items()) == [("c", "C"), ("a", "A"), ("b", "B")]
    assert ended == ["a", "b", "c"]

    def fails() -> str:
        raise QcError("planted")

    def slow() -> None:
        time.sleep(0.2)
        ended.append("slow")

    ended.clear()
    with pytest.raises(QcError, match="planted"):
        multiqc_io.at_once({"fails": fails, "slow": slow}, 2)
    assert ended == ["slow"]


# --- nf-core's own MultiQC tables, read in place (D-116) ------------------------------------

STAR = ROOT / "tests" / "fixtures" / "rnaseq_star"
NFCORE_DATA = STAR / "multiqc" / "star_salmon" / "multiqc_report_data"


def _report_copy(tmp_path: Path) -> Path:
    """The aligned fixture's MultiQC tables in a copy of the outdir's layout."""
    target = tmp_path / "out" / "multiqc" / "star_salmon" / "multiqc_report_data"
    shutil.copytree(NFCORE_DATA, target)
    return tmp_path / "out"


def test_nfcore_reads_the_tables_where_the_pipeline_wrote_them():
    """No MultiQC runs: the run's folder is nf-core's own, its version from MultiQC's log."""
    run = multiqc_io.nfcore(STAR)
    assert run.data == NFCORE_DATA and run.version == tools.TOOLS["multiqc"].version
    star = multiqc_io.table(run, "star")
    assert set(star["Sample"]) == set(SAMPLES)
    assert not (STAR / "multiqc" / "star_salmon" / "multiqc_data").exists()


def test_nfcore_refuses_another_multiqc(tmp_path):
    """Mutation: the version line ignored; MultiQC 1.34's columns are read as 1.35's."""
    out = _report_copy(tmp_path)
    log = out / "multiqc" / "star_salmon" / "multiqc_report_data" / "multiqc.log"
    log.write_text(log.read_text().replace("This is MultiQC v1.35", "This is MultiQC v1.34"))
    with pytest.raises(QcError, match=r"MultiQC 1\.34 wrote"):
        multiqc_io.nfcore(out)
    log.unlink()
    with pytest.raises(QcError, match="names no MultiQC version"):
        multiqc_io.nfcore(out)


def test_nfcore_without_tables_names_the_copy_list(tmp_path):
    with pytest.raises(QcError, match="ALIGNED_FILES"):
        multiqc_io.nfcore(tmp_path)


def test_points_read_every_cell_as_measured():
    """Qualimap's normalised profile: 100 positions per sample, each "(x, y)" cell as written."""
    run = multiqc_io.nfcore(STAR)
    found = multiqc_io.points(run, "qualimap_gene_coverage_profile_Normalised")
    assert found.group_by("Sample").len()["len"].to_list() == [100] * len(SAMPLES)
    raw = pl.read_csv(NFCORE_DATA / "qualimap_gene_coverage_profile_Normalised.txt",
                      separator="\t", infer_schema=False)  # fmt: skip
    first = raw.row(0)
    x, y = (float(v) for v in first[1].strip("()").split(","))
    row = found.filter(pl.col("Sample") == first[0]).row(0, named=True)
    assert (row["x"], row["y"]) == (x, y)


def test_points_refuse_a_cell_that_is_no_point(tmp_path):
    """Mutation: a cell parsed with findall drops "(1.0, )" without a word."""
    out = _report_copy(tmp_path)
    path = (
        out
        / "multiqc"
        / "star_salmon"
        / "multiqc_report_data"
        / "preseq_complexity_plot_molecules.txt"
    )
    path.write_text(path.read_text().replace("(1.0, 0.3569073)", "(1.0, )"))
    with pytest.raises(QcError, match="not an"):
        multiqc_io.points(multiqc_io.nfcore(out), "preseq_complexity_plot_molecules")


def test_ragged_takes_each_sample_its_own_x():
    """dupRadar's table: a header of x values, each sample filling its own cells only."""
    found = multiqc_io.ragged(multiqc_io.nfcore(STAR), "dupradar")
    raw = pl.read_csv(NFCORE_DATA / "multiqc_dupradar.txt", separator="\t", infer_schema=False)
    filled = sum(v is not None and v != "" for row in raw.rows() for v in row[1:])
    assert found.height == filled and set(found["Sample"]) == set(SAMPLES)
    per_row = {row[0]: sum(v not in (None, "") for v in row[1:]) for row in raw.rows()}
    assert dict(found.group_by("Sample").len().rows()) == per_row
