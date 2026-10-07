"""`ctx process` end to end: HISAT2 on chr22_se, its page and files, and STAR where it maps.

The HISAT2 run is ``conftest.process_root``, once per session; the refusals at the door are
tests/test_process_door.py. Tests marked needs_star skip only when the STAR probe fails on
darwin-arm64 with CI unset, and fail everywhere else (R6.2).
"""

import gzip
import html as markup
import json
import os
import re
import shutil
from pathlib import Path

import polars as pl
import pytest
from conftest import CAP, ctx_run
from process_root import (
    DATASET,
    ROOT,
    SAMPLES,
    SE,
    SLICE,
    dataset_file,
    plant,
)
from test_cli import _offline

from ctx.io.naming import lint
from ctx.qc import registry
from ctx.report import prose, text
from ctx.runtime import BUDGETS

PAIRED = ROOT / "tests" / "fixtures" / "chr22"
# The genes no same-strand neighbour overlaps: featureCounts gives each its true reads.
CLEAN = ("SF3A1", "CCDC157", "RNF215", "ENSG00000273350", "ENSG00000303189")
# A STAR that maps, where this machine has one off the pin's path (D-036's route).
STAR_ENV = {k: v for k, v in os.environ.items() if k == "CTX_TOOLS"}
FIGURES = ("mapping", "featured_settings", "strand_assigned",
           *(f"{k}_{c}" for k in ("counters", "agreement")
             for c in ("strand", "overlap", "fraction", "multimap")))  # fmt: skip


def out(root: Path) -> Path:
    """The module's output folder."""
    return root / "outputs" / "03_process"


def table(root: Path, what: str) -> pl.DataFrame:
    """One of the run's tables, from its Parquet file."""
    return pl.read_parquet(next(out(root).glob(f"*_{DATASET}_{what}_v1.parquet")))


def page(root: Path) -> str:
    """The run's overview page."""
    return next(out(root).glob(f"*_process_{DATASET}_overview_v1.html")).read_text()


def visible(html: str) -> str:
    """The text a reader sees: no script, style or tag, entities decoded, spaces collapsed."""
    shown = re.sub(r"<(script|style)\b.*?</\1>", " ", html, flags=re.S)
    return re.sub(r"\s+", " ", markup.unescape(re.sub(r"<[^>]+>", " ", shown)))


def captions(html: str) -> list[str]:
    """Every figure caption on the page, in page order, as a reader sees it."""
    return [re.sub(r" ([,.;])", r"\1", visible(c)).strip()
            for c in re.findall(r"<figcaption>(.*?)</figcaption>", html, flags=re.S)]  # fmt: skip


def template(key: str) -> re.Pattern[str]:
    """A caption template of content/process.toml as a pattern: marks as their display text,
    each placeholder any text."""
    from ctx.report.marks import plain  # noqa: PLC0415

    escaped = re.escape(plain(text.raw("process", f"figures.{key}.caption")))
    return re.compile(re.sub(r"\\\{[a-z_]+\\\}", "(.+?)", escaped) + "$")


def test_every_output_is_named_described_and_offline(process_root):
    """Eleven figures of seven files each; no file written twice; a README beside every CSV and
    Parquet; the naming lint clean outside the tool folders; tools/README.md; the page offline
    and under 15 MB.

    Mutations caught: one comparison dropped from tables.COMPARISONS (eight figures);
    tools/README.md not written; the strand figure named like the strandedness table (its CSV
    and README overwritten by the table's, the path printed twice)."""
    root, _, done = process_root
    written = done.stdout.split()
    assert len(written) == len(set(written))
    folder = out(root)
    for name in FIGURES:
        assert len(list(folder.glob(f"*_{DATASET}_{name}_v1*"))) == 7, name
    assert not list(folder.glob("*_counters_tool_v1*"))  # no STAR, no tool comparison
    for data in [*folder.glob("*.csv"), *folder.glob("*.parquet")]:
        assert (folder / f"{data.name.rsplit('.', 1)[0]}_README.md").is_file(), data.name
    assert lint(folder)[1] == []
    assert (folder / "tools" / "README.md").read_text().startswith("# tools/")
    html = page(root)
    _offline(html)
    # The text a reader sees: no script or style, no tag (a filter's line template is one).
    shown = re.sub(r"<[^>]+>", " ", re.sub(r"<(script|style)\b.*?</\1>", "", html, flags=re.S))
    assert "built-in method" not in shown and not re.search(r"\{[a-z_]+\}", shown)


def test_each_sample_is_counted_at_its_own_call(process_root):
    """se_a to se_d are called reverse and se_e unstranded, each counted at its call, and the
    Salmon twin agrees; gene_counts holds each sample's count at that strand.

    Mutation caught: gene_counts read every sample at -s 2 (se_e's column is its -s 2 count,
    not the called -s 0 the strand comparison holds)."""
    root, _, _ = process_root
    strand = table(root, "strandedness")
    calls = dict(strand.select("sample", "call").rows())
    assert calls == {**dict.fromkeys(SAMPLES, "reverse"), "se_e": "unstranded"}
    assert (strand["used"] == strand["call"]).all() and strand["agrees"].all()
    counts, long = table(root, "gene_counts"), table(root, "gene_comparison")
    for sample in SAMPLES:
        called = long.filter((pl.col("comparison") == "strand") & (pl.col("sample") == sample))
        joined = counts.select("gene_id", sample).join(called.select("gene_id", "b"), on="gene_id")
        assert (joined[sample].cast(pl.Float64) == joined["b"]).all(), sample


def test_counts_match_the_planted_truth(process_root):
    """On the five genes no same-strand gene overlaps, a stranded sample falls short of its
    truth by at most the reads fastp dropped and the aligner left unaligned, so by nothing for
    se_a, se_b and se_d; se_e falls short by at most its -s 0 ambiguous reads.

    Mutation caught: gene_counts from the -s 0 matrix (se_a's CCDC157 loses the reads it shares
    with SF3A1 on the other strand)."""
    root, _, _ = process_root
    truth = pl.read_csv(SE / "truth.tsv", separator="\t").group_by("sample", "gene_id").agg(
        pl.col("n_reads").sum())  # fmt: skip
    counts = table(root, "gene_counts").filter(pl.col("gene_name").is_in(CLEAN))
    long = counts.unpivot(index=["gene_id", "gene_name", "gene_type"], variable_name="sample",
                          value_name="count")  # fmt: skip
    short = long.join(truth, on=["sample", "gene_id"], how="left").group_by("sample").agg(
        (pl.col("n_reads").fill_null(0) - pl.col("count")).abs().sum().alias("off"))  # fmt: skip
    off = dict(short.rows())
    alignment = table(root, "alignment")
    lost = {s: f - u - m for s, f, u, m in alignment.select(
        "sample", "fastp_passed", "unique", "multi").rows()}  # fmt: skip
    dropped = {s: 1000 - f for s, f in alignment.select("sample", "fastp_passed").rows()}
    for sample in ("se_a", "se_b", "se_c", "se_d"):
        assert off[sample] <= dropped[sample] + lost[sample], sample
    assert off["se_a"] == off["se_b"] == off["se_d"] == 0
    ambiguous = dict(table(root, "strandedness").select("sample", "ambiguous_s0").rows())
    assert 0 < off["se_e"] <= ambiguous["se_e"]


def test_half_the_genome_keeps_its_genes_whole(process_root):
    """se_a on the first half of the slice: SF3A1 and CCDC157 keep all 246 reads, nothing gained
    or lost; 8 genes go with the 366 reads they held; uniquely mapped falls from 100 %.

    Mutations caught: the cut alignment made against the whole index (unique share unchanged);
    the half counted against the whole annotation (I4 stops the run)."""
    root, _, _ = process_root
    row = table(root, "incomplete").row(0, named=True)
    assert (row["genes_kept"], row["genes_removed"]) == (2, 8)
    assert row["kept_reads_main"] == row["kept_reads_cut"] == 246
    assert (row["gained"], row["lost"], row["removed_reads_main"]) == (0, 0, 366)
    assert row["unique_pct_cut"] < row["unique_pct_main"] == 100.0
    settings = table(root, "settings")
    assert settings["setting"].to_list() == ["main", "cut", "untrimmed", "unspliced"]
    assert not settings["derived"].any()


def test_star_only_columns_keep_their_types_under_hisat2(process_root):
    """The HISAT2 run's alignment and reference tables carry STAR's columns as typed nulls, so a
    later reader (M2b) reads one shape whatever the aligner. The reference table's columns are
    the M3 plan's (section 9), less build_peak_bytes (D-041), under the plan's names.

    Mutations caught: the casts dropped (the columns' type is Null, the Parquet holds no type);
    the index column named ``index`` again, not ``index_key``."""
    root, _, _ = process_root
    alignment, reference = table(root, "alignment"), table(root, "reference")
    assert reference.columns == [
        "reference", "file", "sha256", "bytes", "contigs", "length_nt", "n_bases",
        "sa_index_bases", "genes_subset", "protein_coding_subset", "genes_genome",
        "annotation_release", "annotation_sha256", "start_codons", "atg_share", "index_key",
        "built", "build_s",
    ]  # fmt: skip
    assert alignment.schema["too_many"] == pl.Int64 and alignment["too_many"].is_null().all()
    assert alignment.schema["splices"] == pl.Float64 and alignment["splices"].is_null().all()
    assert reference.schema["sa_index_bases"] == pl.Int64
    assert reference["sa_index_bases"].is_null().all()


def test_salmon_s_share_is_context_not_a_ceiling(process_root):
    """The planted twin puts 900 of 1,000 counts on the slice's genes for se_a to se_d (0.9, and
    0.72 of every read at 80 % mapped) and 500 for se_e (0.5, 0.4). Observed is featureCounts'
    Assigned at the called strand over input (se_a 612 of 1,000); se_e's is above its Salmon
    share, as on a real partial reference, and the page says so beside the table.

    Mutations caught: observed taken over the Assigned reads of every strand (-s 0 for all);
    the ratio inverted (se_e below one, the note saying 4 of 5 above)."""
    root, _, _ = process_root
    found = table(root, "context")
    assert found["expected_share"].to_list() == pytest.approx([0.9] * 4 + [0.5])
    assert found["expected_read_share"].to_list() == pytest.approx([0.72] * 4 + [0.4])
    observed = dict(found.select("sample", "observed_assigned_share").rows())
    assert observed["se_a"] == pytest.approx(0.612)
    ratio = found["observed_assigned_share"] / found["expected_read_share"]
    assert found["observed_over_expected"].to_list() == pytest.approx(ratio.to_list())
    above = dict(found.select("sample", pl.col("observed_over_expected") > 1).rows())
    assert above == {**dict.fromkeys(SAMPLES, False), "se_e": True}
    shown = visible(page(root))
    assert "above one in 1 of 5 samples" in shown
    assert "Salmon's share is no ceiling" in shown
    assert "median 72.0 % of each sample's reads" in shown


def test_a_moved_rerun_reuses_every_tool_folder(process_rerun, process_root):
    """A copy of the root rerun at the default settings: every index reused, no file under
    tools/ written beside the first run's (a rebuild goes into a new keyed folder, never over an
    old one), and not one of the first run's files rewritten. The folder table says so for each
    index, the trim, every alignment and every count, where the first run built them all.

    Mutations caught: the trim key fed the FASTQ paths, not their SHA-256 (``align.trim``): the
    copied root trims and aligns again into new folders; the folder table reading every
    folder as built."""
    root, html, before = process_rerun
    reference = table(root, "reference")
    assert reference["built"].to_list() == [False, False]
    tools = out(root) / "tools"
    after = {str(p.relative_to(tools)): p.stat().st_mtime_ns for p in tools.rglob("*")
             if p.is_file()}  # fmt: skip
    assert sorted(set(after) - set(before)) == []
    rewritten = [p for p, t in before.items() if after.get(p) != t and p != "README.md"]
    assert rewritten == []
    folders, first = table(root, "folders"), table(process_root[0], "folders")
    assert first["built"].all() and not folders["built"].any()
    assert folders["step"].to_list() == ["index_main", "index_cut", "trim", "align_main",
                                         "align_cut", "align_unspliced", "count_plain",
                                         "count_variants", "count_cut"]  # fmt: skip
    assert set(folders["folder"]) <= set(first["folder"])
    assert visible(html).count(text.word("process", "reused")) >= folders.height


def test_a_subsets_page_says_it_is_the_laptop_exercise(process_rerun, process_root):
    """The rerun's dataset file declares subset_of: its page says first that it is the optional
    laptop exercise on that parent's runs (D-050); the full run's page does not.

    Mutation caught: the note left out of ``page.notes`` (no subset page says what it is)."""
    _, html, _ = process_rerun
    note = "This page is the optional laptop exercise, on a subset of chrse_cohort"
    assert note in html and html.index(note) < html.index("The reference each alignment used")
    assert "optional laptop exercise" not in page(process_root[0])


def test_the_run_is_within_track_a(process_root):
    """The page's own wall time and peak RSS, self plus its largest child, within Track A. Each
    is measured, not left at zero: the run's tools (fastp, HISAT2, featureCounts) are children.

    Mutations caught: the Reproduction given wall_seconds 0, or peak_rss_children 0 (the page
    states no cost, or none of its tools', and Track A's budget is never checked)."""
    root, _, _ = process_root
    html = page(root)
    wall = float(re.search(r'<td class="wall">([\d.]+) s</td>', html).group(1))
    rss = [int(re.search(rf'<td class="{c}">[\d.]+ MiB \((\d+) bytes\)</td>', html).group(1))
           for c in ("rss-self", "rss-children")]  # fmt: skip
    assert 0 < wall < BUDGETS["A"].minutes * 60
    assert all(r > 0 for r in rss) and sum(rss) < BUDGETS["A"].ram_gb * 10**9


# Each figure's caption key, in page order: F1, F1b, F2, then the scatters and the plots.
CAPTION_KEYS = ("mapping", "settings", "strandedness", *("counters",) * 4, *("agreement",) * 4)
REGISTRY = registry.load(ROOT / "conf" / "qc_metrics.toml")


def test_captions_come_from_content_and_state_no_test(process_root):
    """R6.1: every caption is its content key filled, says "no test" and names no test it did
    not run; the strand caption gives the registry's own thresholds; none is over CONVENTIONS'
    40 words as the reader sees it, units and numbers included, even with the preview's
    "Made-up data." put before it (D-050).

    Mutations caught: a caption written by ``ctx.viz.annotate`` (Mann-Whitney, no "no test");
    the strand caption's two thresholds swapped; the Salmon sentence appended to the mapping
    caption again (46 words); the agreement caption at its first length (44 words); the
    counters caption at 40 words, 42 on the preview."""
    root, _, _ = process_root
    found = captions(page(root))
    assert len(found) == len(CAPTION_KEYS)
    made_up = prose.words(text.word("process", "preview.caption"))
    for key, caption in zip(CAPTION_KEYS, found, strict=True):
        assert template(key).match(caption), (key, caption)
        assert "no test" in caption.lower() and "Mann-Whitney" not in caption, key
        assert prose.words(caption) + made_up <= prose.CAPTION_WORDS, (key, prose.words(caption))
    metric = REGISTRY.metrics["gene_strand_bias"]
    assert f"stranded {metric.stranded}, unstranded {metric.unstranded}." in found[2]


def test_the_key_number_counts_the_genes_multimappers_change(process_root):
    """The key number "Genes counted differently": the multimap comparison's genes less its
    identical ones.

    Mutation caught: the key number read from the overlap comparison (7 of 9 on the fixture)
    while its detail still says multimappers."""
    root, _, _ = process_root
    every = (
        table(root, "concordance")
        .filter((pl.col("comparison") == "multimap") & (pl.col("unit") == "all"))
        .row(0, named=True)
    )
    shown = re.search(r"<dt>Genes counted differently</dt><dd class=\"key-value\">([\d,]+)<",
                      page(root))  # fmt: skip
    assert shown and int(shown.group(1).replace(",", "")) == every["n_genes"] - every["n_identical"]


def test_a_hisat2_page_names_no_star_column_or_option(process_root):
    """STAR's ReadsPerGene columns are named only on a STAR page, in the strand figure's reading
    guide, and the term "Multimapper" names no STAR option, since HISAT2 takes none.

    Mutations caught: the STAR reading line added whatever the aligner; the term naming
    --outFilterMultimapNmax again."""
    root, _, _ = process_root
    html = page(root)
    star = text.word("process", "figures.strandedness.reading_star")
    assert "STAR column" not in html and star not in visible(html)
    terms = re.search(r'<details class="terms">(.*?)</details>', html, flags=re.S)
    assert terms and "outFilterMultimapNmax" not in terms.group(1)


def test_records_stand_beside_reads_and_splices_per_aligner(process_root):
    """The page's alignment table gives each sample's flagstat records beside its primary
    mapped %, from the alignment file; STAR's splice columns read "not defined" under HISAT2.

    Mutation caught: the alignment table left off the page (the lesson on records against
    reads has nothing to show)."""
    root, _, _ = process_root
    shown = visible(page(root))
    title = text.word("process", "tables.alignment.page.title")
    assert title in shown
    row = table(root, "alignment").filter(pl.col("sample") == "se_c").row(0, named=True)
    undefined = text.word("process", "undefined")
    cells = (
        f"se_c {row['unique_pct']:.1f} {undefined} {undefined} {undefined} "
        f"{row['records']:,} {row['mapped_records_pct']:.1f} {row['primary_mapped_pct']:.1f}"
    )
    assert cells in shown[shown.index(title) :]


def test_the_cap_and_the_reads_each_run_kept_are_reported(process_root):
    """Under --reads-per-run every run's realised count is fastp's reads before filtering, in
    the alignment table and in the note that opens the results.

    Mutation caught: fastp_raw summed from fastp's passed reads (se_c's 997, not 1,000)."""
    root, _, _ = process_root
    alignment = table(root, "alignment")
    assert alignment["fastp_raw"].to_list() == [int(CAP)] * len(alignment)
    assert (alignment["fastp_raw"] >= alignment["fastp_passed"]).all()
    shown = visible(page(root))
    assert f"Each run was capped at {int(CAP):,} reads" in shown
    assert f"kept {int(CAP):,} to {int(CAP):,} reads" in shown


def test_hover_and_filters_name_every_column_in_words(process_root):
    """No hover line and no filter names a column by its identifier; the strand figure's
    forward share is printed to three decimals, not seventeen digits.

    The bar filters' live lines print no median: on bars stacked to 100 % it is always 100 %,
    and on the strand bars it pools three settings (P2S.6 (d)).

    Mutations caught: the strand bars built without ``labels`` (hover "forward_share:", filter
    "call"); the scatters without them (hover "a:" and "b:"); the strand bars built with
    ``summary_median`` left on (the line says "median 56.3 %")."""
    root, _, _ = process_root
    html = page(root)
    lines = re.findall(r'class="calls-summary filter-summary"[^>]*data-template="([^"]*)"'
                       r'[^>]*>([^<]*)<', html)  # fmt: skip
    assert lines and not [line for line in lines if "median" in " ".join(line)], lines
    hovers = re.findall(r'"hovertemplate":"(.*?)"', html)
    assert hovers
    names = (
        "a|b|forward_share|no_feature|ambiguous|setting|call|reads|input_reads|gene_\\w+|outcome"
    )
    raw = rf"(?:^|\\u003cbr\\u003e)({names}): "
    assert not [h for h in hovers if re.search(raw, h)]
    assert "Forward share: " in html and "Plain count: " in html
    labels = re.findall(r'<label>([^<]*)<select data-column="([^"]+)"', html)
    assert ("Strand call", "call") in labels and all(word != column for word, column in labels)
    shares = pl.read_csv(next(out(root).glob(f"*_{DATASET}_strand_assigned_v1.csv")))
    assert (shares["forward_share"] == shares["forward_share"].round(3)).all()


def test_a_partial_reference_is_said_to_be_one_and_a_whole_one_is_not(process_root):
    """The fixture's genome lacks the annotation's contig "elsewhere": the caption and the
    Interpretation say genome-wide bands do not apply. A reference holding every gene of the
    annotation reads the whole wording, and both versions keep the prose rules.

    Mutation caught: ``page.partial`` always true (a whole genome's page says its reference is
    partial)."""
    from ctx.process import Aligner  # noqa: PLC0415
    from ctx.process import page as pages  # noqa: PLC0415

    root, _, _ = process_root
    shown = visible(page(root))
    assert text.word("process", "figures.mapping.bands.partial") in shown
    assert "so no genome-wide band applies" in shown
    reference = table(root, "reference")
    whole = reference.with_columns(genes_subset=pl.col("genes_genome"))
    for given, partial in ((reference, True), (whole, False)):
        frames = pages.Frames(*([None] * 8), reference=given, incomplete=None, context=None)
        assert pages.partial(frames) is partial
    values = {"samples": 5, "runs": 5, "genome": "genome.fa", "length": "3,100,000,000",
              "release": "GENCODE 50", "low": "90.1", "high": "95.2"}  # fmt: skip
    for partial in (True, False):
        shown_as = pages.Shown(Aligner.HISAT2, "se_a", pl.DataFrame({"sample": ["se_a"]}),
                               "25 contigs", REGISTRY.metrics["gene_strand_bias"], 11,
                               partial=partial)  # fmt: skip
        interpretation = pages.prose(shown_as, values)["interpretation"]
        assert ("so no genome-wide band applies" in interpretation) is partial
        # Filled, its numbers are the placeholders' own: every other rule still holds.
        checked = prose.check_text("interpretation", interpretation, 120, numbers=False)
        assert list(checked) == []


def test_capped_counts_a_pair_once():
    """fastp counts both mates, the cap one pair: a paired run's 1,000 reads are 500 kept.

    Mutation caught: the halving dropped (the note says 1,000 kept of a 500-pair cap)."""
    from ctx.process import page as pages  # noqa: PLC0415

    runs = pl.DataFrame({"run_accession": ["r1", "r2"], "raw_reads": [1000, 998],
                         "passed_reads": [990, 990]})  # fmt: skip
    note = pages.capped(500, runs, paired=True)
    assert "kept 499 to 500 reads, median 499.5" in note


@pytest.mark.needs_star
def test_star_counts_equal_featurecounts_on_single_end_reads(tmp_path):
    """--aligner star on se_a to se_d: twelve figures, the tool comparison among them; per
    sample, STAR's column at the called strand equals featureCounts on every gene (rho 1), and
    the derived limit-1 bar stands beside the real settings. I2 held, since the run finished.

    Mutation caught: STAR's forward column read for a reverse sample (``tables._column`` taking
    the strand's neighbour): every gene differs."""
    root = tmp_path / "root"
    root.mkdir()
    plant(root, SAMPLES[:4])
    spec = dataset_file(root, tmp_path)
    done = ctx_run("--root", str(root), "process", "--dataset", DATASET, "--dataset-file",
                   str(spec), "--genome", str(SLICE), cwd=root,
                   env=STAR_ENV)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    for name in ("counters_tool", "agreement_tool"):
        assert len(list(out(root).glob(f"*_{DATASET}_{name}_v1*"))) == 7, name
    concordance = table(root, "concordance").filter(pl.col("comparison") == "tool")
    per_sample = concordance.filter(pl.col("unit") != "all")
    assert per_sample["spearman"].to_list() == [1.0] * 4
    assert (per_sample["n_identical"] == per_sample["n_genes"]).all()
    settings = table(root, "settings")
    assert settings.filter(pl.col("derived"))["setting"].to_list() == ["limit_1"]


@pytest.mark.needs_star
def test_star_and_featurecounts_differ_on_pairs(tmp_path):
    """The paired fixture through --aligner star: STAR counts a pair only when both mates fall on
    one gene, featureCounts --countReadPairs also when one mate is ambiguous. MTFP1 70 against
    264, SEC14L2 199 against 262, ENSG00000249590 0 against 3, CCDC157 506 in both; rho below 1.

    Mutation caught: --countReadPairs dropped from -p (CCDC157 counts 1,012 reads)."""
    root = tmp_path / "root"
    root.mkdir()
    sheet = plant(root, ("se_a",))
    fastq = root / "data" / "raw" / DATASET / "fastq"
    for mate in (1, 2):
        with (PAIRED / f"reads_R{mate}.fq").open("rb") as src, gzip.open(
                fastq / f"pe_{mate}.fastq.gz", "wb") as dst:  # fmt: skip
            shutil.copyfileobj(src, dst)
    sheet_path = next((root / "outputs" / "01_gather").glob("*_samplesheet_v1.parquet"))
    sheet.with_columns(fastq_1=pl.lit("fastq/pe_1.fastq.gz"), fastq_2=pl.lit("fastq/pe_2.fastq.gz"),
                       status_2=pl.lit("verified")).write_parquet(sheet_path)  # fmt: skip
    spec = dataset_file(root, tmp_path)
    done = ctx_run("--root", str(root), "process", "--dataset", DATASET, "--dataset-file",
                   str(spec), "--genome", str(SLICE), "--settings", "", cwd=root,
                   env=STAR_ENV)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    tool = table(root, "gene_comparison").filter(pl.col("comparison") == "tool")
    pairs = dict(tool.select("gene_name", pl.concat_list("a", "b")).rows())
    assert pairs["MTFP1"] == [70.0, 264.0] and pairs["SEC14L2"] == [199.0, 262.0]
    assert pairs["ENSG00000249590"] == [0.0, 3.0] and pairs["CCDC157"] == [506.0, 506.0]
    whole = (pl.col("comparison") == "tool") & (pl.col("unit") == "all")
    every = table(root, "concordance").filter(whole).row(0, named=True)
    assert every["spearman"] < 1


def test_the_settings_keep_their_labels_however_many_samples(tmp_path):
    """35 samples, more bars than an axis labels: the per-sample figure drops its labels, while
    the featured sample's settings, a figure of their own, name every bar in the SVG.

    Mutation caught: the setting bars drawn in the per-sample figure again (39 bars, no label
    on any setting; the course preview's 33 samples lost theirs so)."""
    from datetime import date  # noqa: PLC0415
    from types import SimpleNamespace  # noqa: PLC0415

    from ctx.de.run import Writer  # noqa: PLC0415
    from ctx.process import Aligner, blocks  # noqa: PLC0415
    from ctx.process import page as pages  # noqa: PLC0415

    def reads(unique: int) -> dict[str, int]:
        return {"unique": unique, "multi": 10, "unmapped": 990 - unique, "input": 1000}

    samples = [f"s{i:02d}" for i in range(35)]
    found = SimpleNamespace(
        outcomes={s: reads(900) for s in samples},
        settings={"main": reads(900), "cut": reads(400), "untrimmed": reads(890),
                  "unspliced": reads(800)},
        alignment=pl.DataFrame({"sample": samples, "input_reads": [1000] * 35,
                                "unique_pct": [90.0] * 35}),
    )  # fmt: skip
    shown = pages.Shown(Aligner.HISAT2, "s00", pl.DataFrame({"sample": samples}), "chr22",
                        REGISTRY.metrics["gene_strand_bias"], 10)  # fmt: skip
    writer = Writer(tmp_path / "03_process", "demo", date(2026, 10, 1), {})
    blocks.mapping_figure(writer, found, shown)
    blocks.settings_figure(writer, found, shown)
    svg = {what: next(writer.out.glob(f"*_demo_{what}_v1_white.svg")).read_text()
           for what in ("mapping", "featured_settings")}  # fmt: skip
    assert ">s01<" not in svg["mapping"]
    for key in ("main", "cut", "untrimmed", "unspliced"):
        assert f">{text.word('process', f'settings.{key}')}<" in svg["featured_settings"], key


def _y_label(svg: Path) -> list[str]:
    """The y-axis label of a figure file, line by line: its texts at the axis-title size, turned
    a quarter, in the order they are drawn."""
    import xml.etree.ElementTree as ET  # noqa: PLC0415, N817

    from ctx.viz.style import SIZES  # noqa: PLC0415

    size, texts = f"font-size: {SIZES['axis_label']:g}px", "{http://www.w3.org/2000/svg}text"
    return ["".join(t.itertext()) for t in ET.parse(svg).getroot().iter(texts)
            if "rotate(-90" in t.get("transform", "") and size in t.get("style", "")]  # fmt: skip


def test_axis_titles_hold_a_phone_and_y_titles_are_the_files(process_root):
    """Every x-axis title on the page, line by line, within ``axes.SCREEN_TITLE_CHARS``: plotly
    never wraps one, and at 375 px a longer line ran past the figure (M3.7, measured in
    headless Chromium). Every y title is the figure file's y label, line by line (item 2, D10):
    the page's own PDF and SVG of a view copy it, and a y title broken by the phone rule
    carried that break into them.

    Mutations caught: ``screen_title`` disabled (lines of 29 to 38 characters); the overlap
    comparison named "Count with every overlapped gene (-O)" again; ``screen_title`` applied
    to a y title."""
    from ctx.viz.axes import SCREEN_TITLE_CHARS  # noqa: PLC0415

    root, _, _ = process_root
    for what in FIGURES:
        (html,) = out(root).glob(f"*_{DATASET}_{what}_v1.html")
        (svg,) = out(root).glob(f"*_{DATASET}_{what}_v1_white.svg")
        titles = {axis: json.loads(f'"{found}"').split("<br>") for axis, found in re.findall(
            r'"([xy])axis":\{[^{}]*?"title":\{"text":"(.*?)"', html.read_text())}  # fmt: skip
        assert max(len(line) for line in titles["x"]) <= SCREEN_TITLE_CHARS, (what, titles)
        assert titles["y"] == _y_label(svg), what
