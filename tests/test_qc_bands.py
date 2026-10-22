"""D-052's QC bands (P2Q.9), P2Q.10's definitions and ctx.qc.identity: each test names the
mutation it was seen to catch."""

import json
import tomllib
from pathlib import Path

import polars as pl
import pytest

from ctx.qc import checkpoint, fastqc, figures, identity, multiqc_io, registry

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "conf" / "qc_metrics.toml"
PROPOSAL = ROOT / "docs" / "plans" / "2026-10-01_qc_bands_proposal.json"
REG = registry.load(REGISTRY)


def with_metric(tmp_path: Path, key: str, **fields: object) -> Path:
    """The repo registry with fields of one metric replaced or added, as a new file."""
    raw = tomllib.loads(REGISTRY.read_text())
    raw["metrics"][key].update(fields)
    lines = []
    for name, entry in raw["metrics"].items():
        lines.append(f"[metrics.{name}]")
        lines += [f"{k} = {json.dumps(v)}" for k, v in entry.items()]
    head = "\n".join(f"{k} = {json.dumps(v)}" for k, v in raw["registry"].items())
    path = tmp_path / "r.toml"
    path.write_text(f"[registry]\n{head}\n" + "\n".join(lines) + "\n")
    return path


# Fields the registry words otherwise than the proposal, each with the text it keeps instead.
DEVIATIONS = {
    # The proposal claimed a control the page lacks: "the page lets you move" the fastp grid.
    ("fastp_surviving", "rationale"): "the page maps one run's -q by -l, without a control.",
}


def test_the_approved_proposal_is_applied_field_by_field():
    """Every add and replace of the approved proposal is in the registry, every removal gone.

    Mutation caught: fastqc_adapter's plot_bands set to [5.0, 20.0]."""
    metrics = tomllib.loads(REGISTRY.read_text())["metrics"]
    for entry in json.loads(PROPOSAL.read_text()):
        metric = metrics[entry["id"]]
        for field, value in {**entry["registry_add"], **entry["registry_replace"]}.items():
            if (kept := DEVIATIONS.get((entry["id"], field))) is not None:
                assert kept in metric[field] and value != metric[field], (entry["id"], field)
                continue
            assert metric.get(field) == value, (entry["id"], field)
        assert not [f for f in entry["registry_remove"] if f in metric], entry["id"]


@pytest.mark.parametrize(("key", "bands"), [("fastqc_adapter", [10.0, 5.0]),
                                            ("globin_fraction", [0.01, 0.001]),
                                            ("fastp_surviving", [90.0, 95.0])])  # fmt: skip
def test_a_band_pair_against_its_direction_is_refused(tmp_path, key, bands):
    """lower wants good below usable, higher above it; direction none takes either order.

    Mutation caught: the order check dropped (both lower and higher swaps load)."""
    raw = tomllib.loads(REGISTRY.read_text())["metrics"][key]
    path = with_metric(tmp_path, key, plot_bands=bands)
    if raw["direction"] == "none":
        assert registry.load(path).metrics[key].plot_spans()[0][0] == "good"
        return
    with pytest.raises(registry.RegistryError, match=rf"{key}\] plot_bands .* good must lie"):
        registry.load(path)


def test_direction_none_reads_the_side_from_the_pair():
    """FastQC's adapter limits, direction none: low is good. Quality's, high is good.

    Mutation caught: plot_spans always shading the high side good."""
    adapter = REG.metrics["fastqc_adapter"].plot_spans()
    assert adapter == [("good", None, 5.0), ("usable", 5.0, 10.0), ("bad", 10.0, None)]
    quality = REG.metrics["fastqc_per_base_quality"].plot_spans()
    assert quality == [("good", 28.0, None), ("usable", 20.0, 28.0), ("bad", None, 20.0)]


@pytest.mark.parametrize(("fields", "message"), [
    ({"plot_lines": 50.0}, "unknown key 'plot_lines'"),
    ({"plot_source": ""}, "plot_line needs a plot_source"),
    ({"plot_unstranded_tolerance": 0.01}, "plot_unstranded_tolerance needs a strandedness"),
    ({"plot_expected_bias_positions": 13}, "plot_expected_bias_positions needs a status"),
])  # fmt: skip
def test_unknown_or_unsourced_plot_keys_are_refused(tmp_path, fields, message):
    """A misspelt key once loaded and drew nothing; a line needs its source.

    Mutations caught: the unknown-key check dropped; the plot_source check dropped."""
    with pytest.raises(registry.RegistryError, match=message):
        registry.load(with_metric(tmp_path, "read_length", **fields))


def strip_of(metric_id: str, values: list[float]):
    """The strip of one metric on made-up samples, one stratum, nothing flagged."""
    metric = REG.metrics[metric_id]
    names = [f"s{i}" for i in range(len(values))]
    frame = pl.DataFrame({"sample": names, metric.column: values})
    units = pl.DataFrame({"sample": names, "cell_type": ["x"] * len(values)})
    judged = (pl.DataFrame(schema=checkpoint.BAND_COLUMNS),
              pl.DataFrame(schema=checkpoint.FLAG_COLUMNS))  # fmt: skip
    return figures.strip(metric, frame, judged, units, ["cell_type"], "sample")


def test_strips_draw_display_bands_lines_and_call_regions():
    """plot_bands shade, plot_line is a dashed line, a strand metric shows its call regions.

    Mutations caught: strip() drawing Metric.spans() alone (no shading on fastp_surviving);
    plot_line left undrawn; the tolerance strip left out."""
    fastp = strip_of("fastp_surviving", [99.0, 98.5, 97.0])
    assert [(s.band, s.low, s.high) for s in fastp.spans] == [
        ("good", 95.0, None),
        ("usable", 90.0, 95.0),
        ("bad", None, 90.0),
    ]
    gc = strip_of("percent_gc", [48.0, 50.0, 51.0])
    assert [v for _, v in gc.references] == [49.8] and not gc.spans
    strand = strip_of("salmon_strand_bias", [0.5, 0.51, 0.49])
    edges = [v for s in strand.spans for v in (s.low, s.high)]
    assert edges == pytest.approx([0.0, 0.2, 0.45, 0.55, 0.8, 1.0, 0.49, 0.51])
    assert figures.shading_kind(REG.metrics["mt_rrna_fraction"]) == "line"


def test_a_reason_never_prints_a_value_equal_to_its_bound():
    """At 3 digits 99.96 and 99.98 both read 100; the reason takes more digits.

    Mutation caught: shown_apart returning 3 digits whatever they print."""
    assert checkpoint.shown_apart(99.96, 99.98) == ("99.96", "99.98")
    assert checkpoint.shown_apart(0.5, 0.25) == ("0.5", "0.25")
    assert checkpoint.shown_apart(20_000_000.4, 20_000_000.1) == ("20,000,000.4", "20,000,000.1")


ANTISENSE = REG.metrics["top_gene_share"].exclude_gene_ids
SCHEMA = ["gene_id", "gene_name", "gene_type", "chrom"]


def test_top_gene_and_rrna_follow_p2q10():
    """The 18S-antisense lncRNA and Mt_rRNA are not the top gene; the lncRNA counts as rRNA.

    Mutations caught: the antisense gene left in the top-gene search (top 0.6); it left out of
    rrna_fraction (0.03); the version kept when matching (neither excluded)."""
    genes = [("A", "MT-RNR2", "Mt_rRNA", "chrM"), ("B", "RNA5-8SN1", "rRNA", "chr21"),
             ("ENSG00000280614.1", "ENSG00000280614", "lncRNA", "chr21"),
             ("E", "ACTB", "protein_coding", "chr7")]  # fmt: skip
    types = pl.DataFrame(genes, schema=SCHEMA, orient="row")
    counts = pl.DataFrame({"gene_id": [g[0] for g in genes], "s1": [300, 30, 600, 70]})
    row = checkpoint.count_metrics(counts, types, ANTISENSE).row(0, named=True)
    assert row["top_gene_share"] == pytest.approx(0.07)
    assert row["rrna_fraction"] == pytest.approx(0.63)
    assert row["mt_rrna_fraction"] == pytest.approx(0.3)


def test_without_gene_types_no_top_gene_share_is_measured():
    """Without an annotation no Mt_rRNA gene can be told apart, so the metric is left out.

    Fails before P2Q.10's correction, which reported MT-RNR2's 0.3 under the same name.
    Mutation caught: the untyped branch searching every gene again (0.3 returned)."""
    counts = pl.DataFrame({"gene_id": ["A", "E"], "s1": [300, 70]})
    untyped = checkpoint.count_metrics(counts, None, ANTISENSE)
    assert untyped.columns == ["sample", "library_size", "genes_detected"]


def planted_annotation() -> tuple[pl.DataFrame, pl.DataFrame]:
    """A fixture annotation with every class P2Q.10 names, and one sample's counts.

    Mt_rRNA is planted under an Ensembl ID and a name with no MT- prefix, so only its type
    marks it; the two 18S-antisense lncRNAs carry versions; an 18S and a 28S gene, as an
    annotation with a 45S unit would type them (GENCODE 50 has none), sit on a nuclear contig;
    a rRNA-typed gene on chrM stays out of rrna_fraction.
    """
    genes = [
        ("ENSG00000211459.2", "RNR_PLANTED", "Mt_rRNA", "chrM"),
        ("ENSG00000280614.1", "ENSG00000280614", "lncRNA", "chr21"),
        ("ENSG00000281383.1", "ENSG00000281383", "lncRNA", "chr21"),
        ("G18S.1", "RNA18SN5", "rRNA", "GL000220.1"),
        ("G28S.1", "RNA28SN5", "rRNA", "GL000220.1"),
        ("G5S.1", "RNA5S1", "rRNA_pseudogene", "chr1"),
        ("GMT.1", "MT-ODD", "rRNA", "chrM"),
        ("GTOP.1", "DEFA1B", "protein_coding", "chr8"),
        ("GLOW.1", "ACTB", "protein_coding", "chr7"),
    ]
    counts = [500, 400, 350, 200, 100, 10, 20, 150, 270]
    types = pl.DataFrame(genes, schema=SCHEMA, orient="row")
    return types, pl.DataFrame({"gene_id": [g[0] for g in genes], "s1": counts})


def test_the_planted_fixture_measures_each_corrected_definition():
    """On the planted annotation: top gene ACTB, 270 of 2,000; nuclear rRNA 18S + 28S + 5S +
    the two lncRNAs, 1,060 of 2,000; Mt_rRNA 500 of 2,000.

    Mutations caught: Mt_rRNA left in the top-gene search (0.25); either lncRNA left in it (0.2,
    0.175); the 18S and 28S genes not counted as rRNA (0.38); the chrM rRNA gene counted as
    nuclear rRNA (0.54)."""
    types, counts = planted_annotation()
    row = checkpoint.count_metrics(counts, types, ANTISENSE).row(0, named=True)
    assert row["library_size"] == 2000
    assert row["top_gene_share"] == pytest.approx(270 / 2000)
    assert row["rrna_fraction"] == pytest.approx(1060 / 2000)
    assert row["mt_rrna_fraction"] == pytest.approx(500 / 2000)


def test_the_exclusion_list_is_read_from_the_registry_alone(tmp_path):
    """The two IDs live in conf/qc_metrics.toml, sourced; checkpoint holds no copy of them.

    Fails before: the registry had no exclude_gene_ids and checkpoint a constant of its own.
    Mutations caught: a constant reinstated in checkpoint.py; a list without a source; a
    versioned ID, which would never match an unversioned gene."""
    metric = REG.metrics["top_gene_share"]
    assert metric.exclude_gene_ids == ("ENSG00000280614", "ENSG00000281383")
    assert "TODO.md P2Q.10" in metric.exclude_source
    code = (ROOT / "src" / "ctx" / "qc" / "checkpoint.py").read_text()
    assert "ENSG00000280614" not in code and "ENSG00000281383" not in code
    types, counts = planted_annotation()
    one = with_metric(tmp_path, "top_gene_share", exclude_gene_ids=["ENSG00000280614"])
    listed = registry.load(one).metrics["top_gene_share"].exclude_gene_ids
    row = checkpoint.count_metrics(counts, types, listed).row(0, named=True)
    assert row["top_gene_share"] == pytest.approx(350 / 2000)
    shorter = checkpoint.count_metrics(counts, types, ()).row(0, named=True)
    assert shorter["top_gene_share"] == pytest.approx(400 / 2000)
    unsourced = tomllib.loads(REGISTRY.read_text())["metrics"]["top_gene_share"]
    unsourced.pop("exclude_source")
    with pytest.raises(registry.RegistryError, match="needs an exclude_source"):
        registry.load(with_entry(tmp_path, "top_gene_share", unsourced))
    with pytest.raises(registry.RegistryError, match="unversioned Ensembl gene IDs"):
        registry.load(with_metric(tmp_path, "top_gene_share",
                                  exclude_gene_ids=["ENSG00000280614.1"]))  # fmt: skip


def with_entry(tmp_path: Path, key: str, entry: dict) -> Path:
    """The repo registry with one metric's table replaced whole, as a new file."""
    raw = tomllib.loads(REGISTRY.read_text())
    raw["metrics"][key] = entry
    lines = []
    for name, fields in raw["metrics"].items():
        lines.append(f"[metrics.{name}]")
        lines += [f"{k} = {json.dumps(v)}" for k, v in fields.items()]
    head = "\n".join(f"{k} = {json.dumps(v)}" for k, v in raw["registry"].items())
    path = tmp_path / "whole.toml"
    path.write_text(f"[registry]\n{head}\n" + "\n".join(lines) + "\n")
    return path


def test_the_sex_check_has_one_definition():
    """M5 imports M2's sex check (CONVENTIONS: it moves to ctx.qc.identity after the record).

    Mutation caught: a copy of sex_check left in ctx.de.run."""
    from ctx.de import run  # noqa: PLC0415 - heavy

    assert run.sex_check is identity.sex_check
    assert "def sex_check" not in (ROOT / "src" / "ctx" / "de" / "run.py").read_text()


def multiqc_data(tmp_path: Path) -> multiqc_io.MultiqcRun:
    """A MultiQC data folder written by hand, in MultiQC 1.35's plot-data layout."""
    data = tmp_path / "multiqc_data"
    data.mkdir()
    (data / "fastqc_adapter_content_plot.txt").write_text(
        "Sample\t0\t1\n"
        "r1_1 - illumina\t(1, 0.0)\t(2, 3.0)\n"
        "r1_1 - polya\t(1, 1.0)\t(2, 2.0)\n"
        "r1_2 - illumina\t(1, 9.0)\t(2, 9.0)\n")  # fmt: skip
    (data / "fastqc_per_sequence_gc_content_plot_Percentages.txt").write_text(
        "Sample\t0\t1\nr1_1\t(40.0, 1.0)\t(50.0, 2.0)\nr1_2\t(40.0, 5.0)\t(50.0, 5.0)\n"
        "Theoretical GC Content\t(40.0, 0.5)\t(50.0, 3.0)\n")  # fmt: skip
    (data / "fastqc_overrepresented_sequences_plot.txt").write_text(
        "Sample\tTop overrepresented sequence\tSum of remaining overrepresented sequences\n"
        "r1_1\t0.2\t1.0\nr1_2\t0.9\t0.0\n")  # fmt: skip
    return multiqc_io.MultiqcRun(data, "1.35")


def test_fastqc_plots_take_read_1_and_the_largest_adapter(tmp_path):
    """Read 2 is never drawn; the adapter curve is the largest adapter at each position, as
    FastQC judges; MultiQC's theoretical row is the reference, not a run.

    Mutations caught: the first adapter row kept instead of the largest (1 at position 1
    reads 0); read 2 drawn; the theoretical row mapped as a run (QcError)."""
    plots = fastqc.read(multiqc_data(tmp_path), {"r1_1": "R1", "r1_2": "R1"}, ["r1_1"],
                        "run_accession")  # fmt: skip
    adapter = plots.curves["fastqc_adapter"]
    assert adapter.rows() == [("R1", 1.0, 1.0), ("R1", 2.0, 3.0)]
    assert plots.curves["fastqc_gc_distribution"]["y"].to_list() == [1.0, 2.0]
    assert plots.theoretical is not None and plots.theoretical["y"].to_list() == [0.5, 3.0]
    assert plots.overrepresented.rows() == [("R1", 0.2)]


def test_read_length_groups_come_from_the_data():
    """The declared batch where the dataset file declares one, else the run's own length.

    Mutation caught: groups from a fixed list of Step-1's lengths (50 and 150 nt unknown)."""
    sheet = pl.DataFrame({"sample": ["a", "b", "c"], "batch": ["short", "short", "long"],
                          "read_length": [50.0, 52.0, 150.0]})  # fmt: skip
    own = checkpoint.length_groups(sheet, declared=False)[checkpoint.LENGTH_GROUP].to_list()
    assert own == ["50 nt", "52 nt", "150 nt"]
    declared = checkpoint.length_groups(sheet, declared=True)[checkpoint.LENGTH_GROUP]
    assert declared.to_list() == ["short", "short", "long"]
    alone = checkpoint.length_groups(sheet.drop("batch", "read_length"), declared=True)
    assert alone[checkpoint.LENGTH_GROUP].to_list() == ["all"] * 3


def test_read_length_is_judged_within_its_group_only():
    """Two batches mixed by design: across them the MAD rule flags the whole long batch;
    within them it flags only the run unlike its own batch (Rytis 2026-10-09).

    Mutation caught: read_length judged across the --strata (the five long runs flagged)."""
    lengths = {"l1": 100, "l2": 101, "l3": 102, "l4": 101, "l5": 100,
               "s1": 75, "s2": 76, "s3": 77, "s4": 76, "s5": 75, "s6": 76, "odd": 60}  # fmt: skip
    names = list(lengths)
    column = REG.metrics["read_length"].column
    values = pl.DataFrame({"sample": names, column: [float(v) for v in lengths.values()]})
    batch = ["long" if n.startswith("l") else "short" for n in names]
    samples = pl.DataFrame({"sample": names, "cell_type": ["x"] * len(names), "batch": batch})
    samples = checkpoint.length_groups(samples, declared=True)
    _, flags = checkpoint.judge(values, samples, REG, ["cell_type"])
    assert flags.filter(pl.col("metric") == "read_length")["sample"].to_list() == ["odd"]
    assert REG.metrics["read_length"].judged_within == checkpoint.LENGTH_GROUP
