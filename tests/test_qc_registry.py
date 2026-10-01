"""conf/qc_metrics.toml and ctx.qc.registry: the bands, their sources, and MultiQC's colours."""

import json
import re
import shutil
from pathlib import Path

import pytest

from ctx.qc import multiqc_io, registry

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "conf" / "qc_metrics.toml"
FIXTURE = ROOT / "tests" / "fixtures" / "rnaseq"


@pytest.fixture(scope="module")
def reg():
    return registry.load(REGISTRY)


def test_the_repo_registry_loads(reg):
    assert len(reg.metrics) >= 30
    assert {m.stage for m in reg.metrics.values()} == set(registry.STAGES)
    for m in reg.metrics.values():
        assert m.source and m.source_type in registry.SOURCE_TYPES, m.id
        if m.source_type == "course_convention" or (m.kind == "numeric" and not m.banded):
            assert m.rationale, m.id


def test_nothing_is_excluded_before_the_record(reg):
    """P2M M2.8: flags only until the Step-1 record, so no exclusion can change it."""
    if not list((ROOT / "deliverables").glob("*_de_tier1_sanity_v*.csv")):
        assert not [m.id for m in reg.metrics.values() if m.action == "exclude"]


def test_every_problem_is_reported_together(tmp_path):
    bad = REGISTRY.read_text()
    bad = bad.replace('source_type = "tool_code"', 'source_type = "hearsay"', 1)
    bad = bad.replace("good = 30.0\nusable = 50.0", "good = 30.0", 1)
    bad = bad.replace('stage = "counts"', 'stage = "later"', 1)
    (tmp_path / "r.toml").write_text(bad)
    with pytest.raises(registry.RegistryError) as raised:
        registry.load(tmp_path / "r.toml")
    message = str(raised.value)
    assert "'hearsay'" in message and "needs both good and usable" in message
    assert "stage = 'later'" in message


@pytest.mark.parametrize("direction", ["higher", "lower"])
def test_bands_are_closed_on_the_good_side(reg, direction):
    metric = next(m for m in reg.metrics.values() if m.banded and m.direction == direction)
    assert metric.band(metric.good) == "good"
    step = 1 if direction == "higher" else -1
    assert metric.band(metric.good - step * 1e-9) in {"usable", "bad"}
    assert metric.band(None) is None


def multiqc_class(rules: dict, value: float) -> str | None:
    """MultiQC 1.35's rule: classes in the order pass, warn, fail; the last that matches wins."""
    ops = {"gt": float.__gt__, "lt": float.__lt__, "eq": float.__eq__}
    found = None
    for level in ("pass", "warn", "fail"):
        if any(ops[k](value, v) for cond in rules.get(level, []) for k, v in cond.items()):
            found = level
    return found


def test_multiqc_rules_colour_as_the_registry_bands(reg):
    names = {"pass": "good", "warn": "usable", "fail": "bad"}
    rules = registry.multiqc_config(reg, "")["table_cond_formatting_rules"]
    assert rules
    for column, rule in rules.items():
        metric = next(m for m in reg.metrics.values() if m.column == column)
        edges = {metric.good, metric.usable}
        for value in sorted(edges | {e + d for e in edges for d in (-0.5, 0.5)} | {0.0, 1e9}):
            assert names[multiqc_class(rule, value)] == metric.band(value), (column, value)


@pytest.mark.parametrize("metric_id", ["salmon_strand_bias", "gene_strand_bias"])
def test_strandedness_follows_nf_core(reg, metric_id):
    """Mutation caught: gene_strand_bias's unstranded 0.1 set to 0.01 (0.52 undetermined)."""
    metric = reg.metrics[metric_id]
    assert [metric.strandedness(b) for b in (0.9, 0.1, 0.5, 0.52, 0.6)] == [
        "forward", "reverse", "unstranded", "unstranded", "undetermined",
    ]  # fmt: skip


@pytest.mark.parametrize("metric_id", ["salmon_strand_bias", "gene_strand_bias"])
@pytest.mark.parametrize(
    ("share", "call"),
    [(0.8, "forward"), (0.79, "undetermined"), (0.2, "reverse"), (0.21, "undetermined"),
     (0.549, "unstranded"), (0.55, "undetermined")],
)  # fmt: skip
def test_strand_cut_offs_are_nf_cores(reg, metric_id, share, call):
    """Each side of stranded 0.8 and of unstranded 0.1, as TODO.md M3.2 sets them.

    Mutations caught: gene_strand_bias's stranded at 0.85 (0.8 undetermined) or 0.75 (0.79
    forward)."""
    assert reg.metrics[metric_id].strandedness(share) == call


def test_salmon_strand_bias_is_computed_as_stated(reg):
    """Salmon's strand_mapping_bias is SF / (SF + SR) for single-end runs and ISF / (ISF + ISR)
    for paired-end ones, in every fixture report; the registry's text says both.

    Mutation caught: the old text, SF / (SF + SR) alone."""
    computed = reg.metrics["salmon_strand_bias"].computed
    assert "SF / (SF + SR)" in computed and "ISF / (ISF + ISR)" in computed
    reports = sorted((FIXTURE / "salmon").glob("*/lib_format_counts.json"))
    assert reports
    for path in reports:
        counts = json.loads(path.read_text())
        paired = counts["expected_format"].startswith("I")
        forward, reverse = ("ISF", "ISR") if paired else ("SF", "SR")
        share = counts[forward] / (counts[forward] + counts[reverse])
        assert share == pytest.approx(counts["strand_mapping_bias"], abs=1e-12), path.parent.name


def test_reference_start_codons_judge_the_measured_shares(reg):
    """GENCODE 50 on chr22 is good (8,852 of 8,866), a genome shifted by one base is bad.

    Mutations caught: usable at 0.96 (0.95 turns bad); good at 0.995 (0.99 turns usable)."""
    metric = reg.metrics["reference_start_codons"]
    assert (metric.stage, metric.table, metric.column) == ("align", "reference", "atg_share")
    assert metric.band(8852 / 8866) == "good"
    assert metric.band(0.95) == "usable"
    assert metric.band(0 / 8866) == "bad"
    # Each side of good 0.99 and usable 0.9, as TODO.md M3.2 sets them.
    assert [metric.band(v) for v in (0.99, 0.989, 0.9, 0.899)] == ["good", "usable", "usable",
                                                                   "bad"]  # fmt: skip
    assert metric.action == "exclude"


def test_multiqc_reads_the_generated_config(reg, tmp_path):
    """One real run: every trim and quant column the registry names exists, and MultiQC paints
    the fixture's Trim Galore totals, all at least 10,000 reads, in the good band's colour."""
    inputs = tmp_path / "inputs"
    for part in ("trimgalore", "salmon"):
        shutil.copytree(FIXTURE / part, inputs / part)
    config = tmp_path / "config.json"
    config.write_text(json.dumps(registry.multiqc_config(reg, "generated by ctx qc")))
    run = multiqc_io.run([inputs], tmp_path / "multiqc", config=config,
                         modules=("trim_galore", "salmon"))  # fmt: skip
    general = multiqc_io.table(run, "general_stats")
    for m in reg.metrics.values():
        if m.stage in {"trim", "quant"}:
            assert m.column in general.columns, m.id
    html = (tmp_path / "multiqc" / "multiqc_report.html").read_text()
    cells = re.findall(r'class="data-coloured trim_galore-tg_total_reads ">.*?background-color:'
                       r"(#[0-9A-F]{6})", html, re.S)  # fmt: skip
    assert cells and set(cells) == {"#009E73"}
    assert "generated by ctx qc" in html


def test_spans_are_the_ranges_band_assigns(reg):
    """A value inside a span gets that span's band, at the edges too; one cut-off, two bands."""
    banded = [m for m in reg.metrics.values() if m.banded]
    assert banded
    for m in banded:
        spans = m.spans()
        edges = sorted({v for _, lo, hi in spans for v in (lo, hi) if v is not None})
        probes = [v + d for v in edges for d in (-1e-6, 0.0, 1e-6)]
        for value in probes:
            inside = [b for b, lo, hi in spans if (lo is None or value >= lo)
                      and (hi is None or value <= hi)]  # fmt: skip
            assert m.band(value) in inside, (m.id, value, inside)
    q30 = reg.metrics["q30_rate"]
    assert q30.good == q30.usable
    assert [b for b, *_ in q30.spans()] == ["good", "bad"]
