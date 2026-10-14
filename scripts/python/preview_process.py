#!/usr/bin/env python3
# /// codebase-meta
# exempt: preview builder; draws the M3 page on made-up data where no reads or genome exist
# ///
"""The course preview of M3 Trim, align, count (``ctx process``), drawn on made-up data.

``ctx process`` needs reads and a genome, which the preview's made-up dataset does not have. This
page draws what the real page draws, through the real page's own builders and words
(``ctx.process.page``, ``content/process.toml``): the aligner's outcome per sample and for one
featured sample under each setting, reads assigned at each featureCounts ``-s``, and plain
against ``-M`` gene counts as a scatter and a Bland-Altman plot. Every count comes from a
seeded generator shaped like STAR's Log.final.out and featureCounts' summaries; the statistics
on those counts are computed. Samples, groups and the plain gene counts are the dataset's gather
outputs, so this page agrees with the other preview pages. ``build_preview.py`` calls ``page``.
Every caption starts "Made-up data.", and the build refuses a figure file that does not say so,
since a figure is downloaded and shown apart from the page (D-050).
"""

from __future__ import annotations

import re
import time
from datetime import date
from pathlib import Path

import numpy as np
import polars as pl
import typer

from ctx import runtime
from ctx.cli import app
from ctx.de.run import Writer, newest
from ctx.io.naming import MODULES, stem
from ctx.process import CATEGORIES, Aligner, figures
from ctx.process import page as pages
from ctx.process import tables as tab
from ctx.qc import registry as reg
from ctx.qc.checkpoint import gene_types
from ctx.report import text
from ctx.report.builder import ModulePage, Reproduction, write
from ctx.report.cards import parameter_cards

# The made-up dataset's name on the page, one word for every preview (content/ui.toml).
SHOWN = ("ui", "preview.dataset")

PAGE = MODULE = "process"
SEED = 3
# The planned default (--aligner star); its outcomes are STAR's Log.final.out rows.
ALIGNER = Aligner.STAR
GENOME = "chr22"
# GRCh38 chr22's length in bases (D-037, measured on the real chr22 FASTA).
GENOME_LENGTH = 50_818_468
STRAND = ("s0", "s1", "s2")
# The annotation build_preview writes beside the gather outputs; without it no gene has a type.
ANNOTATION = Path("resources") / "gencode" / "mini.gtf"
# What every preview figure file must carry, in its caption.
MADE_UP = re.compile(r"Made-up", re.I)


def made_up_runs(
    rng: np.random.Generator, sheet: pl.DataFrame, assigned: np.ndarray
) -> pl.DataFrame:
    """Made-up Log.final.out shares and featureCounts strand shares, one row per sample.

    Each batch shifts the unique share; every library is reverse-stranded. A sample's gather
    counts are taken as its reads assigned at ``-s 2``, which fixes its input reads.
    """
    n = sheet.height
    levels = sheet["batch"].unique(maintain_order=True).to_list()
    shift = dict(zip(levels, rng.normal(0.0, 1.5, len(levels)), strict=True))
    batch = np.array([shift[b] for b in sheet["batch"]])
    shares = {
        "multi": np.clip(rng.normal(5.0, 0.6, n) - batch / 2, 1.0, None),
        "too_many": rng.uniform(0.1, 0.4, n),
        "mismatches": rng.uniform(0.05, 0.3, n),
        "too_short": np.clip(rng.normal(5.0, 1.2, n) - batch / 2, 1.0, None),
        "other": rng.uniform(0.3, 1.0, n),
    }
    unique = 100.0 - sum(shares.values())
    genic = unique / 100 * rng.uniform(0.74, 0.82, n)  # unique reads that fall on a gene
    forward = rng.uniform(0.04, 0.08, n)
    strand = {"s0": genic * (1 - rng.uniform(0.01, 0.02, n)), "s1": genic * forward,
              "s2": genic * (1 - forward)}  # fmt: skip
    return pl.DataFrame({
        "sample": sheet["sample"],
        "input_reads": np.round(assigned / strand["s2"]).astype(np.int64),
        "unique": unique, **shares, **strand,
    })  # fmt: skip


def outcome_reads(shares: dict[str, float], reads: int) -> dict[str, int]:
    """Reads per STAR outcome from percent shares, unique taking the rounding so they add up."""
    rest = {k: round(shares[k] / 100 * reads) for k in CATEGORIES[ALIGNER][1:]}
    return {"unique": reads - sum(rest.values()), **rest, "input": reads}


def settings_shares(run: dict[str, float]) -> dict[str, dict[str, float]]:
    """Made-up outcome shares of one sample as run, on half the genome, untrimmed, and at limit 1.

    Limit 1 is derived as the real page derives it: multimapped reads become too many loci.
    """
    main = {o: float(run[o]) for o in CATEGORIES[ALIGNER]}
    # About half the reads come from genes on the removed half; a few align to paralogues.
    cut = main | {"unique": main["unique"] * 0.52, "multi": main["multi"] * 0.6,
                  "mismatches": main["mismatches"] + 0.4, "other": main["other"] + 0.5}  # fmt: skip
    cut["too_short"] = 100.0 - sum(v for k, v in cut.items() if k != "too_short")
    # Adapter tails that fastp would cut leave reads too short to align, or mismatched.
    untrimmed = main | {"unique": main["unique"] - 1.4, "mismatches": main["mismatches"] + 0.3,
                        "too_short": main["too_short"] + 1.1}  # fmt: skip
    limit = main | {"multi": 0.0, "too_many": main["too_many"] + main["multi"]}
    return {"main": main, "limit_1": limit, "cut": cut, "untrimmed": untrimmed}


def with_multimappers(rng: np.random.Generator, plain: np.ndarray) -> np.ndarray:
    """Made-up ``-M`` counts: the plain counts plus multimapped reads, most on gene families."""
    genes = plain.shape[0]
    # Per gene, extra reads per plain read: none, a few, or a gene family's many.
    rates = np.stack([np.zeros(genes), rng.lognormal(np.log(0.01), 0.8, genes),
                      rng.lognormal(np.log(0.3), 0.8, genes)])  # fmt: skip
    kind = rng.choice(len(rates), size=genes, p=(0.45, 0.45, 0.10))
    extra: np.ndarray = rng.poisson(plain * rates[kind, np.arange(genes)][:, None])
    return np.asarray(plain + extra)


def strand_table(runs: pl.DataFrame, outcomes: dict[str, dict[str, int]]) -> pl.DataFrame:
    """The strandedness table's columns the page reads, from the made-up shares."""
    metric = reg.load(reg.resolve(None)).metrics["gene_strand_bias"]
    rows = []
    for r in runs.iter_rows(named=True):
        n, unique = r["input_reads"], outcomes[r["sample"]]["unique"]
        assigned = {s: round(r[s] * n) for s in STRAND}
        share = assigned["s1"] / (assigned["s1"] + assigned["s2"])
        rows.append({"sample": r["sample"],
                     **{f"assigned_{s}": float(v) for s, v in assigned.items()},
                     **{f"no_feature_{s}": float(unique - v) for s, v in assigned.items()},
                     **{f"ambiguous_{s}": 0.0 for s in STRAND}, "forward_share": share,
                     "call": metric.strandedness(share), "compare_bias": None,
                     "compare_call": None})  # fmt: skip
    return pl.DataFrame(rows).with_columns(used=pl.col("call"))


def gene_table(root: Path, dataset: str, counts: pl.DataFrame) -> pl.DataFrame:
    """gene_id, gene_name, gene_type: the made-up annotation's types where it was written."""
    genes = newest(root / "outputs" / "01_gather", dataset, "genes")
    names = pl.read_parquet(genes)
    if (root / ANNOTATION).is_file():
        types = gene_types(root / ANNOTATION).select("gene_id", "gene_type")
        names = names.join(types, on="gene_id", how="left")
    else:
        names = names.with_columns(gene_type=pl.lit(None, pl.String))
    return counts.select("gene_id").join(names.select("gene_id", "gene_name", "gene_type"),
                                         on="gene_id", how="left")  # fmt: skip


def frames(root: Path, dataset: str, sheet: pl.DataFrame, counts: pl.DataFrame) -> pages.Frames:
    """The made-up run as the real page's frames: outcomes, settings, strand, one comparison."""
    samples = sheet["sample"].to_list()
    plain = counts.select(samples).to_numpy().astype(np.int64)
    rng = np.random.default_rng(SEED)
    runs = made_up_runs(rng, sheet, plain.sum(axis=0))
    multi = with_multimappers(rng, plain)
    outcomes = {r["sample"]: outcome_reads(r, r["input_reads"]) for r in runs.iter_rows(named=True)}
    featured = runs.row(0, named=True)
    settings = {k: outcome_reads(v, featured["input_reads"])
                for k, v in settings_shares(featured).items()}  # fmt: skip
    genes = gene_table(root, dataset, counts)
    long = pl.DataFrame({
        "comparison": "multimap", "sample": np.repeat(samples, plain.shape[0]),
        "gene_id": np.tile(counts["gene_id"].to_numpy(), len(samples)),
        "a": plain.T.ravel().astype(float), "b": multi.T.ravel().astype(float),
    })  # fmt: skip
    alignment = alignment_rows(rng, outcomes)
    called = counts.select("gene_id", *samples)
    return pages.Frames(
        alignment=alignment, outcomes=outcomes, settings=settings,
        strand=strand_table(runs, outcomes), counting=pl.DataFrame(),
        gene_counts=genes.join(called, on="gene_id"), comparison=long,
        agreements=tab.agreements(long, genes), reference=reference_row(genes.height),
        incomplete=None, context=None,
    )  # fmt: skip


def alignment_rows(rng: np.random.Generator, outcomes: dict[str, dict[str, int]]
                   ) -> pl.DataFrame:  # fmt: skip
    """The alignment table's columns the page prints, made up as STAR would report them.

    Splices per 100 uniquely mapped reads, their annotated share and non-canonical count;
    flagstat's records, a multimapper once per alignment, beside primary mapped reads.
    """
    rows = []
    for sample, r in outcomes.items():
        unique, multi, total = r["unique"], r["multi"], r["input"]
        splices = round(unique * rng.uniform(0.18, 0.30))
        records = total + round(multi * rng.uniform(1.5, 3.0))  # secondary alignments
        rows.append({
            "sample": sample, "input_reads": total, "unique_pct": 100 * unique / total,
            "splices_per_100_unique": 100 * splices / unique,
            "splices_annotated_share": rng.uniform(0.97, 0.995),
            "splices_noncanonical": round(splices * rng.uniform(0.001, 0.004)),
            "records": records, "mapped_records_pct": 100 * (records - total + unique + multi)
            / records, "primary_mapped_pct": 100 * (unique + multi) / total,
        })  # fmt: skip
    return pl.DataFrame(rows)


def reference_row(genes: int) -> pl.DataFrame:
    """The reference table's row the page prints: chr22 by name, the made-up dataset's genes."""
    return pl.DataFrame([{
        "reference": text.word(PAGE, "references.main"), "contigs": GENOME,
        "length_nt": GENOME_LENGTH, "genes_subset": genes, "genes_genome": genes,
        "atg_share": None,
    }])  # fmt: skip


def check_made_up(files: list[Path]) -> None:
    """Refuse a preview figure whose SVG does not say its data are made up.

    Raises:
        RuntimeError: Naming each figure that does not.
    """
    if bare := [f.name for f in files if f.suffix == ".svg" and not MADE_UP.search(f.read_text())]:
        raise RuntimeError(f"preview figures without the made-up mark: {bare}")


def page(root: Path, dataset: str, day: date) -> Path:
    """Write the M3 preview page and its figures under ``root/outputs/03_process``.

    Args:
        root: A ctx root holding the dataset's gather outputs.
        dataset: The dataset slug.
        day: The date every output name carries.

    Returns:
        The page.
    """
    started = time.perf_counter()
    gathered = root / "outputs" / "01_gather"
    inputs = [newest(gathered, dataset, what) for what in ("samplesheet", "counts")]
    sheet = pl.read_parquet(inputs[0]).unique("sample", keep="first", maintain_order=True)
    counts = pl.read_parquet(inputs[1])
    found = frames(root, dataset, sheet, counts)
    releases = {"input:samplesheet": inputs[0].name, "input:counts": inputs[1].name,
                "made_up": f"numpy default_rng({SEED})"}  # fmt: skip
    writer = Writer(root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}", dataset, day,
                    releases)  # fmt: skip
    writer.table(
        found.agreements.concordance,
        "concordance",
        pages.described("tables.concordance", found.agreements.concordance.columns),
    )
    # chr22 alone of a whole annotation: a partial reference, whatever the made-up genes say.
    shown = pages.Shown(ALIGNER, sheet["sample"][0], pages.groups_of(sheet, sheet["sample"]),
                        GENOME, reg.load(reg.resolve(None)).metrics["gene_strand_bias"],
                        found.gene_counts.height, partial=True,
                        prefix=f"{text.word(PAGE, 'preview.caption')} ")  # fmt: skip
    body = pages.results(writer, found, shown, [text.word(PAGE, "preview.note")])
    check_made_up(writer.files)
    low, high = (figures.stat(found.alignment["unique_pct"], h) for h in ("min", "max"))
    filled = pages.prose(shown, {
        "samples": sheet.height, "runs": sheet.height, "genome": text.word(PAGE, "preview.genome"),
        "length": f"{GENOME_LENGTH:,}", "release": text.word(PAGE, "preview.release"),
        "low": f"{low:.1f}", "high": f"{high:.1f}",
    })  # fmt: skip
    group = typer.main.get_command(app)
    usage, budget = runtime.measure(started), runtime.budget(MODULE)
    command = f"ctx process --dataset {dataset} --dataset-file {dataset}.toml --genome genome.fa.gz"
    built = ModulePage(
        title=text.fill(PAGE, "title", dataset=text.word(*SHOWN)),
        **(text.skeleton(PAGE) | filled),
        parameters=parameter_cards(group, group.commands[MODULE]),  # type: ignore[attr-defined]
        results=body, key_numbers=pages.key_numbers(found),
        reproduce=Reproduction(
            command=command, wall_seconds=usage.wall_seconds, peak_rss_self=usage.peak_rss_self,
            # page() starts no child process; the build's own ctx children are not this page's.
            peak_rss_children=0, profile=text.word(PAGE, "preview.profile"),
            budget=text.fill(PAGE, "preview.budget", minutes=budget.minutes, ram_gb=budget.ram_gb),
            versions=runtime.versions(), inputs=tuple(p.name for p in inputs),
        ),
    )  # fmt: skip
    return write(built, writer.out / f"{stem(MODULE, f'{dataset}_overview', 1, day=day)}.html")
