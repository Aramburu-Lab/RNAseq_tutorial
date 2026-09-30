"""The checkpoint on counts: per-sample metrics from the gene matrix, judged, and tested.

Three steps, in this order, because each reads the one before:

1. ``sample_metrics``: library size, genes detected, top-gene share, the four composition
   shares (disjoint gene classes from the GENCODE gene types), agreement with the peer group,
   distance from it on the PCA, and the sex check, imported from M5.
2. ``judge``: every registry metric the values hold, against its absolute band where it has
   one and against its stratum's median ± k scaled MAD always; a flag names its reason.
3. ``confounding``: whether a metric, or a flag, differs between the two groups of each
   comparison M5 makes, with M5's own sample selection. Nothing here excludes a sample.
"""

import gzip
import math
import re
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import fisher_exact, median_abs_deviation, rankdata

from ctx.dataset import Contrast
from ctx.de.deseq import pca, select_samples
from ctx.de.run import sex_check
from ctx.de.stats import adjust
from ctx.qc import QcError
from ctx.qc.registry import Metric, Registry
from ctx.report import text
from ctx.viz.annotate import TEST_NAME, compare_groups

PAGE = "qc"
MT_CHROM = "chrM"
MT_RRNA = ("Mt_rRNA",)
# GENCODE 50 types the 5S and 5.8S genes rRNA and their copies rRNA_pseudogene; it has no
# 45S gene, so 18S and 28S reads land on no gene at all (said on the page).
RRNA = ("rRNA", "rRNA_pseudogene")
GLOBIN = ("HBA1", "HBA2", "HBB")
ATTRIBUTE = re.compile(r'(gene_id|gene_name|gene_type) "([^"]*)"')
# ponytail: sex_check is the registry's one category metric; a per-metric flag value in the
# registry when a second one arrives.
FLAGGED_CATEGORY = "disagrees"
STRATUM = "stratum"
# One sample selection of M5: (sample, "reference" or "tested") pairs, sorted.
Selection = tuple[tuple[str, str], ...]
FLAG_COLUMNS = ["sample", "metric", "reason"]


def gene_types(gtf: Path) -> pl.DataFrame:
    """The gene records of a GENCODE GTF, gzipped or not: gene_id, gene_name, gene_type, chrom.

    Streams the file, so the 1.6 GB of text in GENCODE 50 never sits in memory.
    """
    opener = gzip.open if gtf.suffix == ".gz" else open
    rows = []
    with opener(gtf, "rt") as handle:
        for line in handle:
            if "\tgene\t" not in line:
                continue
            fields = line.split("\t", 8)
            if len(fields) == 9 and fields[2] == "gene":  # noqa: PLR2004 - GTF's 9 columns
                found = dict(ATTRIBUTE.findall(fields[8]))
                rows.append((found.get("gene_id"), found.get("gene_name"),
                             found.get("gene_type"), fields[0]))  # fmt: skip
    return pl.DataFrame(
        rows, schema=["gene_id", "gene_name", "gene_type", "chrom"], orient="row"
    ).unique("gene_id", keep="first", maintain_order=True)


def count_metrics(counts: pl.DataFrame, types: pl.DataFrame | None) -> pl.DataFrame:
    """Library size, genes detected, top-gene share and the four composition shares per sample.

    The shares are of disjoint gene classes: mitochondrial rRNA (Mt_rRNA), the other chrM
    genes, nuclear rRNA and its pseudogenes, and the three globin genes. Without gene types
    the shares are left out.

    Raises:
        QcError: If a gene of the counts is not in the annotation; the counts then came from
            another annotation, and every share would be wrong.
    """
    if types is None:
        matrix = counts.drop("gene_id").to_numpy().astype(float)
        size = matrix.sum(axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            top = matrix.max(axis=0) / size
        return pl.DataFrame(
            {
                "sample": counts.columns[1:],
                "library_size": size,
                "genes_detected": (matrix >= 1).sum(axis=0),
                "top_gene_share": top,
            }
        )
    typed = counts.join(types, on="gene_id", how="left")
    missing = typed.filter(pl.col("gene_type").is_null())["gene_id"]
    if missing.len():
        raise QcError(f"{missing.len():,} gene IDs of the counts are not in the annotation, "
                      f"for example {missing.head(10).to_list()}")  # fmt: skip
    names = counts.columns[1:]
    matrix = typed.select(names).to_numpy().astype(float)
    size = matrix.sum(axis=0)
    classes = {
        "mt_fraction": (pl.col("chrom") == MT_CHROM) & ~pl.col("gene_type").is_in(MT_RRNA),
        "mt_rrna_fraction": pl.col("gene_type").is_in(MT_RRNA),
        "rrna_fraction": pl.col("gene_type").is_in(RRNA) & (pl.col("chrom") != MT_CHROM),
        "globin_fraction": pl.col("gene_name").is_in(GLOBIN),
    }
    with np.errstate(divide="ignore", invalid="ignore"):
        shares = {
            name: matrix[typed.select(mask).to_series().to_numpy()].sum(axis=0) / size
            for name, mask in classes.items()
        }
        top = matrix.max(axis=0) / size
    return pl.DataFrame({
        "sample": names, "library_size": size, "genes_detected": (matrix >= 1).sum(axis=0),
        "top_gene_share": top, **shares,
    })  # fmt: skip


def replicate_spearman(
    counts: pl.DataFrame, samples: pl.DataFrame, peers: Sequence[str]
) -> pl.DataFrame:
    """Each sample's median Spearman correlation with the other samples of its peer group.

    Over the genes read in any sample of the group; a sample with no peer has none.
    """
    rows: list[tuple[str, float | None]] = []
    for _, group in samples.group_by(list(peers), maintain_order=True):
        names = group["sample"].to_list()
        if len(names) < 2:  # noqa: PLR2004 - a correlation needs a second sample
            rows += [(n, None) for n in names]
            continue
        matrix = counts.select(names).to_numpy()
        rho = np.corrcoef(rankdata(matrix[matrix.sum(axis=1) > 0], axis=0), rowvar=False)
        np.fill_diagonal(rho, np.nan)
        rows += [(n, float(np.nanmedian(r))) for n, r in zip(names, rho, strict=True)]
    return pl.DataFrame(rows, schema={"sample": pl.String, "replicate_spearman": pl.Float64},
                        orient="row")  # fmt: skip


def pca_distance(
    counts: pl.DataFrame, samples: pl.DataFrame, peers: Sequence[str], threads: int
) -> pl.DataFrame:
    """Each sample's distance from its peer group's median on PC1 and PC2 of one blind VST.

    One PCA of every sample, the one M5 draws, so the distance is read on the plot students
    see; a sample with no peer has none.
    """
    scores = pca(counts, samples["sample"].to_list(), threads=threads).scores
    scores = scores.join(samples.select("sample", *peers), on="sample", how="left")
    centre = [pl.col(c).median().over(peers) for c in ("pc1", "pc2")]
    distance = ((pl.col("pc1") - centre[0]) ** 2 + (pl.col("pc2") - centre[1]) ** 2).sqrt()
    alone = pl.len().over(peers) < 2  # noqa: PLR2004 - a centre needs a second sample
    return scores.select("sample", pca_distance=pl.when(alone).then(None).otherwise(distance))


def sample_metrics(
    counts: pl.DataFrame,
    sheet: pl.DataFrame,
    genes: pl.DataFrame,
    types: pl.DataFrame | None,
    peers: Sequence[str],
    threads: int,
) -> tuple[pl.DataFrame, pl.DataFrame | None]:
    """Every count-stage metric, one row per sample, and the sex check's own table.

    Args:
        counts: ``gene_id`` then one column per sample, as gather writes it.
        sheet: The harmonised sheet, one row per run.
        genes: Gather's ``gene_id`` and ``gene_name``, which the sex check reads.
        types: ``gene_types`` of the annotation the counts were made with; None leaves the
            composition shares out.
        peers: Sheet columns whose shared values make a peer group.
        threads: CPUs for the VST.

    Returns:
        The metrics, and the sex check's table (None where it could not run).
    """
    samples = sheet.unique("sample", keep="first", maintain_order=True)
    sexes, _ = sex_check(sheet, counts, genes)
    frame = count_metrics(counts, types)
    for part in (replicate_spearman(counts, samples, peers),
                 pca_distance(counts, samples, peers, threads)):  # fmt: skip
        frame = frame.join(part, on="sample", how="left")
    agrees = pl.Series("sex_agrees", [None] * frame.height, pl.String)
    if sexes is not None:
        label = (
            pl.when(pl.col("agrees").is_null())
            .then(pl.lit("not recorded"))
            .when(pl.col("agrees"))
            .then(pl.lit("agrees"))
            .otherwise(pl.lit(FLAGGED_CATEGORY))
        )
        frame = frame.join(sexes.select("sample", sex_agrees=label), on="sample", how="left")
    else:
        frame = frame.with_columns(agrees)
    return frame, sexes


def _show(value: float) -> str:
    """A value as a reason prints it: counts with thousands separators, others to 3 digits."""
    return f"{value:,.0f}" if abs(value) >= 1000 else f"{value:.3g}"  # noqa: PLR2004


def mad_band(values: np.ndarray, metric: Metric, k: float, scaled: bool) -> tuple[float, ...]:
    """Median, MAD and the band median ± k MAD, on the metric's scale, back in its units."""
    with np.errstate(divide="ignore"):
        x = np.log10(values) if metric.scale == "log10" else values
    median = float(np.median(x))
    mad = float(median_abs_deviation(x, scale="normal" if scaled else 1.0))
    centre_low_high = np.array([median, median - k * mad, median + k * mad])
    if metric.scale == "log10":
        centre_low_high = 10**centre_low_high
    centre, low, high = (float(v) for v in centre_low_high)
    return centre, mad, low, high


def _stratum_label(row: Mapping[str, object], strata: Sequence[str]) -> str:
    """A stratum as a reason names it: its values joined, e.g. "PB / NDN"."""
    return " / ".join(str(row[c]) for c in strata)


def _numeric(
    frame: pl.DataFrame, metric: Metric, registry: Registry, label: str
) -> tuple[dict[str, object], list[tuple[str, str]]]:
    """One stratum's MAD band for one numeric metric, and the samples outside it."""
    present = frame.filter(pl.col(metric.column).is_not_null() & pl.col(metric.column).is_not_nan())
    values = present[metric.column].to_numpy().astype(float)
    band = {"metric": metric.id, STRATUM: label, "n": len(values), "median": None, "mad": None,
            "low": None, "high": None, "note": ""}  # fmt: skip
    if len(values) < registry.min_samples:
        band["note"] = text.fill(PAGE, "notes.small_stratum", n=len(values),
                                 minimum=registry.min_samples)  # fmt: skip
        return band, []
    median, mad, low, high = mad_band(values, metric, registry.mad_multiplier, registry.mad_scaled)
    band.update(median=median, mad=mad, low=low, high=high)
    if mad == 0:
        band["note"] = text.word(PAGE, "notes.no_spread")
        return band, []
    flags = []
    for sample, v in present.select("sample", metric.column).rows():
        side = "low" if v < low else "high" if v > high else None
        if side and metric.direction in ("none", {"low": "higher", "high": "lower"}[side]):
            bound = low if side == "low" else high
            flags.append((sample, text.fill(PAGE, f"reasons.{side}", title=metric.title,
                          value=_show(v), bound=_show(bound), k=f"{registry.mad_multiplier:g}",
                          n=len(values), stratum=label)))  # fmt: skip
    return band, flags


def _absolute(frame: pl.DataFrame, metric: Metric, label: str) -> list[tuple[str, str]]:
    """The samples a metric flags without a MAD band.

    A value in the bad band, a disagreeing category, or a strand call that is undetermined or
    differs from the stratum's usual one.
    """
    rows = frame.select("sample", metric.column).drop_nulls().rows()
    if metric.kind == "category":
        return [(s, text.fill(PAGE, "reasons.category", title=metric.title, value=v))
                for s, v in rows if v == FLAGGED_CATEGORY]  # fmt: skip
    if metric.kind == "strandedness":
        calls = [(s, metric.strandedness(float(v))) for s, v in rows]
        usual = max(sorted({c for _, c in calls}), key=[c for _, c in calls].count, default="")
        return [
            (
                s,
                text.fill(
                    PAGE, "reasons.strand", title=metric.title, call=c, usual=usual, stratum=label
                ),
            )  # fmt: skip
            for s, c in calls
            if c == "undetermined" or c != usual
        ]
    return [
        (
            s,
            text.fill(
                PAGE,
                "reasons.band",
                title=metric.title,
                value=_show(float(v)) if metric.kind == "numeric" else v,
                source=metric.source,
            ),
        )  # fmt: skip
        for s, v in rows
        if metric.band(v) == "bad"
    ]


def judge(
    values: pl.DataFrame, samples: pl.DataFrame, registry: Registry, strata: Sequence[str]
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Every registry metric the values hold, judged per stratum.

    A numeric metric is judged against its absolute band where it has one, and always against
    its stratum's median ± k MAD, on its bad side only (both sides for direction "none"). A
    metric whose action is keep is judged and never flagged.

    Args:
        values: ``sample`` then one column per metric, named as the registry's ``column``.
        samples: One row per sample with the strata columns.
        registry: The loaded registry.
        strata: Sheet columns whose shared values make a stratum.

    Returns:
        The bands (one row per metric and stratum) and the flags (sample, metric, reason).
    """
    frame = values.join(samples.select("sample", *strata), on="sample", how="left")
    bands, flags = [], []
    for key, stratum in frame.group_by(list(strata), maintain_order=True):
        label = _stratum_label(dict(zip(strata, key, strict=True)), strata)
        for metric in registry.metrics.values():
            if metric.column not in values.columns:
                continue
            found = (
                _absolute(stratum, metric, label)
                if metric.kind != "numeric" or metric.banded
                else []
            )
            if metric.kind == "numeric":
                band, outside = _numeric(stratum, metric, registry, label)
                bands.append(band)
                found += outside
            if metric.action != "keep":
                flags += [(s, metric.id, r) for s, r in found]
    return (pl.DataFrame(bands, infer_schema_length=None),
            pl.DataFrame(flags, schema=FLAG_COLUMNS, orient="row"))  # fmt: skip


def status(runs: pl.DataFrame, flags: pl.DataFrame, registry: Registry) -> pl.DataFrame:
    """One row per run: qc_status pass, flag or exclude, the reasons, and how many.

    A run takes its sample's flags. Exclude only where a flagged metric's action says so;
    none does before the Step-1 record, so no run is excluded here.
    """
    excluding = [m.id for m in registry.metrics.values() if m.action == "exclude"]
    per = flags.group_by("sample", maintain_order=True).agg(
        n_flags=pl.col("metric").n_unique(),
        qc_reason=pl.col("reason").str.join("; "),
        excluded=pl.col("metric").is_in(excluding).any(),
    )
    joined = runs.select("run_accession", "sample").join(per, on="sample", how="left")
    return joined.select(
        "run_accession", "sample",
        qc_status=pl.when(pl.col("excluded")).then(pl.lit("exclude"))
        .when(pl.col("n_flags") > 0).then(pl.lit("flag")).otherwise(pl.lit("pass")),
        qc_reason=pl.col("qc_reason").fill_null(""),
        n_flags=pl.col("n_flags").fill_null(0).cast(pl.Int32),
    )  # fmt: skip


def selections(
    sheet: pl.DataFrame, contrasts: Sequence[Contrast], order: Mapping[str, tuple[str, ...]]
) -> tuple[dict[Selection, list[Contrast]], list[Contrast]]:
    """M5's samples for each unpaired contrast, contrasts with the same samples merged.

    A paired contrast is returned apart: its groups share individuals, which Mann-Whitney
    and Fisher's test assume they do not.
    """
    merged: dict[Selection, list[Contrast]] = {}
    paired = []
    for contrast in contrasts:
        if contrast.paired:
            paired.append(contrast)
            continue
        used, _ = select_samples(sheet, contrast, order)
        used = used.unique("sample", keep="first", maintain_order=True)
        merged.setdefault(tuple(sorted(used.select("sample", "group").rows())), []).append(contrast)
    return merged, paired


def fisher(tested: tuple[int, int], reference: tuple[int, int]) -> tuple[float, float]:
    """Fisher's exact test of (flagged, not flagged) between the groups, and the odds ratio.

    The odds ratio is tested over reference, with 0.5 added to every cell when one is 0
    (Haldane), so it is finite whenever the test is informative.
    """
    table = np.array([tested, reference], dtype=float)
    p = float(fisher_exact(table.astype(int)).pvalue)
    if (table == 0).any():
        table += 0.5
    return p, float(table[0, 0] * table[1, 1] / (table[0, 1] * table[1, 0]))


def _numeric_test(values: dict[str, list[float]]) -> dict[str, object]:
    """Mann-Whitney on one metric's values, through the M5 figures' own function."""
    reference, tested = values["reference"], values["tested"]
    if not reference or not tested:
        return {"note": text.word(PAGE, "notes.no_values")}
    if len(set(reference + tested)) == 1:
        return {"note": text.word(PAGE, "notes.constant")}
    long = pl.DataFrame({"f": "m", "group": ["reference"] * len(reference) + ["tested"] *
                         len(tested), "value": reference + tested})  # fmt: skip
    test = compare_groups(long, feature="f", group="group", value="value",
                          groups=("reference", "tested")).tests[0]  # fmt: skip
    return {"test": TEST_NAME, "statistic": test.u, "effect_name": text.word(PAGE, "effects.r"),
            "effect": test.r, "p": test.p}  # fmt: skip


def _flag_test(flagged: dict[str, list[bool]]) -> dict[str, object]:
    """Fisher's exact test on one flag between the groups."""
    cells = {g: (sum(v), len(v) - sum(v)) for g, v in flagged.items()}
    if not all(len(v) for v in flagged.values()):
        return {"note": text.word(PAGE, "notes.no_values")}
    if cells["tested"][0] + cells["reference"][0] == 0:
        return {"note": text.word(PAGE, "notes.no_flags")}
    if cells["tested"][1] + cells["reference"][1] == 0:
        return {"note": text.word(PAGE, "notes.all_flagged")}
    p, odds = fisher(cells["tested"], cells["reference"])
    return {"test": text.word(PAGE, "tests.fisher"), "statistic": None,
            "effect_name": text.word(PAGE, "effects.odds"), "effect": odds, "p": p}  # fmt: skip


def _tests(
    group: dict[str, str], values: pl.DataFrame, flags: pl.DataFrame, registry: Registry
) -> list[dict[str, object]]:
    """Every test of one sample selection.

    Numeric metrics by their values; the other metrics and the overall status by flag.
    """
    rows = []
    flagged = {(s, m) for s, m in flags.select("sample", "metric").rows()}
    any_flag = set(flags["sample"])
    chosen = values.filter(pl.col("sample").is_in(list(group)))
    for metric in registry.metrics.values():
        if metric.column not in values.columns:
            continue
        if metric.kind == "numeric":
            per: dict[str, list[float]] = {"reference": [], "tested": []}
            for s, v in chosen.select("sample", metric.column).rows():
                if v is not None and not math.isnan(v):
                    per[group[s]].append(float(v))
            result = _numeric_test(per)
        else:
            present = chosen.filter(pl.col(metric.column).is_not_null())["sample"]
            result = _flag_test({g: [(s, metric.id) in flagged for s in present if group[s] == g]
                                 for g in ("reference", "tested")})  # fmt: skip
        rows.append({"metric": metric.id, "read_length_tied": metric.read_length_tied, **result})
    overall = _flag_test({g: [s in any_flag for s in group if group[s] == g]
                          for g in ("reference", "tested")})  # fmt: skip
    return [*rows, {"metric": "qc_status", "read_length_tied": False, **overall}]


def confounding(
    values: pl.DataFrame,
    flags: pl.DataFrame,
    sheet: pl.DataFrame,
    contrasts: Sequence[Contrast],
    order: Mapping[str, tuple[str, ...]],
    registry: Registry,
) -> pl.DataFrame:
    """Does a QC metric, or a flag, differ between the two groups of a comparison M5 makes?

    Each distinct sample selection of the unpaired contrasts is tested once: numeric metrics
    by Mann-Whitney with the rank-biserial r, the other metrics' flags and the overall
    status by Fisher's exact test with the odds ratio; BH over every test in the table.

    Returns:
        One row per selection and metric, with the contrasts it serves and n per group; a
        paired contrast has one row saying why it is not tested.
    """
    merged, paired = selections(sheet, contrasts, order)
    rows = []
    for key, served in merged.items():
        group = dict(key)
        head = {"contrasts": "; ".join(c.name for c in served),
                "reference": served[0].reference, "tested": served[0].tested,
                "n_reference": sum(g == "reference" for g in group.values()),
                "n_tested": sum(g == "tested" for g in group.values())}  # fmt: skip
        rows += [{**head, **r} for r in _tests(group, values, flags, registry)]
    rows += [{"contrasts": c.name, "reference": c.reference, "tested": c.tested,
              "metric": "", "note": text.word(PAGE, "notes.paired")} for c in paired]  # fmt: skip
    table = pl.DataFrame(rows, infer_schema_length=None)
    for column, dtype in (("test", pl.String), ("statistic", pl.Float64),
                          ("effect_name", pl.String), ("effect", pl.Float64), ("p", pl.Float64),
                          ("note", pl.String)):  # fmt: skip
        if column not in table.columns:
            table = table.with_columns(pl.lit(None, dtype).alias(column))
    tested = table["p"].is_not_null().to_numpy()
    p_adj = np.full(table.height, np.nan)
    p_adj[tested] = adjust(table["p"].to_numpy()[tested], "padj")
    return table.with_columns(
        pl.Series("p_adj", p_adj).fill_nan(None), note=pl.col("note").fill_null("")
    )
