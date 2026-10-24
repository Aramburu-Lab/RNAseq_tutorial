"""M3's figure contracts, built from ``ctx.process.tables``' frames; every word from content.

Four kinds: where each sample's reads align (F1), and the featured sample's under each setting
(F1b, a figure of its own so its few bars keep their labels however many samples F1 draws);
what featureCounts assigns at each -s (F2); one counter against another gene by gene (F3) and
how far apart they are (F4). The sample bars filter by the sheet's groups and the strand bars
also by call (D-051); the setting bars, at most five of one sample, and the scatters carry no
control, for the reasons ``TODO.md`` P2S.6 records. Captions come from
``content/process.toml``: these figures estimate agreement or describe, and test nothing
(D-040, R6.1), so ``ctx.viz.annotate`` writes none. Hover and filters name every column in
words from the same file.
"""

from collections.abc import Mapping, Sequence
from typing import Any

import polars as pl

from ctx.concordance import LEVEL
from ctx.process import CATEGORIES, Aligner
from ctx.process.tables import GENE_TYPES, OTHER
from ctx.report import text
from ctx.viz.bar import BarData
from ctx.viz.scatter import Line, ScatterData
from ctx.viz.style import BAND_COLOURS, DATA_COLOURS, category_colours

PAGE = "process"
# Each outcome's colour by its role, never by its place in a list: a read placed once is good,
# one left unaligned bad, as the QC bands colour them; the two kinds of multimapper are blues;
# STAR's reasons for leaving a read unaligned take the remaining colours. No outcome takes the
# yellow a projector washes out.
OUTCOME_COLOURS = {
    "unique": BAND_COLOURS["good"], "multi": DATA_COLOURS[1], "too_many": DATA_COLOURS[4],
    "mismatches": DATA_COLOURS[6], "too_short": BAND_COLOURS["usable"],
    "other": BAND_COLOURS["bad"], "unmapped": BAND_COLOURS["bad"],
}  # fmt: skip
TYPE_LEVELS = (*GENE_TYPES, OTHER)
# The empirical percentiles of d beside the limits of agreement (D-040).
PERCENTILES = (2.5, 97.5)
# Each caption placeholder and the concordance column it prints.
SCATTER_FIELDS = {
    "rho": "spearman", "rho_low": "spearman_low", "rho_high": "spearman_high", "r": "pearson",
    "r_low": "pearson_low", "r_high": "pearson_high", "bias": "bias", "bias_low": "bias_low",
    "bias_high": "bias_high", "loa_low": "loa_low", "loa_high": "loa_high", "d_low": "d_p2_5",
    "d_high": "d_p97_5", "fold": "fold",
}  # fmt: skip


def stat(series: pl.Series, how: str) -> float:
    """A column's min, max or median as a float; Polars types the result loosely."""
    value: Any = getattr(series, how)()
    return float(value)


def reads(value: float) -> str:
    """A number of reads as text: whole with thousands separators, a median's half kept."""
    return f"{value:,.0f}" if value == int(value) else f"{value:,.1f}"


def number(value: float | None, spec: str) -> str:
    """A number as a caption prints it, or the content's word for a value that is undefined."""
    return text.word(PAGE, "undefined") if value is None else format(value, spec)


def _columns(key: str, frame: pl.DataFrame, extra: Mapping[str, str]) -> dict[str, tuple[str, str]]:
    """(unit, meaning) of every column, from ``figures.<key>``; ``extra`` describes the rest."""
    units = text.raw(PAGE, f"figures.{key}.units")
    return {c: (units.get(c, "none"), extra[c] if c in extra
                else text.word(PAGE, f"figures.{key}.columns.{c}"))
            for c in frame.columns}  # fmt: skip


def group_columns(groups: Sequence[str]) -> dict[str, str]:
    """The meaning of each sheet column a figure filters by."""
    return {g: text.fill(PAGE, "figures.group_column", column=g) for g in groups}


def labels(key: str, columns: Sequence[str]) -> dict[str, str]:
    """Words for columns on hover and on filters: ``figures.<key>.labels``, then the filters'."""
    words = {**text.raw(PAGE, "figures.filter_labels"), **text.raw(PAGE, f"figures.{key}.labels")}
    return {c: words.get(c, c) for c in columns}


def outcome_bars(frame: pl.DataFrame, aligner: Aligner, key: str, groups: Sequence[str],
                 caption: str) -> BarData:  # fmt: skip
    """F1 and F1b: one stacked bar per sample (``mapping``) or per setting (``settings``).

    The filter's live line prints no median: every bar is stacked to 100 % (P2S.6 (d)).

    Args:
        frame: label, outcome, reads, input_reads, percent, and the ``groups`` columns.
        aligner: Its outcomes are the series, in the order it reports them.
        key: The figure's words in content, ``mapping`` or ``settings``.
        groups: Sheet columns to filter by; none for the settings.
        caption: The filled caption.
    """
    levels = CATEGORIES[aligner]
    names = {o: text.word(PAGE, f"outcomes.{o}") for o in levels}
    return BarData(
        frame=frame, columns=_columns(key, frame, group_columns(groups)), category="label",
        value="percent", category_label=text.word(PAGE, f"figures.{key}.category"),
        value_label=text.word(PAGE, "figures.mapping.value"), series="outcome",
        colours={o: OUTCOME_COLOURS[o] for o in levels}, caption=caption,
        horizontal=True, stacked=True, hover=("reads", "input_reads"), width="double",
        names=names, filters=tuple(groups),
        labels=labels("mapping", ["outcome", "reads", "input_reads", *groups]),
        summary_median=False,
    )  # fmt: skip


def strand_bars(frame: pl.DataFrame, groups: Sequence[str], caption: str) -> BarData:
    """F2: reads featureCounts assigns at -s 0, 1 and 2, as % of input, per sample.

    The filter's live line prints no median: one would pool the three settings (P2S.6 (d)).

    Args:
        frame: sample, setting, percent, reads, no_feature, ambiguous, forward_share, call,
            and the ``groups`` columns.
        groups: Sheet columns to filter by, beside the call.
        caption: The filled caption.
    """
    settings = ("s0", "s1", "s2")
    hover = ("reads", "no_feature", "ambiguous", "forward_share", "call")
    return BarData(
        frame=frame, columns=_columns("strandedness", frame, group_columns(groups)),
        category="sample", value="percent",
        category_label=text.word(PAGE, "figures.strandedness.category"),
        value_label=text.word(PAGE, "figures.strandedness.value"), series="setting",
        colours=category_colours(list(settings)), caption=caption, hover=hover, width="double",
        names=text.raw(PAGE, "figures.strandedness.names"), filters=("call", *groups),
        labels=labels("strandedness", ["setting", *hover, *groups]), summary_median=False,
    )  # fmt: skip


def comparison_words(name: str) -> dict[str, str]:
    """A comparison's words: its name, and what its a and b are."""
    return {part: text.word(PAGE, f"comparisons.{name}.{part}") for part in ("name", "a", "b")}


def scatter_values(every: Mapping[str, Any], per_sample: pl.DataFrame,
                   samples: int) -> dict[str, Any]:  # fmt: skip
    """The placeholders of both agreement figures, each number formatted as captions print it.

    Args:
        every: The comparison's ``all`` row of the concordance table.
        per_sample: Its per-sample rows.
        samples: How many samples were summed.
    """
    rho = per_sample["spearman"].drop_nulls()
    value = {
        "n": f"{every['n_genes']:,}", "n_same": f"{every['n_identical']:,}",
        "n_diff": f"{every['n_genes'] - every['n_identical']:,}", "samples": samples,
        "ci": f"{LEVEL * 100:g}", "p_low": f"{PERCENTILES[0]:g}", "p_high": f"{PERCENTILES[1]:g}",
        "mid": f"{PERCENTILES[1] - PERCENTILES[0]:g}",
        "log_unit": text.word(PAGE, "log_unit"),
        **{f"rho_{how[:3]}": number(stat(rho, how) if rho.len() else None, ".3f")
           for how in ("median", "min", "max")},
    }  # fmt: skip
    return value | {k: number(every[f], ".3f") for k, f in SCATTER_FIELDS.items()}


def _typed(sums: pl.DataFrame) -> pl.DataFrame:
    """A comparison's genes with a read in either count, their logs, and whether they differ.

    A count whose every value is whole is typed as one, so hover prints 1,640, not 1.64e+3.
    """
    kept = sums.filter(pl.col("a") + pl.col("b") > 0)
    whole = [c for c in ("a", "b") if (kept[c] == kept[c].round()).all()]
    return kept.with_columns(pl.col(whole).cast(pl.Int64)).with_columns(
        log2_a=(pl.col("a") + 1).log(2), log2_b=(pl.col("b") + 1).log(2),
        differs=pl.col("a") != pl.col("b"))  # fmt: skip


def _hover_words(words: Mapping[str, str]) -> dict[str, str]:
    """The scatters' hover columns in words: the gene, and each count by its comparison."""
    return {**labels("counters", ["gene_name", "gene_id"]), "a": words["a"], "b": words["b"]}


def limits_of(row: Mapping[str, Any], empirical: bool = False) -> tuple[float | None, float | None]:
    """A concordance row's limits: bias plus and minus 1.96 SD, or d's 2.5th and 97.5th."""
    low, high = ("d_p2_5", "d_p97_5") if empirical else ("loa_low", "loa_high")
    return row[low], row[high]


def beyond(limits: tuple[float | None, float | None]) -> pl.Expr:
    """Whether a gene's d lies beyond the limits, as the leads count it; null with no limits."""
    low, high = limits
    d = pl.col("log2_b") - pl.col("log2_a")
    if low is None or high is None:
        return pl.lit(None, pl.Boolean)
    return ~d.is_between(low, high)


def _density(frame: pl.DataFrame) -> dict[str, Any]:
    """The page's densities (D-118): genes beyond the limits as dots, the rest as regions."""
    if frame["outside"].null_count():
        return {}
    # A template: ScatterData fills it from what the page draws.
    return {"outside": "outside", "dense_words": str(text.raw(PAGE, "figures.dense"))}


def counter_scatter(sums: pl.DataFrame, words: Mapping[str, str], caption: str,
                    limits: tuple[float | None, float | None] = (None, None)
                    ) -> ScatterData:  # fmt: skip
    """F3: one gene per point, its two counts summed over samples, on log2(count + 1) axes.

    Points are coloured by gene type. With ``limits``, the view's limits of agreement, the page
    draws each gene type as its densities and every gene beyond the limits as a dot (D-118);
    without them every gene whose two counts differ is kept on the page whatever the thinning,
    so no disagreement is hidden.
    """
    frame = _typed(sums).select("gene_id", "gene_name", "type_group", "a", "b", "log2_a",
                                "log2_b", "differs", outside=beyond(limits))  # fmt: skip
    span = (0.0, max(stat(frame["log2_a"], "max"), stat(frame["log2_b"], "max")))
    extra = {"a": text.fill(PAGE, "figures.counters.columns.a", a=words["a"]),
             "b": text.fill(PAGE, "figures.counters.columns.b", b=words["b"])}  # fmt: skip
    return ScatterData(
        frame=frame, columns=_columns("counters", frame, extra), x="log2_a", y="log2_b",
        x_label=words["a"], y_label=words["b"], colour="type_group",
        colours=category_colours(list(TYPE_LEVELS)),
        names={t: text.word(PAGE, f"gene_types.{t}") for t in TYPE_LEVELS},
        caption=caption, hover=("gene_name", "gene_id", "a", "b"), keep="differs",
        lines=(Line(text.word(PAGE, "figures.counters.identity"), span, span),),
        labels=_hover_words(words), **_density(frame),
    )  # fmt: skip


def agreement_plot(sums: pl.DataFrame, every: Mapping[str, Any], words: Mapping[str, str],
                   caption: str, empirical: bool = False) -> ScatterData:  # fmt: skip
    """F4: Bland and Altman's plot, d = log2(b + 1) - log2(a + 1) against the pair's mean.

    The bias is solid, the limits dashed: bias plus and minus 1.96 SD, or with ``empirical``
    the 2.5th and 97.5th percentiles of d (WO8, D-051). A line whose value is undefined is not
    drawn. The page draws each gene type as its densities and every gene beyond the limits
    drawn as a dot (D-118).
    """
    limits = ("d_p2_5", "d_p97_5") if empirical else ("loa_low", "loa_high")
    frame = _typed(sums).select(
        "gene_id", "gene_name", "type_group", "a", "b",
        mean_log2=(pl.col("log2_a") + pl.col("log2_b")) / 2,
        diff_log2=pl.col("log2_b") - pl.col("log2_a"), differs="differs",
        outside=beyond(limits_of(every, empirical)))  # fmt: skip
    span = (stat(frame["mean_log2"], "min"), stat(frame["mean_log2"], "max"))
    lines = []
    if every["bias"] is not None:
        lines.append(Line(text.word(PAGE, "figures.agreement.bias"), span, (every["bias"],) * 2))
    label = text.word(PAGE, "figures.agreement.percentiles" if empirical
                      else "figures.agreement.limits")  # fmt: skip
    for limit in limits:
        if every[limit] is not None:
            lines.append(Line(label, span, (every[limit],) * 2, dash=True))
    extra = {"a": text.fill(PAGE, "figures.counters.columns.a", a=words["a"]),
             "b": text.fill(PAGE, "figures.counters.columns.b", b=words["b"])}  # fmt: skip
    return ScatterData(
        frame=frame, columns=_columns("agreement", frame, extra), x="mean_log2", y="diff_log2",
        x_label=text.word(PAGE, "figures.agreement.x"),
        y_label=text.word(PAGE, "figures.agreement.y"), colour="type_group",
        colours=category_colours(list(TYPE_LEVELS)),
        names={t: text.word(PAGE, f"gene_types.{t}") for t in TYPE_LEVELS},
        caption=caption, hover=("gene_name", "gene_id", "a", "b"), keep="differs",
        lines=tuple(lines), labels=_hover_words(words), **_density(frame),
    )  # fmt: skip
