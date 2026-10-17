"""M5's lesson figures and the views a reader looks up on them (WO6, D-051, D-073).

The call counts of a lesson are drawn at every setting of the volcano's rule, the adjustment,
the cutoff and the fold change, each counted by ``ctx.de.stats.classify`` over the genes the
volcano tests; the page looks a setting up and never recomputes one. A pair of contrasts'
fold changes can be coloured by how their calls agree, and the named genes' boxes can carry
DESeq2's adjusted p over every tested gene instead of the Mann-Whitney p over the boxes drawn.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from typing import Any

import polars as pl
from scipy.stats import spearmanr

from ctx.de.deseq import RESULT_COLUMNS, Fit
from ctx.de.page import Words, figure
from ctx.de.stats import ADJUSTMENTS, ALPHA, LFC_MIN, classify
from ctx.report import text
from ctx.viz.annotate import Caption, wrap
from ctx.viz.bar import BarData
from ctx.viz.box import BoxData
from ctx.viz.scatter import CUTOFF_PRESETS, Line, ScatterData
from ctx.viz.style import DATA_COLOURS, ROLES

PAGE = "de"
# The volcano's bases, as its select lists them, raw p last (D-012).
BASES = (*ADJUSTMENTS, "pvalue")
# The volcano's cutoff presets and the fold changes the controls plan names (D-051 amendment,
# M5 batch lesson: "0, 0.5, 1 (default), 1.5, 2 log2").
ALPHAS = CUTOFF_PRESETS
FOLDS = (0.0, 0.5, 1.0, 1.5, 2.0)
# Agreement of two contrasts' calls, by role: called the same way takes the colour "called in
# either" takes in the default view; the yellow, 1.32:1 on white, is left out.
AGREEMENT = ("same", "x_only", "y_only", "opposite", "neither")
AGREEMENT_COLOURS = {"same": ROLES["comparison"], "x_only": DATA_COLOURS[0],
                     "y_only": DATA_COLOURS[1], "opposite": DATA_COLOURS[6],
                     "neither": None}  # fmt: skip


@dataclass(frozen=True)
class Setting:
    """One setting of the call rule: the p basis, the cutoff (strict), the fold change."""

    basis: str = "padj"
    alpha: float = ALPHA
    lfc_min: float = LFC_MIN


BUILD = Setting()


def tested(result: Fit) -> pl.DataFrame:
    """The genes the volcano draws: those with an adjusted p (``run.volcano``)."""
    return result.results.filter(pl.col("padj").is_not_null())


def counted(result: Fit, setting: Setting) -> tuple[int, int]:
    """Up and down calls of one contrast at one setting, by the one rule (D-012)."""
    calls = (
        tested(result)
        .select(
            classify(
                pl.col(setting.basis), pl.col("log2FoldChange"), setting.alpha, setting.lfc_min
            )
        )
        .to_series()
    )
    return int((calls == "up").sum()), int((calls == "down").sum())


def basis_words(basis: str) -> str:
    """A basis as a caption names it: ``BH-adjusted p``, ``raw p``."""
    return text.word(PAGE, f"figures.calls.bases.{basis}")


def calls_bars(title: str, said: Mapping[str, Words], fits: Mapping[str, Fit],
               sizes: str, setting: Setting = BUILD) -> BarData:  # fmt: skip
    """Up and down calls per contrast, for one lesson, at one setting of the rule.

    ``said`` holds the lesson's contrasts in words; ``sizes`` is the caption's group sizes.
    """
    rows = [(n, said[n].label, call, k) for n in said
            for call, k in zip(("up", "down"), counted(fits[n], setting), strict=True)]  # fmt: skip
    frame = pl.DataFrame(rows, schema=["contrast", "label", "call", "genes"], orient="row")
    rule = figure("calls", "rule", basis=basis_words(setting.basis), alpha=f"{setting.alpha:g}",
                  lfc_min=f"{setting.lfc_min:g}")  # fmt: skip
    columns = {
        "contrast": ("none", "Contrast."),
        "label": ("none", "The contrast in words, as the page names it."),
        "call": ("none", "up or down."),
        "genes": ("genes", f"Genes called: {rule}."),
    }
    return BarData(
        frame=frame, columns=columns, category="label", value="genes",
        category_label=figure("calls", "x"), value_label=figure("calls", "y"), series="call",
        colours={"up": ROLES["comparison"], "down": ROLES["reference"]},
        caption=figure("calls", "caption", lesson=title, rule=rule, sizes=sizes),
        horizontal=True, stacked=True, width="double",
    )  # fmt: skip


def settings() -> list[tuple[Setting, tuple[int, int, int]]]:
    """Every setting of the lesson grid with its choice positions, the build's first."""
    grid = [(Setting(b, a, f), (i, j, k)) for i, b in enumerate(BASES)
            for j, a in enumerate(ALPHAS) for k, f in enumerate(FOLDS)]  # fmt: skip
    return sorted(grid, key=lambda s: s[0] != BUILD)


def axes() -> list[tuple[str, list[str]]]:
    """The grid's selects: the p basis, the cutoff and the fold change, in words."""
    short = [text.word(PAGE, f"figures.calls.short.{b}") for b in BASES]
    return [(figure("calls", "axis_basis"), short),
            (figure("calls", "axis_alpha"), [f"{a:g}" for a in ALPHAS]),
            (figure("calls", "axis_fold"), [f"{f:g}" for f in FOLDS])]  # fmt: skip


def option_name(setting: Setting) -> str:
    """How the live line names a setting."""
    return figure("calls", "option", short=text.word(PAGE, f"figures.calls.short.{setting.basis}"),
                  alpha=f"{setting.alpha:g}", lfc_min=f"{setting.lfc_min:g}")  # fmt: skip


def grid_table(said: Mapping[str, Words], fits: Mapping[str, Fit]) -> pl.DataFrame:
    """Every count the lesson's grid draws: basis, cutoff, fold change, contrast, up, down."""
    rows = [(s.basis, s.alpha, s.lfc_min, n, *counted(fits[n], s))
            for s, _ in settings() for n in said]  # fmt: skip
    return pl.DataFrame(rows, schema=["basis", "alpha", "lfc_min", "contrast", "up", "down"],
                        orient="row")  # fmt: skip


GRID_COLUMNS = {
    "basis": ("none", "The p the rule reads: padj (BH), padj_by, padj_holm, padj_bonferroni or "
                      "pvalue (raw)."),
    "alpha": ("dimensionless", "The cutoff on that p, strict."),
    "lfc_min": ("log2 FC", "The smallest absolute raw log2 fold change called, inclusive."),
    "contrast": ("none", "Contrast."),
    "up": ("genes", "Genes called up at this setting, of those with an adjusted p."),
    "down": ("genes", "Genes called down at this setting, of those with an adjusted p."),
}  # fmt: skip


def lfc_scatter(x: Fit, y: Fit, genes: pl.DataFrame, said: tuple[Words, Words]
                ) -> tuple[ScatterData, float]:  # fmt: skip
    """One contrast's raw LFC against another's, on genes tested in both."""
    called, neither = figure("lfc", "called"), figure("lfc", "neither")
    frame = (
        x.results.select(
            "gene_id", pl.col("log2FoldChange").alias("lfc_x"), pl.col("call").alias("call_x")
        )
        .join(
            y.results.select(
                "gene_id", pl.col("log2FoldChange").alias("lfc_y"), pl.col("call").alias("call_y")
            ),
            on="gene_id",
        )
        .drop_nulls()
        .join(tested(x).select("gene_id").join(tested(y).select("gene_id"), on="gene_id"),
              on="gene_id")
        .join(genes, on="gene_id", how="left")
        .with_columns(
            called=pl.when((pl.col("call_x") != "ns") | (pl.col("call_y") != "ns"))
            .then(pl.lit(called))
            .otherwise(pl.lit(neither))
        )
        .with_columns(keep=pl.col("called") == called)
    )  # fmt: skip
    rho = float(spearmanr(frame["lfc_x"], frame["lfc_y"]).statistic)
    shared = set(x.samples["donor"]) & set(y.samples["donor"])
    note = figure("lfc", "shared", n=len(shared)) if shared else figure("lfc", "apart")
    ends = (float(frame.select(pl.min_horizontal("lfc_x", "lfc_y").min()).item()),
            float(frame.select(pl.max_horizontal("lfc_x", "lfc_y").max()).item()))  # fmt: skip
    columns = {
        "gene_id": RESULT_COLUMNS["gene_id"],
        "gene_name": ("none", "Gene symbol."),
        "lfc_x": ("log2 FC", f"Raw log2 fold change in {x.contrast.name}."),
        "lfc_y": ("log2 FC", f"Raw log2 fold change in {y.contrast.name}."),
        "call_x": ("none", f"Call in {x.contrast.name}."),
        "call_y": ("none", f"Call in {y.contrast.name}."),
        "called": ("none", "Called in either contrast, or neither."),
        "keep": ("none", "Called in either; always drawn on the page."),
    }
    data = ScatterData(
        frame=frame, columns=columns, x="lfc_x", y="lfc_y",
        x_label=figure("lfc", "axis", label=said[0].label),
        y_label=figure("lfc", "axis", label=said[1].label),
        colour="called", colours={called: ROLES["comparison"], neither: None},
        caption=figure("lfc", "caption", genes=frame.height, rho=rho, note=note),
        hover=("gene_name",), lines=(Line(figure("lfc", "line"), ends, ends),),
        width="double", keep="keep", thinned=text.word(PAGE, "figures.lfc.thinned"),
    )  # fmt: skip
    return data, rho


def agreement_level(call_x: pl.Expr, call_y: pl.Expr) -> pl.Expr:
    """How two calls agree: ``AGREEMENT``'s keys, the one definition the view and its table use."""
    x, y = call_x != "ns", call_y != "ns"
    return (pl.when(x & y & (call_x == call_y)).then(pl.lit("same"))
            .when(x & y).then(pl.lit("opposite"))
            .when(x).then(pl.lit("x_only"))
            .when(y).then(pl.lit("y_only"))
            .otherwise(pl.lit("neither")))  # fmt: skip


def agreement(data: ScatterData, said: tuple[Words, Words]
              ) -> tuple[ScatterData, dict[str, Any]]:  # fmt: skip
    """The fold-change scatter coloured by how the two calls agree, and its numbers.

    The numbers: each level's genes, the genes called in each contrast, and Spearman's rho
    over the genes called in either, which the caption's rho over every gene does not show.
    Each legend entry carries its count, so the view's files state them.
    """
    level = agreement_level(pl.col("call_x"), pl.col("call_y"))
    found = dict(data.frame.select(level.alias("a"))["a"].value_counts().rows())
    counts = {k: int(found.get(k, 0)) for k in AGREEMENT}
    words = {k: figure("lfc", "agree_entry", n=counts[k],
                       level=figure("lfc", f"agree_{k}", x=said[0].label, y=said[1].label))
             for k in AGREEMENT}  # fmt: skip
    frame = data.frame.with_columns(agreement=level.replace_strict(words))
    either = frame.filter(pl.col("keep"))
    rho = None
    if either.height > 2:  # noqa: PLR2004 - a rank correlation needs three
        rho = float(spearmanr(either["lfc_x"], either["lfc_y"]).statistic)
    numbers = {**counts, "called_x": counts["same"] + counts["opposite"] + counts["x_only"],
               "called_y": counts["same"] + counts["opposite"] + counts["y_only"],
               "either": either.height, "rho_called": rho}  # fmt: skip
    rho_words = (figure("lfc", "rho_called", either=either.height, rho=rho) if rho is not None
                 else figure("lfc", "rho_none"))  # fmt: skip
    caption = f"{data.caption} {rho_words}"
    columns = {**data.columns, "agreement": ("none", "How the two calls agree.")}
    shown = replace(data, frame=frame, columns=columns, colour="agreement",
                    colours={words[k]: AGREEMENT_COLOURS[k] for k in AGREEMENT},
                    caption=caption)  # fmt: skip
    return shown, numbers


AGREEMENT_TABLE = {
    "agreement": ("none", "same, x_only, y_only, opposite or neither: how the two calls agree."),
    "genes": ("genes", "Genes tested in both contrasts at that agreement."),
    "rho_called": ("dimensionless", "Spearman rho of the raw log2 FC over genes called in either "
                                    "(every row the same)."),
}  # fmt: skip


@dataclass(frozen=True)
class DeseqCaption(Caption):
    """The named genes' boxes labelled with DESeq2's padj, BH over every tested gene.

    Each gene's label is its padj and raw log2 FC from the contrast's results; the caption
    names that family, so the boxes read as the test the volcano draws.
    """

    padj: Mapping[str, float | None] = field(default_factory=dict)
    lfc: Mapping[str, float | None] = field(default_factory=dict)
    n_tested: int = 0

    def label(self, feature: str, key: str = "label") -> str:
        """DESeq2's padj and log2 FC above one gene; ``label_line`` on one line."""
        p, lfc = self.padj.get(feature), self.lfc.get(feature)
        if p is None or lfc is None:
            return figure("box", f"deseq_{key}_none")
        return figure("box", f"deseq_{key}", p=f"{p:.2g}", lfc=lfc)

    def text(self, width: int = 0) -> str:
        """The caption: the family BH ran over, n per group, the sign of the fold change."""
        reference, comparison = self.groups
        ns = sorted({(t.n_reference, t.n_comparison) for t in self.tests})
        sizes = "; ".join(
            text.fill("annotate", "sizes", n_reference=a, reference=reference, n_comparison=b,
                      comparison=comparison)
            for a, b in ns
        )  # fmt: skip
        body = figure("box", "deseq_caption", tested=self.n_tested, sizes=sizes,
                      comparison=comparison)  # fmt: skip
        return wrap(body, width)


def deseq_box(box: BoxData, result: Fit, symbols: Sequence[tuple[str, str]]) -> BoxData:
    """The same boxes, each gene labelled with DESeq2's padj and raw log2 FC (WO6, F12).

    ``symbols`` pairs each drawn symbol with its gene identifier, as the boxes chose them.
    """
    found = {g: (p, f) for g, p, f in
             result.results.select("gene_id", "padj", "log2FoldChange").rows()}  # fmt: skip
    padj = {s: found.get(g, (None, None))[0] for s, g in symbols}
    lfc = {s: found.get(g, (None, None))[1] for s, g in symbols}
    caption = DeseqCaption(box.caption.groups, box.caption.tests, padj, lfc,
                           int(result.results["padj"].is_not_null().sum()))  # fmt: skip
    return replace(box, caption=caption)
