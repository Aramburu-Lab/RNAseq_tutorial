"""The drawn own-sample panels of a card (E.1c2): the featured contrast's box and F1.

Both are ``BoxData``, drawn by ``ctx.viz.box``, with M9's own caption: its family is the card
list (the box) or the contrasts drawn (F1), not the features of one figure, and a paired
contrast's test is the signed-rank, which ``annotate.compare_groups`` does not run. The values
and tests come from ``samples``.
"""

import textwrap
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from typing import Any

import polars as pl

from ctx.dataset import Dataset
from ctx.de.page import Words, every
from ctx.de.run import Writer
from ctx.de.stats import ALPHA
from ctx.genecard import samples
from ctx.genecard.cards import Card
from ctx.report import text
from ctx.viz.annotate import Caption, FeatureTest, wrap
from ctx.viz.box import STAT_COLUMNS, BoxData, least_width
from ctx.viz.box import interactive as box_interactive
from ctx.viz.style import HEADROOM

PAGE = "genecard"
# Characters per line of a contrast's name under its box, at PER_FIGURE boxes in 170 mm.
LABEL_CHARS = 14
# Contrasts per F1 figure: six such labels fit 170 mm with no tick label hidden; eleven hid
# every other one (Step-1, 2026-10-17).
PER_FIGURE = 6
# F1 offers a view per batch where the sample sheet has two to this many (Step-1 has two).
MAX_BATCHES = 3
# The margin ``box.interactive`` leaves under the lowest value, as a fraction of the values' span.
BELOW = 0.05


@dataclass(frozen=True)
class OwnCaption(Caption):
    """A box caption written by M9: its family and test words, and n over each box for F1."""

    words: str = ""
    with_n: bool = False

    def label(self, feature: str, key: str = "label") -> str:
        """The p and r over a feature, as every box prints them; F1 puts each group's n first."""
        found = super().label(feature, key)
        if not self.with_n:
            return found
        t = next(t for t in self.tests if t.feature == feature)
        n = text.fill(PAGE, "figures.contrasts.n", reference=t.n_reference, tested=t.n_comparison)
        return f"{n}\n{found}" if key == "label" else f"{n}, {found}"

    def text(self, width: int = 0) -> str:
        """The caption M9 filled, wrapped to ``width`` characters; 0 leaves one line."""
        return wrap(self.words, width)


@dataclass(frozen=True)
class OwnBox(BoxData):
    """A box whose data file names its own test's statistic, not Mann-Whitney's alone."""

    stats: Mapping[str, tuple[str, str]] = field(default_factory=dict)

    def table_columns(self) -> dict[str, tuple[str, str]]:
        """(unit, meaning) of every column of ``table()``, the statistics in M9's words."""
        return {**self.columns, **self.stats}

    def interactive(self) -> Any:
        """The plotly box with each wrapped feature name on its lines, as the files draw it.

        plotly prints a category's line break as nothing, and ``box.least_width`` measured the
        name as one line: F1's six contrasts asked for a phone's layout at 1,009 px. The ticks
        name each category with ``<br>`` between its lines, and the upright layout's least width
        is measured on each name's widest line. The phone's layout, turned labels, is kept only
        where it is the narrower of the two; ``box.narrow_width`` counts its lines already.
        """
        fig = box_interactive(self)
        if isinstance(self.caption, OwnCaption) and self.caption.with_n:
            # Three lines over each box, n then p and r, where the box leaves room for two:
            # F1's n ran into the legend at 1,009 px.
            low, high = fig.layout.yaxis.range
            span = (high - low) / (1 + BELOW + HEADROOM)
            fig.update_yaxes(range=[low, high + span * HEADROOM / 2])
        names = self.frame["feature"].unique(maintain_order=True).to_list()
        if not any("\n" in n for n in names):
            return fig
        fig.update_xaxes(tickmode="array", tickvals=names,
                         ticktext=[n.replace("\n", "<br>") for n in names])  # fmt: skip
        meta = fig.layout.meta
        relayout = meta["narrow"]["relayout"]
        relayout["xaxis.ticktext"] = [t.replace("\n", "<br>") for t in relayout["xaxis.ticktext"]]
        widest = {n: max(n.split("\n"), key=len) for n in names}
        one_line = replace(
            self, frame=self.frame.with_columns(pl.col("feature").replace_strict(widest)),
            caption=replace(self.caption, tests=tuple(replace(t, feature=widest[t.feature])
                                                      for t in self.caption.tests)),
        )  # fmt: skip
        upright = least_width(one_line, tuple(fig.layout.yaxis.range))
        if sum(meta["least_width"][0]) < sum(upright[0]):
            fig.update_layout(meta={**meta, "narrow": {**meta["narrow"], "below": sum(upright[0])}})
        else:
            # Turned, F1's labels need more room than upright (1,184 px against 621 for
            # Step-1's six): the upright layout is kept at every width and scrolls below it.
            fig.layout.meta = {"least_width": upright}
        return fig


def _columns(key: str) -> dict[str, tuple[str, str]]:
    """(unit, meaning) of one figure's columns from the content file."""
    return {c: (u, m) for c, (u, m) in text.raw(PAGE, f"figures.{key}.columns").items()}


def _p(value: float) -> str:
    """A p-value as the box labels print it: two significant digits."""
    return f"{value:.2g}"


def _test_words(paired: bool) -> str:
    """The test's name in the caption."""
    return text.word(PAGE, f"figures.own.{'paired' if paired else 'unpaired'}")


def _levels(said: Words) -> tuple[str, str]:
    """The two groups' names: the contrast's levels, or its roles where the levels read alike."""
    if said.reference == said.tested:
        roles = text.raw(PAGE, "figures.contrasts.groups")
        return roles[0], roles[1]
    return said.reference, said.tested


def box(card: Card, values: samples.Values, contrast: str, said: Words,
        tests: Mapping[str, samples.Test], padj: float | None
        ) -> tuple[OwnBox, dict[str, Any]] | None:  # fmt: skip
    """The gene in the featured contrast's samples, its test from the card list's family.

    Returns:
        The figure and its numbers; None where the family holds no test of the gene.
    """
    t = tests.get(card.gene_id)
    if t is None:
        return None
    reference, tested = _levels(said)
    named = dict(zip(samples.GROUPS, (reference, tested), strict=True))
    rows = values.frame.filter((pl.col("gene_id") == card.gene_id)
                               & (pl.col("contrast") == contrast))  # fmt: skip
    frame = rows.select(pl.lit(card.symbol).alias("feature"),
                        group=pl.col("group").replace_strict(named),
                        sample="sample", donor="donor", value="value",
                        deseq2_padj=pl.lit(padj, dtype=pl.Float64))  # fmt: skip
    paired = values.paired.get(contrast, False)
    numbers = {"symbol": card.symbol, "contrast": said.label, "reference": reference,
               "tested": tested, "n_reference": t.n_reference, "n_tested": t.n_tested,
               "family": len(tests), "test": _test_words(paired), "p": _p(t.p_adj), "r": t.r,
               "padj": text.word(PAGE, "none") if padj is None else _p(padj)}  # fmt: skip
    caption = OwnCaption(
        (reference, tested),
        (FeatureTest(card.symbol, t.n_reference, t.n_tested, t.statistic, t.p, t.p_adj, t.r),),
        words=text.fill(PAGE, "figures.box.caption", **numbers),
    )
    data = OwnBox(frame=frame, columns=_columns("box"),
                  feature_label=text.word(PAGE, "figures.box.x"),
                  value_label=text.word(PAGE, "figures.box.y"), unit=_columns("box")["value"][0],
                  groups=(reference, tested), hover=("sample", "donor", "deseq2_padj"),
                  caption=caption, stats=_stats(paired, not paired))  # fmt: skip
    return data, numbers


def _stats(paired: bool, unpaired: bool) -> dict[str, tuple[str, str]]:
    """The statistic columns' meanings: the statistic named for the tests the figure holds."""
    words = text.raw(PAGE, "figures.own.stats")
    kind = "u_both" if paired and unpaired else "u_paired" if paired else "u"
    return {k: (words[kind if k == "u" else k][0], words[kind if k == "u" else k][1])
            for k in STAT_COLUMNS}  # fmt: skip


def names(spec: Dataset) -> dict[str, str]:
    """Each contrast's words by name, in the dataset file's order."""
    labels = every(spec.contrasts, spec.columns)
    return {c.name: labels[c.name].label for c in spec.contrasts}


def write_tables(writer: Writer, values: samples.Values, spec: Dataset,
                 chosen: Sequence[Card]) -> None:  # fmt: skip
    """Every value the own-sample panels of the kept cards draw, and every F1 view's tests."""
    cards = {c.gene_id: c for c in chosen}
    labels = names(spec)
    frame = values.frame.filter(pl.col("gene_id").is_in(list(cards))).select(
        symbol=pl.col("gene_id").replace_strict({g: c.symbol for g, c in cards.items()}),
        gene_id="gene_id", contrast=pl.col("contrast").replace_strict(labels),
        sample="sample", donor="donor", group="group", batch="batch", value="value",
    )  # fmt: skip
    writer.table(frame, "own_values", _columns("own_values"))
    if len(spec.contrasts) > 1:
        levels = samples.batches(values, MAX_BATCHES)
        tests = [t for c in chosen if (t := table(c, values, labels, levels)).height]
        if tests:
            writer.table(pl.concat(tests), "contrast_tests", _columns("contrast_tests"))


def label(words: str) -> str:
    """A contrast's name wrapped under its box."""
    return textwrap.fill(words, LABEL_CHARS, break_on_hyphens=False)


@dataclass(frozen=True)
class Family:
    """One F1 view's tests: the gene in each contrast within a batch or all, BH over them.

    Attributes:
        name: The view's samples in words: all, or one batch.
        batch: The batch level, None for every sample.
        tests: Each contrast tested, in the dataset's order, with p adjusted over all of them.
        left: The contrasts with too few values in a group, or nothing to rank.
    """

    name: str
    batch: str | None
    tests: Mapping[str, samples.Test]
    left: Sequence[str]


def pages(names: Sequence[str], most: int = PER_FIGURE) -> list[list[str]]:
    """Contrasts split into the fewest figures of at most ``most``, as even as they go."""
    count = -(-len(names) // most)
    size = -(-len(names) // count) if count else 0
    return [list(names[i : i + size]) for i in range(0, len(names), size)] if size else []


def contrasts(card: Card, values: samples.Values, labels: Mapping[str, str],
              family: Family, drawn: Sequence[str]) -> tuple[OwnBox, dict[str, Any]]:  # fmt: skip
    """F1: the gene in the ``drawn`` contrasts of one view, p adjusted over the view's family.

    Each contrast's values are its own fit's (``samples``). The family holds every contrast
    the view tests, drawn here or on another of its figures, so a p does not move with paging.
    """
    tests = {c: family.tests[c] for c in drawn}
    roles = text.raw(PAGE, "figures.contrasts.groups")
    named = dict(zip(samples.GROUPS, roles, strict=True))
    shown = {c: label(labels[c]) for c in tests}
    rows = values.frame.filter((pl.col("gene_id") == card.gene_id)
                               & pl.col("contrast").is_in(list(tests)))  # fmt: skip
    if family.batch is not None:
        rows = rows.filter(pl.col("batch") == family.batch)
    frame = rows.select(feature=pl.col("contrast").replace_strict(shown),
                        group=pl.col("group").replace_strict(named),
                        contrast=pl.col("contrast").replace_strict(dict(labels)),
                        sample="sample", donor="donor", batch="batch", value="value")  # fmt: skip
    paired = {values.paired.get(c, False) for c in tests}
    top = max(tests, key=lambda c: (abs(tests[c].r), -tests[c].p_adj))
    names = list(family.tests)
    numbers = {"symbol": card.symbol, "family": len(family.tests), "left": len(family.left),
               "view": family.name, "top": labels[top], "r": tests[top].r,
               "p": _p(tests[top].p_adj), "alpha": ALPHA,
               "below": sum(t.p_adj < ALPHA for t in family.tests.values()),
               "first": names.index(drawn[0]) + 1, "last": names.index(drawn[-1]) + 1}  # fmt: skip
    caption = OwnCaption(
        (roles[0], roles[1]),
        tuple(FeatureTest(shown[c], t.n_reference, t.n_tested, t.statistic, t.p, t.p_adj, t.r)
              for c, t in tests.items()),
        words=text.fill(PAGE, "figures.contrasts.caption", **numbers), with_n=True,
    )  # fmt: skip
    data = OwnBox(frame=frame, columns=_columns("contrasts"),
                  feature_label=text.word(PAGE, "figures.contrasts.x"),
                  value_label=text.word(PAGE, "figures.contrasts.y"),
                  unit=_columns("contrasts")["value"][0], groups=(roles[0], roles[1]),
                  hover=("contrast", "sample", "donor", "batch"), caption=caption,
                  stats=_stats(True in paired, False in paired))  # fmt: skip
    return data, numbers


def families(card: Card, values: samples.Values, contrasts_: Sequence[str],
             levels: Sequence[str]) -> list[Family]:  # fmt: skip
    """F1's views' families: every sample first, then each batch level; two tests or more each.

    Empty where every sample together tests fewer than two contrasts: the box shows the one.
    """
    out = []
    for batch in (None, *levels):
        tests, left = samples.across(values, card.gene_id, contrasts_, batch)
        if len(tests) < 2:  # noqa: PLR2004 - one contrast is the box
            if batch is None:
                return []
            continue
        name = (text.word(PAGE, "figures.contrasts.views.all") if batch is None
                else text.fill(PAGE, "figures.contrasts.views.batch", batch=batch))  # fmt: skip
        out.append(Family(name, batch, tests, left))
    return out


def views(card: Card, values: samples.Values, labels: Mapping[str, str],
          levels: Sequence[str]) -> list[tuple[str, OwnBox, dict[str, Any]]]:  # fmt: skip
    """(label, figure, numbers) of every F1 view: each family, a figure per page of contrasts.

    The first, every sample's first page, is the figure as built, which the files draw.
    """
    out = []
    for family in families(card, values, list(labels), levels):
        for drawn in pages(list(family.tests)):
            data, numbers = contrasts(card, values, labels, family, drawn)
            out.append((text.fill(PAGE, "figures.contrasts.views.page", **numbers), data, numbers))
    return out


def table(card: Card, values: samples.Values, labels: Mapping[str, str],
          levels: Sequence[str]) -> pl.DataFrame:  # fmt: skip
    """Every F1 view's tests for the data file: past the first, a view has no files of its own."""
    rows = [(card.symbol, card.gene_id, f.name, labels[c], t.n_reference, t.n_tested,
             _test_words(values.paired.get(c, False)), t.statistic, t.p, t.p_adj, t.r)
            for f in families(card, values, list(labels), levels)
            for c, t in f.tests.items()]  # fmt: skip
    return pl.DataFrame(rows, schema=list(_columns("contrast_tests")), orient="row")
