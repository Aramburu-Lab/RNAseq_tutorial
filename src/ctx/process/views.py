"""M3's agreement figures with the views a reader looks up (WO8, D-051, D-073).

The counter scatter (F3) can show one gene type, and one sample instead of every sample's
sums; the agreement plot (F4) one gene type, its limits as bias plus and minus 1.96 SD or as
d's empirical percentiles. Each view's numbers are ``concordance.per_unit`` on its genes, so
a near-1 rho over every gene cannot hide a gene type that disagrees. Those rows depend on the
comparison and the gene type alone, so ``Agreed`` computes them once per type and every view
of either figure reads them: recomputed per view, the bootstrap ran 11 times over 103 samples
of Breast-tissue-BRCA, 2,579 of a 3,300 s profile (2026-10-24). A sample axis is offered up to
``SAMPLES_MAX`` samples: each one's view carries its genes again, which 210 runs could not.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import polars as pl

from ctx import concordance
from ctx.process import figures
from ctx.process import tables as tab
from ctx.report import text
from ctx.viz.scatter import ScatterData

PAGE = "process"
# Every gene first, as the build draws them; then each type the colours name.
TYPES = (None, *tab.GENE_TYPES, tab.OTHER)
# The most samples the scatter's sample axis offers, each a view of its own (WO8).
SAMPLES_MAX = 8
# The fewest genes with a read a view of one type needs: three ranks make a correlation.
MIN_GENES = 3


@dataclass(frozen=True)
class Subset:
    """One view's genes: a gene type or every gene, and one sample or every sample's sums."""

    gene_type: str | None
    sample: str

    def words(self) -> str:
        """How the live line and the caption name the subset."""
        every = self.sample == concordance.ALL
        unit = (text.word(PAGE, "figures.counters.every_sample") if every
                else text.fill(PAGE, "figures.counters.one_sample",
                               sample=self.sample))  # fmt: skip
        return text.fill(PAGE, "figures.counters.subset", genes=type_words(self.gene_type),
                         unit=unit)  # fmt: skip


def type_words(gene_type: str | None) -> str:
    """A gene type as the page names it; None is every gene."""
    return (text.word(PAGE, "figures.counters.every_type") if gene_type is None
            else text.word(PAGE, f"gene_types.{gene_type}"))  # fmt: skip


def typed(long: pl.DataFrame, sums: pl.DataFrame) -> pl.DataFrame:
    """The comparison table with each gene's name, type and type group as the sums carry them.

    ``tables.agreements`` gave the sums their types once; the comparison table, which the
    preview builds without them, takes them from there.
    """
    names = sums.select("gene_id", "gene_name", "gene_type", "type_group").unique(
        "gene_id", maintain_order=True)  # fmt: skip
    plain = long.select("comparison", "sample", "gene_id", "a", "b")
    return plain.join(names, on="gene_id", how="left", maintain_order="left")


def genes_of(part: pl.DataFrame, gene_type: str | None) -> pl.DataFrame:
    """A comparison's rows of one gene type; None is every gene."""
    return part if gene_type is None else part.filter(pl.col("type_group") == gene_type)


@dataclass(frozen=True)
class Agreed:
    """One comparison's agreement rows per gene type, each ``per_unit`` once, when first asked.

    Attributes:
        part: One comparison's rows of ``typed``.
        rows: Rows already computed, by gene type; None is every gene, which
            ``tables.agreements`` computed for the concordance table (``from_table``).
    """

    part: pl.DataFrame
    rows: dict[str | None, pl.DataFrame] = field(default_factory=dict)

    def of(self, gene_type: str | None) -> pl.DataFrame:
        """The concordance rows over one type's genes, each sample and ``all``."""
        if gene_type not in self.rows:
            genes = genes_of(self.part, gene_type).select("gene_id", "sample", "a", "b")
            self.rows[gene_type] = concordance.per_unit(genes, "a", "b", unit="sample",
                                                        gene="gene_id")  # fmt: skip
        return self.rows[gene_type]


def from_table(table: pl.DataFrame, name: str) -> pl.DataFrame:
    """The concordance table's rows of one comparison as ``per_unit`` returned them."""
    rows = table.filter(pl.col("comparison") == name).drop("comparison")
    return rows.rename({"unit": "sample"})


def measured(part: pl.DataFrame, subset: Subset) -> pl.DataFrame:
    """A subset's genes' counts, typed for the figures.

    Args:
        part: One comparison's rows of ``typed``.
        subset: The view.

    Returns:
        One row per gene: gene_id, a, b, gene_name, gene_type, type_group.
    """
    genes = genes_of(part, subset.gene_type)
    if subset.sample == concordance.ALL:
        counts = genes.group_by("gene_id", maintain_order=True).agg(pl.col("a").sum(),
                                                                    pl.col("b").sum())  # fmt: skip
    else:
        counts = genes.filter(pl.col("sample") == subset.sample).select("gene_id", "a", "b")
    names = genes.select("gene_id", "gene_name", "gene_type", "type_group").unique(
        "gene_id", maintain_order=True)  # fmt: skip
    return counts.join(names, on="gene_id", how="left", maintain_order="left")


def values(rows: pl.DataFrame, subset: Subset, words: Mapping[str, str]) -> dict[str, Any]:
    """The caption's and lead's placeholders for one view."""
    every = rows.filter(pl.col("sample") == subset.sample).row(0, named=True)
    samples = rows.filter(pl.col("sample") != concordance.ALL)
    shown = samples if subset.sample == concordance.ALL else samples.filter(
        pl.col("sample") == subset.sample)  # fmt: skip
    return figures.scatter_values(every, shown, shown.height) | dict(words) | {"every": every}


def types(part: pl.DataFrame, units: Sequence[str]) -> list[str | None]:
    """Every gene, then each type with ``MIN_GENES`` genes with a read in every unit offered.

    Spearman, Pearson and the limits need a few genes; a type that lacks them in one sample
    is not offered at all, so the grid's every combination holds a figure.
    """
    found: list[str | None] = [None]
    for t in TYPES[1:]:
        genes = part.filter(pl.col("type_group") == t)
        sums = genes.group_by("gene_id").agg(pl.col("a").sum(), pl.col("b").sum())
        per = [sums if u == concordance.ALL else genes.filter(pl.col("sample") == u)
               for u in units]  # fmt: skip
        if all(p.filter(pl.col("a") + pl.col("b") > 0).height >= MIN_GENES for p in per):
            found.append(t)
    return found


def subsets(part: pl.DataFrame) -> tuple[list[Subset], list[tuple[str, list[str]]],
                                         list[tuple[int, int]]]:  # fmt: skip
    """The scatter's views, its grid's axes and each view's place, the build's first."""
    samples = part["sample"].unique(maintain_order=True).to_list()
    units = [concordance.ALL, *(samples if len(samples) <= SAMPLES_MAX else [])]
    offered = types(part, units)
    axes = [(text.word(PAGE, "figures.counters.axis_type"), [type_words(t) for t in offered]),
            (text.word(PAGE, "figures.counters.axis_sample"),
             [text.word(PAGE, "figures.counters.every_sample") if u == concordance.ALL else u
              for u in units])]  # fmt: skip
    found = [(Subset(t, u), (i, j)) for i, t in enumerate(offered) for j, u in enumerate(units)]
    return [s for s, _ in found], axes, [p for _, p in found]


def scatter(agreed: Agreed, subset: Subset, words: Mapping[str, str],
            prefix: str) -> tuple[ScatterData, dict[str, Any]]:  # fmt: skip
    """F3 for one view: its contract and its numbers."""
    filled = values(agreed.of(subset.gene_type), subset, words)
    caption = " ".join([prefix + subset.words() + ".",
                        text.fill(PAGE, "figures.counters.caption", **filled)])  # fmt: skip
    data = figures.counter_scatter(measured(agreed.part, subset), words, caption,
                                   figures.limits_of(filled["every"]))  # fmt: skip
    return data, filled


def plot(agreed: Agreed, gene_type: str | None, empirical: bool, words: Mapping[str, str],
         prefix: str) -> tuple[ScatterData, dict[str, Any]]:  # fmt: skip
    """F4 for one view: its contract and its numbers, ``n_out`` beyond the limits drawn."""
    subset = Subset(gene_type, concordance.ALL)
    filled = values(agreed.of(gene_type), subset, words)
    every = filled["every"]
    caption = " ".join([prefix + subset.words() + ".",
                        text.fill(PAGE, "figures.agreement.caption", **filled)])  # fmt: skip
    data = figures.agreement_plot(measured(agreed.part, subset), every, words, caption,
                                  empirical=empirical)  # fmt: skip
    low, high = ("d_p2_5", "d_p97_5") if empirical else ("loa_low", "loa_high")
    outside = text.word(PAGE, "undefined")
    if every[low] is not None:
        beyond = ~pl.col("diff_log2").is_between(every[low], every[high])
        outside = f"{data.frame.filter(beyond).height:,}"
    return data, {**filled, "n_out": outside}
