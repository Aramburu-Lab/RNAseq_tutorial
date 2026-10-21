"""M9's own samples (E.1c2): the box of the featured contrast, and F1, the gene in every contrast.

A value is M5's normalised value: a contrast's size factors as M5 wrote them
(``<date>_de_<dataset>_size_factors_v1``) over gather's counts, through
``ctx.de.normalised.log2_normalised``, the one definition (D-083). M9 fits no size factor of
its own, so a sample in two contrasts has a value in each, and a sample in none has none.

The tests follow the box convention (TODO.md, Figure and statistics conventions): Mann-Whitney
U for an unpaired contrast, Wilcoxon signed-rank over individuals for a paired one, each with
its rank-biserial r, BH over the family the caption names.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import false_discovery_control, mannwhitneyu, rankdata, wilcoxon

from ctx.dataset import Dataset, DatasetError
from ctx.de.normalised import WHAT as SIZE_FACTORS
from ctx.de.normalised import log2_normalised
from ctx.de.results import newest
from ctx.io.naming import MODULES
from ctx.report import text

# Fewer values than this in a group, or pairs in a paired contrast, leave it untested.
MIN_PER_GROUP = 2
GROUPS = ("reference", "tested")


@dataclass(frozen=True)
class Values:
    """The listed genes' normalised values in every contrast M5 wrote size factors for.

    Attributes:
        frame: ``gene_id, contrast, sample, donor, group, batch, value``; ``group`` is
            reference or tested, as the fit coded the compared levels; ``batch`` is the
            sample sheet's, empty where the sheet has none.
        paired: Each contrast's name to whether its test pairs by individual.
    """

    frame: pl.DataFrame
    paired: Mapping[str, bool]


@dataclass(frozen=True)
class Test:
    """One gene's comparison of a contrast's two groups, before and after BH."""

    n_reference: int
    n_tested: int
    statistic: float
    p: float
    r: float
    p_adj: float = float("nan")


def _out(root: Path, module: str) -> Path:
    """A module's output folder under the root."""
    return root / "outputs" / f"{dict(MODULES)[module]:02d}_{module}"


def read(root: Path, spec: Dataset, gene_ids: Sequence[str]) -> Values | str:
    """Every listed gene's value in each contrast's samples; why not, where a table is missing.

    Only the listed genes' rows of the counts are read: Step-1's table holds 78,470 genes.
    The reason names the table and the module that writes it, never a path: a card printed
    the root's absolute path, which ran past a phone's page (and a kit carries no root).
    """
    found = {}
    for module, what in (("de", SIZE_FACTORS), ("gather", "counts"), ("gather", "samplesheet")):
        try:
            found[what] = newest(_out(root, module), spec.slug, what)
        except DatasetError:
            return text.fill("genecard", "gene.no_table", what=what, module=module)
    factors = pl.read_parquet(found[SIZE_FACTORS])
    counts_path, sheet_path = found["counts"], found["samplesheet"]
    wanted = ["gene_id", *factors["sample"].unique(maintain_order=True)]
    counts = (pl.scan_parquet(counts_path).filter(pl.col("gene_id").is_in(list(gene_ids)))
              .select(wanted).collect())  # fmt: skip
    sheet = pl.read_parquet(sheet_path)
    batch = sheet.select("sample", pl.col("batch").cast(pl.String).fill_null("")) if (
        "batch" in sheet.columns) else sheet.select("sample", batch=pl.lit(""))  # fmt: skip
    frames = []
    for (name,), one in factors.group_by("contrast", maintain_order=True):
        levels = log2_normalised(counts, one, gene_ids)
        frames.append(
            levels.join(one.select("sample", "donor", "group"), on="sample", how="left",
                        maintain_order="left")
            .with_columns(contrast=pl.lit(name))
        )  # fmt: skip
    frame = pl.concat(frames).join(batch, on="sample", how="left", maintain_order="left")
    frame = frame.with_columns(pl.col("batch").fill_null("")).select(
        "gene_id", "contrast", "sample", "donor", "group", "batch", "value"
    )
    paired = {c.name: c.paired for c in spec.contrasts}
    return Values(frame, {c: paired.get(c, False) for c in frame["contrast"].unique()})


def rank_biserial(differences: np.ndarray) -> float:
    """Matched-pairs rank-biserial r: rank sums of positive minus negative differences, over all.

    Zero differences are dropped, as scipy's signed-rank test drops them by default.
    """
    d = differences[differences != 0]
    ranks = rankdata(np.abs(d))
    return float((ranks[d > 0].sum() - ranks[d < 0].sum()) / ranks.sum())


def compare(rows: pl.DataFrame, paired: bool) -> Test | None:
    """The two groups of one gene in one contrast compared; None where no test can rank them.

    Unpaired: Mann-Whitney U of tested against reference, r = 2U / (n1 n2) - 1. Paired: the
    individuals with a value in both groups, Wilcoxon signed-rank on tested minus reference.
    None with fewer than ``MIN_PER_GROUP`` values or pairs, or with nothing to rank: every
    value equal, or every difference zero.
    """
    if paired:
        both = (rows.filter(pl.col("group") == "reference").select("donor", ref="value")
                .join(rows.filter(pl.col("group") == "tested").select("donor", tst="value"),
                      on="donor"))  # fmt: skip
        d = (both["tst"] - both["ref"]).to_numpy()
        if both.height < MIN_PER_GROUP or not np.any(d != 0):
            return None
        result = wilcoxon(both["tst"].to_numpy(), both["ref"].to_numpy())
        return Test(both.height, both.height, float(result.statistic), float(result.pvalue),
                    rank_biserial(d))  # fmt: skip
    ref = rows.filter(pl.col("group") == "reference")["value"].to_numpy()
    tst = rows.filter(pl.col("group") == "tested")["value"].to_numpy()
    if min(ref.size, tst.size) < MIN_PER_GROUP or np.ptp(np.concatenate([ref, tst])) == 0:
        return None
    u, p = mannwhitneyu(tst, ref, alternative="two-sided")
    return Test(ref.size, tst.size, float(u), float(p), 2 * float(u) / (ref.size * tst.size) - 1)


def adjusted(tests: Mapping[str, Test]) -> dict[str, Test]:
    """The tests with BH-adjusted p over all of them, in their order."""
    if not tests:
        return {}
    p_adj = false_discovery_control([t.p for t in tests.values()], method="bh")
    return {k: Test(t.n_reference, t.n_tested, t.statistic, t.p, t.r, float(q))
            for (k, t), q in zip(tests.items(), p_adj, strict=True)}  # fmt: skip


def family(values: Values, contrast: str, gene_ids: Sequence[str]) -> dict[str, Test]:
    """The box's family: each listed gene with a test in ``contrast``, BH over them, list order."""
    one = values.frame.filter(pl.col("contrast") == contrast)
    raw = {}
    for gene in gene_ids:
        found = compare(one.filter(pl.col("gene_id") == gene), values.paired.get(contrast, False))
        if found is not None:
            raw[gene] = found
    return adjusted(raw)


def across(values: Values, gene_id: str, contrasts: Sequence[str],
           batch: str | None = None) -> tuple[dict[str, Test], list[str]]:  # fmt: skip
    """F1's tests: the gene in each contrast, within one batch or all, BH over those tested.

    Returns:
        The tests by contrast, in ``contrasts``' order, and the contrasts left untested.
    """
    rows = values.frame.filter(pl.col("gene_id") == gene_id)
    if batch is not None:
        rows = rows.filter(pl.col("batch") == batch)
    raw, left = {}, []
    for name in contrasts:
        found = compare(rows.filter(pl.col("contrast") == name), values.paired.get(name, False))
        if found is None:
            left.append(name)
        else:
            raw[name] = found
    return adjusted(raw), left


def batches(values: Values, most: int) -> list[str]:
    """The sheet's batch levels among the contrasts' samples, when there are 2 to ``most``."""
    found = sorted(set(values.frame["batch"].to_list()) - {""})
    return found if 2 <= len(found) <= most else []  # noqa: PLR2004 - a choice needs two
