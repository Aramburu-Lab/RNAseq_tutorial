"""F5 on the page: each tumour's rRNA-depleted library against its other one.

Also the rule's metrics, and, where the sheet carries no read-group record, why they are absent.
"""

from typing import TYPE_CHECKING, Any

import polars as pl

from ctx.de.run import Writer
from ctx.ncrna import biotypes, libraries
from ctx.report import text
from ctx.report.builder import FigureText, Group, Result, ResultTable
from ctx.viz.scatter import Line, ScatterData
from ctx.viz.style import ROLES

if TYPE_CHECKING:
    from ctx.ncrna.run import Tcga

PAGE = "ncrna"
RULE_METRICS = ("library.histone_pairs", "library.histone_higher")
COLUMNS = {
    "sample_id": ("none", "Sample submitter ID: one tumour, two libraries."),
    "sample": ("none", "Aliquot of the library."),
    "library_selection": ("none", "GDC read groups' library selection."),
    "histone": ("%", "Share of counted reads on replication-dependent histone genes."),
    biotypes.SMALL: ("%", "Share of counted reads on small-RNA genes."),
    biotypes.MITO: ("%", "Share of counted reads on mitochondrial rRNA genes."),
    **{f"{c}_{k}": ("%", f"{c} share in the {w} library.") for c in libraries.CLASSES
       for k, w in (("rrna", "rRNA-depleted"), ("other", "other"))},
    "pair": ("none", "One tumour sequenced both ways."),
}  # fmt: skip


def word(key: str) -> str:
    """A fixed text of this page's content file."""
    return text.word(PAGE, key)


def fill(key: str, **values: Any) -> str:
    """A template of this page's content file, filled."""
    return text.fill(PAGE, key, **values)


def group(writer: Writer, tcga: "Tcga", genes: pl.DataFrame, tumour: str) -> tuple[
        Result | None, dict[str, tuple[float, str]], dict[str, str]]:  # fmt: skip
    """F5 and its table, the rule's metrics, and the reason they are absent where they are."""
    if not libraries.measured(tcga.files):
        why = word("tcga.no_read_groups")
        return None, {}, dict.fromkeys(RULE_METRICS, why)
    paired = libraries.pairs(tcga.files, tcga.counts, genes, tumour)
    if paired.is_empty():
        return None, {}, dict.fromkeys(RULE_METRICS, word("tcga.no_pairs"))
    found = libraries.compare(paired)
    n = int(found["histone"]["pairs"])
    metrics = {
        "library.histone_pairs": (float(n), "tumours with both library types"),
        "library.histone_higher": (found["histone"]["higher"],
                                   "of them, rRNA-depleted holds more histone"),
    }  # fmt: skip
    writer.table(paired, "library_pairs", {c: COLUMNS[c] for c in paired.columns})
    frame = libraries.wide(paired).with_columns(pair=pl.lit("tumour"))
    other = ", ".join(sorted(set(paired["library_selection"]) - {libraries.RRNA_DEPLETION}))
    h = found["histone"]
    highest = pl.max_horizontal(pl.col("histone_rrna", "histone_other")).max()
    top = float(frame.select(highest).item())
    caption = fill("figures.libraries.caption", n=n, other=other, higher=int(h["higher"]),
                   p=h["p"])  # fmt: skip
    data = ScatterData(
        frame=frame, columns={c: COLUMNS[c] for c in frame.columns}, x="histone_other",
        y="histone_rrna", x_label=fill("figures.libraries.x", other=other),
        y_label=word("figures.libraries.y"), colour="pair", colours={"tumour": ROLES["comparison"]},
        caption=caption, hover=("sample_id",),
        lines=(Line(word("figures.libraries.equal"), (0.0, top), (0.0, top), dash=True),),
    )  # fmt: skip
    shown = FigureText(title=word("figures.libraries.title"),
                       lead=fill("figures.libraries.lead", higher=int(h["higher"]), n=n,
                                 rrna=h["rrna"], other_median=h["other"]),
                       reading=tuple(line.format(other=other) for line in text.raw(
                           PAGE, "figures.libraries.reading")))  # fmt: skip
    block = writer.figure(data, "library_types", shown)
    rows = tuple((c, f"{v['rrna']:.3f}", f"{v['other']:.3f}", f"{int(v['higher'])} of {n}",
                  f"{v['p']:.2g}") for c, v in found.items())  # fmt: skip
    table = ResultTable(fill("figures.libraries.table_title", other=other),
                        text.words(PAGE, "figures.libraries.columns"), rows,
                        numeric=(1, 2, 4), results=True)  # fmt: skip
    metrics |= {f"library.{c}.p": (v["p"], "Wilcoxon signed-rank, rRNA-depleted against other")
                for c, v in found.items()}  # fmt: skip
    return Group(word("groups.libraries"), (block, table)), metrics, {}
