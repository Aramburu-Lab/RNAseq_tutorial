"""The declared sample annotation of an across-samples dataset, joined onto gather's sheet.

A GDC dataset with an across-samples ``[proteo]`` table may declare ``[columns.<name>]``: each
column's ``source`` is a column annotation of the protein table (for CPTAC breast, the freeze's
Table S1A), its ``values`` the curated spelling. ``ctx gather --gdc-dir`` writes them onto the
sheet, so ``ctx de`` fits a contrast on them and M7's F4 compares the same tumours. Tumours meet
annotation columns by ``pair_samples``, the one pairing M7 also measures on.
"""

from pathlib import Path

import polars as pl

from ctx.dataset import Column, Dataset, DatasetError
from ctx.gather import gdc
from ctx.proteo import ProteoError, tables


def pair_samples(names: tuple[str, ...], sheet: pl.DataFrame, prefix: str
                 ) -> tuple[list[tuple[int, str]], int]:  # fmt: skip
    """Each protein column's chosen tumour: by its case after ``prefix``, else by GDC sample ID.

    The CPTAC freeze names most columns by case and two by GDC's sample submitter ID, so a join
    by case alone left both out (TODO.md, M7: the two freeze columns named by aliquot).

    Args:
        names: The protein table's sample names, in column order.
        sheet: Gather's chosen tumours, one per case, with ``case``, ``sample_id`` and ``sample``.
        prefix: Text before a case ID in the protein names.

    Returns:
        (protein column, sheet sample) per paired column in column order, and how many were
        paired by sample ID.

    Raises:
        ProteoError: If no column pairs, or two columns pair with one tumour.
    """
    by_case = dict(sheet.select("case", "sample").rows())
    by_id = dict(sheet.select("sample_id", "sample").rows())
    pairs, by_sample_id = [], 0
    for j, name in enumerate(names):
        case = name.removeprefix(prefix)
        if case in by_case:
            pairs.append((j, by_case[case]))
        elif name in by_id:
            pairs.append((j, by_id[name]))
            by_sample_id += 1
    if not pairs:
        raise ProteoError("no protein sample's case has a chosen tumour in gather's sheet")
    if len({s for _, s in pairs}) < len(pairs):
        raise ProteoError("two protein columns pair with one tumour of gather's sheet")
    return pairs, by_sample_id


def _spelled(column: Column, raw: list[str], problems: list[str]) -> list[str | None]:
    """One annotation column in the curated spelling; an unmapped value is a problem."""
    if column.values is None:
        return [v or None for v in raw]
    unmapped = sorted({v for v in raw if v not in column.values})
    problems += [f"columns.{column.name}: {v!r} is not in values" for v in unmapped]
    return [column.values.get(v) for v in raw]


def sheet_columns(sheet: pl.DataFrame, spec: Dataset, root: Path) -> pl.DataFrame:
    """Gather's sheet with every declared column, filled on the chosen tumours paired.

    Order: the chosen tumours (``gdc.chosen``), then ``pair_samples`` against the protein
    table's columns, then each value mapped to its curated spelling. Every other row is null,
    so no normal file of a case takes its tumour's call.

    Raises:
        DatasetError: If the dataset declares columns without an across-samples ``[proteo]``.
        ProteoError: If the protein table is missing, lacks a source annotation, or holds a
            value the declared spelling does not map; every problem is listed.
    """
    p = spec.proteo
    if p is None or spec.gdc is None or p.design != "across_samples" or not p.member:
        raise DatasetError(f"{spec.path.name}: a GDC dataset's [columns] are read from an "
                           "across-samples [proteo] table's annotation")  # fmt: skip
    layer, annotation = tables.read_gct(root / spec.resource(p.protein).dest, p.member)
    problems = [f"columns.{c.name}: no annotation {c.source!r} in {p.member}"
                for c in spec.columns if c.source not in annotation.columns]  # fmt: skip
    if problems:
        raise ProteoError("; ".join(problems))
    pairs, _ = pair_samples(layer.samples, gdc.chosen(sheet, spec.gdc.tumour), p.sample_prefix)
    added = {}
    for column in spec.columns:
        spelled = _spelled(column, annotation[column.source].to_list(), problems)
        found = {s: spelled[j] for j, s in pairs}
        added[column.name] = pl.col("sample").replace_strict(found, default=None,
                                                             return_dtype=pl.String)  # fmt: skip
    if problems:
        raise ProteoError("; ".join(problems))
    return sheet.with_columns(**added)
