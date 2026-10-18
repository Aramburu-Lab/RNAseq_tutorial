"""M7 across samples on a protein table of any study: ``[proteo.table]`` read beside RNA-seq.

RNA is gather's counts of the samples the sample table pairs, as log2(CPM + 1), from any route
gather takes; protein is ``quantities.read`` through the adapter ``[proteo.table] format``
names. The result is the across-samples design's ``Measured``, so the page, the figures, F4 and
the rules are the same as for the CPTAC breast freeze; only the reader and the pairing differ.
"""

import gzip
from pathlib import Path

import numpy as np
import polars as pl

from ctx import runtime
from ctx.dataset import Dataset, ProteoTable
from ctx.de.filters import log2cpm
from ctx.de.run import newest
from ctx.gather import gdc
from ctx.io.naming import MODULES
from ctx.proteo import ProteoError, adapters, identifiers, quantities
from ctx.proteo.figures import COLOURS
from ctx.proteo.measure import per_gene
from ctx.proteo.run import GATHER, Measured, Options

SEVERAL_COLUMNS = {
    "group": ("none", "The protein table's identifier cell: a protein group's members."),
    "genes": ("none", "The distinct genes those members name, joined by ';'; none was chosen."),
}
# The page line on each acquisition, named whole so the content test finds both keys.
DESIGN_NOTES = {"tmt": "table_tmt", "label_free": "table_label_free"}
BATCH_COLUMNS = {
    "batch": ("none", "A TMT plex or a label-free run, from the sample table."),
    "samples": ("samples", "Protein samples in it."),
    "paired": ("samples", "Of those, samples with an RNA sample in gather's counts."),
}


def protein_path(settings: runtime.Settings, spec: Dataset) -> Path:
    """The protein table: beside the dataset file (``table.file``) or where setup put it."""
    p = spec.proteo
    assert p is not None and p.table is not None  # noqa: S101 - run() dispatches here on it
    return (
        spec.path.parent / p.table.file
        if p.table.file
        else settings.root / spec.resource(p.protein).dest
    )


def id_map(settings: runtime.Settings, spec: Dataset, table: ProteoTable) -> dict[str, str]:
    """UniProt accession to gene symbol from the declared map, gzip or not; empty if none.

    Raises:
        ProteoError: If the map is missing or lacks a declared column.
    """
    if not table.id_map:
        return {}
    path = settings.root / spec.resource(table.id_map).dest
    if not path.is_file():
        raise ProteoError(f"{path.name} is missing; run ctx setup --dataset-file {spec.path}")
    data = path.read_bytes()
    data = gzip.decompress(data) if path.suffix == ".gz" else data
    separator = "," if path.name.removesuffix(".gz").endswith(".csv") else "\t"
    frame = quantities.text_frame(data, separator)
    accession, symbol = table.id_map_columns
    if missing := [c for c in (accession, symbol) if c not in frame.columns]:
        raise ProteoError(f"{path.name}: no column {missing}")
    # UniProt lists one primary symbol per entry; an entry with none maps nowhere.
    return {a: s for a, s in frame.select(accession, symbol).rows() if a and s and " " not in s}


def text_of(
    path: Path, member: str | None, table: ProteoTable, samples: quantities.SampleTable
) -> tuple[quantities.Text, str]:
    """The adapter's text columns and the identifier column it read."""
    data = quantities.raw_bytes(path, member)
    id_col = adapters.id_column(table.format, table.id_kind, table.id_column or None)
    if table.format == "maxquant":
        return adapters.maxquant(data, id_col, table.quantity), id_col
    if table.format == "diann":
        return adapters.diann(data, id_col), id_col
    separator = "," if (member or path.name).endswith(".csv") else "\t"
    return adapters.table(data, separator, id_col, samples.protein), id_col


def library_sizes(sheet: pl.DataFrame, counts: pl.DataFrame, samples: list[str]) -> np.ndarray:
    """Each sample's library size: the sheet's where it has one (GDC), else its counts' sum."""
    if "library_size" in sheet.columns:
        sizes = dict(sheet.select("sample", "library_size").rows())
        return np.array([sizes[s] for s in samples])
    return counts.select(samples).sum().to_numpy().ravel()


def _batches(samples: quantities.SampleTable, paired: set[str]) -> pl.DataFrame:
    """One row per batch: its samples and how many have RNA, in first-seen order."""
    frame = pl.DataFrame({"batch": list(samples.batch), "sample": list(samples.protein)})
    return (
        frame.group_by("batch", maintain_order=True)
        .agg(samples=pl.len(), paired=pl.col("sample").is_in(list(paired)).sum())
        .with_columns(pl.col("samples", "paired").cast(pl.Int64))
    )


def across_table(settings: runtime.Settings, spec: Dataset, options: Options) -> Measured:
    """One cohort at both layers from a ``[proteo.table]``: per-gene correlation across samples.

    Only the samples the sample table pairs with an RNA sample enter a correlation.

    Raises:
        ProteoError: Listing every problem of the sample table, the protein table and the
            pairing together.
    """
    p = spec.proteo
    assert p is not None and p.table is not None  # noqa: S101 - run() dispatches here on it
    table = p.table
    gather = settings.root / "outputs" / f"{dict(MODULES)['gather']:02d}_gather"
    found = {w: pl.read_parquet(newest(gather, spec.slug, w)) for w in GATHER}
    sheet = gdc.chosen(found["samplesheet"]).unique("sample", keep="first", maintain_order=True)
    problems: list[str] = []
    samples = quantities.sample_table(spec.path.parent / table.samples, table, problems)
    if samples is not None:
        problems += _pairing(samples, table.samples, found["counts"], sheet)
        named = [c for c in (*p.tracks, p.amplified) if c and c not in samples.tracks.columns]
        problems += (
            [f"proteo: tracks or amplified name {named}, not sample table columns"] if named else []
        )
    if problems or samples is None:
        raise ProteoError("; ".join(problems))
    path = protein_path(settings, spec)
    text, id_col = text_of(path, p.member, table, samples)
    read = quantities.read(
        text, samples, table, found["genes"], id_col, id_map(settings, spec, table)
    )
    protein = read.layer
    assert read.samples is not None  # noqa: S101 - read() always sets it
    samples = read.samples
    pairs = [(j, r) for j, r in enumerate(samples.rna) if r]
    if len(pairs) < 3:  # noqa: PLR2004 - a correlation needs three samples (measure.MIN_PAIRS)
        raise ProteoError(f"{table.samples}: {len(pairs)} protein samples have an RNA sample")
    rna_samples = [r for _, r in pairs]
    rna_all = log2cpm(
        found["counts"].select(rna_samples).to_numpy(),
        library_sizes(sheet, found["counts"], rna_samples),
    )
    rows = {n: i for i, n in enumerate(identifiers.rna_names(found["genes"]))}
    genes = [g for g in protein.genes if g in rows]
    rna = rna_all[[rows[g] for g in genes]]
    at = protein.row()
    values = protein.values[[at[g] for g in genes]][:, [j for j, _ in pairs]]
    frame = per_gene(tuple(genes), rna, values, options.min_quantified)
    annotation = samples.tracks.with_columns(pl.Series(p.plex, list(samples.batch)))
    out = Measured(protein, list(samples.batch), frame, rna.mean(axis=1))
    batches = _batches(samples, {samples.protein[j] for j, _ in pairs})
    out.notes += _notes(read, table, batches, (len(pairs), len(samples.protein), path.name))
    out.extra |= {
        "rna": rna,
        "values": values,
        "n_samples": len(pairs),
        "annotation": annotation[[j for j, _ in pairs]],
        "pairs": pairs,
        "sheet": sheet,
        "inputs": {"protein table": path.name, "sample table": f"{table.samples} {samples.digest}"},
        "tables": [
            (
                "protein_groups_several_genes",
                pl.DataFrame(list(read.several), schema=list(SEVERAL_COLUMNS), orient="row"),
                SEVERAL_COLUMNS,
            ),
            ("batches", batches, BATCH_COLUMNS),
        ],
    }
    out.extra["colour"] = colour(annotation, [*p.tracks, p.plex], table.acquisition)
    return out


def colour(annotation: pl.DataFrame, columns: list[str | None], fallback: str) -> list[str]:
    """F1's bar colours: the first column with no more levels than data colours, else one level.

    A label-free study of one run per sample has as many batches as samples, more than F1's
    colours (figures.COLOURS); its bars then share the acquisition's name, the run in the hover.
    """
    for name in columns:
        if name and annotation[name].n_unique() <= len(COLOURS):
            return [str(v) for v in annotation[name]]
    return [fallback] * annotation.height


def _pairing(
    samples: quantities.SampleTable, name: str, counts: pl.DataFrame, sheet: pl.DataFrame
) -> list[str]:
    """Why the sample table's RNA samples cannot be read, or nothing.

    Each must be a column of gather's counts and a sample of its sheet as ``gdc.chosen`` passes it.
    """
    known = set(counts.columns[1:]) & set(sheet["sample"])
    unknown = [r for r in samples.rna if r and r not in known]
    if not unknown:
        return []
    return [
        f"{name}: {len(unknown)} RNA samples are not in gather's counts and sheet, such as "
        f"{unknown[: quantities.EXAMPLES]}"
    ]


def _notes(
    read: quantities.Read, table: ProteoTable, batches: pl.DataFrame, n: tuple[int, int, str]
) -> list[tuple[str, dict[str, object]]]:
    """The page's lines on the reader: drops, tokens, mapping, the batch design, the pairing."""
    scale = table.scale if table.scale == "linear" else f"log{table.log_base}"
    tokens = ", ".join(f"{t!r} {k:,}" for t, k in read.tokens.items())
    out: list[tuple[str, dict[str, object]]] = [
        ("table_read", {"file": n[2], "format": table.format, "columns": read.columns,
                        "scale": scale, "tokens": tokens}),
    ]  # fmt: skip
    if table.format == "maxquant":
        dropped = dict(read.dropped)
        out.append(("table_maxquant", {k: dropped.get(k, 0) for k in
                                       ("reverse", "contaminant", "only_by_site")}))  # fmt: skip
    sizes = [int(v) for v in batches["samples"]]
    out += [
        ("table_mapped", {"kind": table.id_kind, **read.mapped}),
        (DESIGN_NOTES[table.acquisition], {"n": len(sizes), "low": min(sizes), "high": max(sizes)}),
        ("table_paired", {"n": n[0], "n_columns": n[1]}),
    ]
    if table.reference:
        out.insert(-1, ("table_reference", {"column": table.reference}))
    return out
