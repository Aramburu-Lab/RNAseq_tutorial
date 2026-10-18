"""The one reader of a protein (or gene) by sample table, beside its sample table.

The dataset file's ``[proteo]`` declares what the table holds: the values as intensities or as
logarithms of a named base, the tokens written for a missing value, the identifier column and
its kind, and the sample table pairing each protein column with an RNA sample and naming its
TMT plex or label-free run. An adapter (``ctx.proteo.adapters``) turns a tool's own export into
the text columns this reader takes. Everything wrong is listed together before any value is
used. The order is fixed and tested: the adapter's row drops, the values parsed (a declared
token is missing, never zero; an undeclared zero or text is a problem), the identifiers mapped
(``ctx.proteo.identifiers``), then one row per gene (``tables.one_per_gene``).
"""

import io
import math
import zipfile
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import polars as pl

from ctx.dataset import ProteoTable
from ctx.gather.own_fastq import table_rows
from ctx.proteo import ProteoError, identifiers, tables

# A logarithm's base to the factor that turns it into log2: log2(x) = log_b(x) * log2(b).
LOG_BASES = {2: 1.0, 10: math.log2(10), "e": math.log2(math.e)}
# A problem list shows this many examples of one kind.
EXAMPLES = 3


@dataclass(frozen=True)
class Text:
    """A table as an adapter hands it over: the identifier and sample columns, as text.

    Attributes:
        frame: One column per identifier and sample, every cell a string, "" for an empty one.
        dropped: (reason, rows) of each row drop the adapter made, in its order.
        columns: How many sample columns the export held, before the sample table chose.
    """

    frame: pl.DataFrame
    dropped: tuple[tuple[str, int], ...] = ()
    columns: int = 0


@dataclass(frozen=True)
class SampleTable:
    """The sample table: one row per protein column, in its order.

    Attributes:
        protein: Protein column names.
        rna: Each one's RNA sample; "" for one with no RNA.
        batch: Each one's plex or run.
        tracks: Every further column, for the page's strips.
        digest: SHA-256 of the file as saved.
        reference: Per row, whether it is its plex's reference channel; empty for none.
    """

    protein: tuple[str, ...]
    rna: tuple[str, ...]
    batch: tuple[str, ...]
    tracks: pl.DataFrame
    digest: str
    reference: tuple[bool, ...] = ()

    def without_references(self) -> "SampleTable":
        """The rows that are not a reference channel, in order."""
        keep = [i for i in range(len(self.protein)) if not (self.reference and self.reference[i])]
        return SampleTable(
            tuple(self.protein[i] for i in keep),
            tuple(self.rna[i] for i in keep),
            tuple(self.batch[i] for i in keep),
            self.tracks[keep],
            self.digest,
        )


@dataclass(frozen=True)
class Read:
    """The reader's result: the layer and every count of what it left out.

    Attributes:
        layer: One row per RNA gene, the sample table's protein columns, log2, NaN missing.
        dropped: The adapter's row drops.
        mapped: Rows per mapping reason (``identifiers.Mapped.counts``).
        several: (identifier cell, the genes it names) of each row naming several genes.
        tokens: Cells read as missing, per declared token.
        columns: Sample columns of the export.
        samples: The sample table's rows the layer holds: every row but the reference channels.
    """

    layer: tables.Layer
    dropped: tuple[tuple[str, int], ...]
    mapped: dict[str, int]
    several: tuple[tuple[str, str], ...]
    tokens: dict[str, int] = field(default_factory=dict)
    columns: int = 0
    samples: SampleTable | None = None


def raw_bytes(path: Path, member: str | None) -> bytes:
    """The table's bytes: the file itself, or one member of a tar or zip archive by base name.

    Raises:
        ProteoError: If the file is missing, or the archive holds the member other than once.
    """
    if not path.is_file():
        raise ProteoError(
            f"{path.name} is missing; run ctx setup with the dataset file, or put "
            "the table where the dataset file names it"
        )
    if member is None:
        return path.read_bytes()
    if not zipfile.is_zipfile(path):
        return tables._member(path, member)  # noqa: SLF001 - the one tar reader of M7
    with zipfile.ZipFile(path) as archive:
        found = [n for n in archive.namelist() if Path(n).name == member and not n.endswith("/")]
        if len(found) != 1:
            raise ProteoError(f"{path.name}: {len(found)} members named {member}")
        return archive.read(found[0])


def text_frame(
    data: bytes, separator: str, columns: list[str] | None = None, *, quoted: bool = True
) -> pl.DataFrame:
    """A tab- or comma-separated table, every cell a string and an empty one "".

    A plain table honours quotes, as a spreadsheet or CPTAC's own tables write them; MaxQuant
    and DIA-NN quote nothing, so their adapters read with ``quoted`` off and keep a cell that
    starts with a quote as written.

    Raises:
        ProteoError: If the text is not UTF-8 or a row's width differs from the header's.
    """
    try:
        frame = pl.read_csv(
            io.BytesIO(data),
            separator=separator,
            infer_schema=False,
            quote_char='"' if quoted else None,
            columns=columns,
        )
    except (pl.exceptions.ComputeError, UnicodeDecodeError) as error:
        raise ProteoError(f"the protein table cannot be read as text: {error}") from error
    return frame.with_columns(pl.all().fill_null(""))


def sample_table(path: Path, declared: ProteoTable, problems: list[str]) -> SampleTable | None:
    """The sample table beside the dataset file, every problem added to ``problems``."""
    if not path.is_file():
        problems.append(f"proteo.samples: {path.name} is not beside the dataset file")
        return None
    rows, digest = table_rows(path)
    header, body = (rows[0], rows[1:]) if rows else ([], [])
    wanted = tuple(c for c in (declared.sample_column, declared.rna_column, declared.batch,
                               declared.reference) if c)  # fmt: skip
    if missing := [c for c in wanted if c not in header]:
        problems.append(f"{path.name}: no column {missing}; its header is {header}")
        return None
    if bad := [i for i, r in enumerate(body, 2) if len(r) != len(header)]:
        problems.append(f"{path.name}: line {bad[0]} has not the header's {len(header)} cells")
        return None
    frame = pl.DataFrame(body, schema=header, orient="row")
    protein, rna = frame[declared.sample_column].to_list(), frame[declared.rna_column].to_list()
    for name, values in ((declared.sample_column, protein), (declared.rna_column, rna)):
        if twice := sorted(v for v, k in Counter(v for v in values if v).items() if k > 1):
            problems.append(f"{path.name}: {name} {twice[:EXAMPLES]} on several rows")
    if empty := [
        i
        for i, r in enumerate(body, 2)
        if not r[header.index(declared.sample_column)] or not r[header.index(declared.batch)]
    ]:
        problems.append(f"{path.name}: line {empty[0]} lacks its protein sample or batch")
    if not body:
        problems.append(f"{path.name}: no sample rows")
    tracks = frame.drop(*wanted).with_columns(pl.Series("sample", protein))
    flags = tuple(v.strip().lower() == "yes" for v in frame[declared.reference].to_list()
                  ) if declared.reference else ()  # fmt: skip
    return SampleTable(
        tuple(protein), tuple(rna), tuple(frame[declared.batch].to_list()), tracks, digest, flags
    )


def to_reference(values: np.ndarray, samples: SampleTable, problems: list[str]) -> np.ndarray:
    """Each value as a log2 ratio to its plex's one reference channel; NaN where that is missing.

    Raises nothing: a plex with no reference or several is added to ``problems``.
    """
    out = values.copy()
    for batch in dict.fromkeys(samples.batch):
        columns = [i for i, b in enumerate(samples.batch) if b == batch]
        refs = [i for i in columns if samples.reference[i]]
        if len(refs) != 1:
            problems.append(f"plex {batch!r}: {len(refs)} reference channels, not one")
            continue
        out[:, columns] = values[:, columns] - values[:, refs]
    return out


def _column(cells: pl.Series, declared: ProteoTable, problems: list[str]) -> np.ndarray:
    """One sample column as log2, NaN where a declared token stands; problems for the rest."""
    token = cells.is_in(list(declared.missing))
    number = cells.str.strip_chars().cast(pl.Float64, strict=False)
    finite = number.is_not_null() & number.is_finite()
    text = cells.filter(~token & ~finite.fill_null(False))
    if text.len():
        problems.append(
            f"column {cells.name!r}: {text.len():,} cells neither a finite number nor a "
            f"declared missing token, such as {text.unique()[:EXAMPLES].to_list()}"
        )
    values = np.where(token.to_numpy(), np.nan, number.fill_null(np.nan).to_numpy())
    if declared.scale == "log":
        return values * LOG_BASES[declared.log_base]
    if (low := int(np.sum(values <= 0))) > 0:
        problems.append(
            f"column {cells.name!r}: {low:,} intensities at or below 0 that are not "
            f"a declared missing token {list(declared.missing)}"
        )
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(values > 0, np.log2(np.where(values > 0, values, 1.0)), np.nan)


def read(
    text: Text,
    samples: SampleTable,
    declared: ProteoTable,
    genes: pl.DataFrame,
    id_col: str,
    id_map: dict[str, str] | None = None,
) -> Read:
    """Parse, map and collapse a table in that order; every problem raised together.

    Args:
        text: The adapter's text columns.
        samples: The sample table; its protein columns are the layer's, in its order.
        declared: What ``[proteo]`` declares.
        genes: Gather's genes table, the RNA side the identifiers meet.
        id_col: The identifier column (``adapters.id_column``).
        id_map: UniProt accession to symbol, for ``uniprot`` identifiers.

    Returns:
        The layer and the counts of what each step left out.

    Raises:
        ProteoError: Listing every problem found.
    """
    frame, problems = text.frame, []
    if id_col not in frame.columns:
        problems.append(f"no identifier column {id_col!r} in the protein table")
    absent = [s for s in samples.protein if s not in frame.columns]
    if absent:
        problems.append(
            f"{len(absent)} sample table protein samples are not columns of the "
            f"protein table, such as {absent[:EXAMPLES]}"
        )
    if problems:
        raise ProteoError("; ".join(problems))
    values = np.column_stack([_column(frame[s], declared, problems) for s in samples.protein])
    tokens = {t: int(sum((frame[s] == t).sum() for s in samples.protein)) for t in declared.missing}
    used = samples
    if samples.reference:
        values = to_reference(values, samples, problems)
        keep = [i for i, r in enumerate(samples.reference) if not r]
        values, used = values[:, keep], samples.without_references()
    try:
        mapped = identifiers.to_genes(frame[id_col].to_list(), declared.id_kind, genes, id_map)
    except ValueError as error:
        problems.append(str(error))
    if problems:
        raise ProteoError("; ".join(problems))
    kept, report = tables.one_per_gene(list(mapped.genes), values)
    counts = mapped.counts()
    collapse = tables.Collapse(
        report.rows,
        counts[identifiers.NO_ID],
        counts[identifiers.SEVERAL],
        report.extra_rows,
        report.genes_twice,
        report.kept,
    )
    layer = tables.Layer(
        tuple(str(mapped.genes[i]) for i in kept),
        used.protein,
        values.reshape(frame.height, len(used.protein))[kept],
        collapse,
    )
    return Read(layer, text.dropped, counts, mapped.several, tokens, text.columns, used)
