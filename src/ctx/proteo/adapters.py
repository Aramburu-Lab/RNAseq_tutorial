"""Thin adapters from a tool's own export to the reader's text columns (``quantities.Text``).

Each adapter does only what its tool's format needs: it picks the identifier and sample
columns, names each sample column as the sample table does, and drops the rows its tool marks
as not to be quantified, counting each reason. Values, missing tokens, identifiers and the
collapse are the reader's, the same for every format.

MaxQuant ``proteinGroups.txt``: tab-separated, one quantity per sample named ``<quantity>
<sample>`` (``LFQ intensity``, ``Intensity``, ``iBAQ`` or ``Reporter intensity corrected``).
Rows marked ``+`` in ``Reverse`` (decoys), ``Potential contaminant`` (``Contaminant`` before
MaxQuant 1.6) or ``Only identified by site`` leave, each row counted once under the first of
those three it carries. MaxQuant writes 0 for a protein not quantified in a sample.

DIA-NN ``pg_matrix.tsv``: tab-separated, the descriptor columns then one column per run, named
by the run's file path; a sample is that file's name without its folder and ending. DIA-NN
leaves a protein not quantified in a run empty, and writes decoys and contaminants to no matrix.
"""

import csv
import io
import re
from pathlib import PurePosixPath, PureWindowsPath

import polars as pl

from ctx.proteo import ProteoError
from ctx.proteo.quantities import Text, text_frame

# MaxQuant's row flags, in the order a row is counted under them; the contaminant column was
# renamed in MaxQuant 1.6, so either name is that flag.
MAXQUANT_FLAGS = (
    ("reverse", ("Reverse",)),
    ("contaminant", ("Potential contaminant", "Contaminant")),
    ("only_by_site", ("Only identified by site",)),
)
MAXQUANT_QUANTITIES = ("LFQ intensity", "Intensity", "iBAQ", "Reporter intensity corrected")
# The identifier column of each tool, by the kind the dataset file declares.
MAXQUANT_IDS = {"symbol": "Gene names", "uniprot": "Majority protein IDs"}
DIANN_IDS = {"symbol": "Genes", "uniprot": "Protein.Group"}
# DIA-NN's descriptor columns, every other column being a run (DIA-NN 1.8 to 2.x).
DIANN_DESCRIPTORS = (
    "Protein.Group",
    "Protein.Ids",
    "Protein.Names",
    "Genes",
    "First.Protein.Description",
    "N.Sequences",
    "N.Proteotypic.Sequences",
)
# Run file endings DIA-NN names its columns by.
RUN_ENDINGS = re.compile(r"\.(raw|d|wiff|wiff2|mzml|dia)$", re.IGNORECASE)
FLAGGED = "+"


def id_column(fmt: str, kind: str, declared: str | None) -> str:
    """The identifier column: the dataset file's, else the tool's own for the kind.

    Raises:
        ProteoError: If neither names one.
    """
    own = {"maxquant": MAXQUANT_IDS, "diann": DIANN_IDS}.get(fmt, {})
    if declared or kind in own:
        return declared or own[kind]
    raise ProteoError(f"proteo.id_column: the {fmt} format has no own column for {kind!r}")


def header(data: bytes, separator: str = "\t", *, quoted: bool = False) -> list[str]:
    """The table's first line, split; quotes honoured for a plain table."""
    first = data.split(b"\n", 1)[0].decode("utf-8", errors="replace").rstrip("\r")
    return (
        next(csv.reader(io.StringIO(first), delimiter=separator))
        if quoted
        else first.split(separator)
    )


def table(data: bytes, separator: str, id_col: str, samples: tuple[str, ...]) -> Text:
    """A plain table: the identifier column and the sample table's columns, nothing dropped."""
    names = header(data, separator, quoted=True)
    wanted = [c for c in dict.fromkeys([id_col, *samples]) if c in names]
    frame = text_frame(data, separator, wanted)
    return Text(frame, (), len([c for c in names if c != id_col]))


def maxquant(data: bytes, id_col: str, quantity: str) -> Text:
    """``proteinGroups.txt``: flagged rows dropped and counted, sample columns renamed.

    A ``<quantity> <sample>`` column becomes ``<sample>``.

    Raises:
        ProteoError: If the quantity is unknown or no column carries it.
    """
    if quantity not in MAXQUANT_QUANTITIES:
        raise ProteoError(f"proteo.quantity {quantity!r}: one of {MAXQUANT_QUANTITIES}")
    names = header(data)
    prefix = f"{quantity} "
    # "Intensity " also starts "Intensity L "; a SILAC table is not read, so they stay apart.
    picked = {c: c.removeprefix(prefix) for c in names if c.startswith(prefix)}
    if not picked:
        raise ProteoError(f"proteinGroups.txt: no column starts {prefix!r}")
    flags = [
        (reason, next((c for c in options if c in names), None))
        for reason, options in MAXQUANT_FLAGS
    ]
    present = [c for _, c in flags if c]
    wanted = [c for c in names if c == id_col or c in picked or c in present]
    frame = text_frame(data, "\t", wanted, quoted=False)
    dropped, keep = [], pl.lit(value=True)
    for reason, column in flags:
        hit = keep & (pl.col(column) == FLAGGED) if column else pl.lit(value=False)
        dropped.append((reason, int(frame.select(hit.sum()).item())))
        if column:
            keep &= pl.col(column) != FLAGGED
    frame = frame.filter(keep).drop(present).rename(picked)
    return Text(frame, tuple(dropped), len(picked))


def run_name(column: str) -> str:
    """A DIA-NN run column as a sample name: its file name, without its folder or ending."""
    name = PureWindowsPath(column).name if "\\" in column else PurePosixPath(column).name
    return RUN_ENDINGS.sub("", name)


def diann(data: bytes, id_col: str) -> Text:
    """``pg_matrix.tsv``: each run column renamed to its file's name; nothing to drop.

    Raises:
        ProteoError: If two runs share a file name.
    """
    names = header(data)
    runs = {c: run_name(c) for c in names if c not in DIANN_DESCRIPTORS}
    if len(set(runs.values())) < len(runs):
        raise ProteoError("pg_matrix.tsv: two run columns share a file name")
    wanted = [c for c in names if c == id_col or c in runs]
    frame = text_frame(data, "\t", wanted, quoted=False).rename(runs)
    return Text(frame, (), len(runs))
