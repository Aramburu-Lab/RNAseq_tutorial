"""The one gene summary of Salmon's transcripts: R tximport's summarizeToGene, in its order.

tximport 1.30.0 (``summarizeToGene.R``, the release in nf-core/rnaseq 3.27.0's image), step by
step: per sample, a gene's abundance is the sum of its transcripts' TPM, its counts the sum of
their NumReads, and its length the TPM-weighted mean of their EffectiveLength. A gene with no
TPM in a sample has no weighted length there: it takes the geometric mean of its lengths in
the samples that have one, or, with none anywhere, the plain mean over its transcripts of
their lengths averaged over samples. Only then are counts made from abundance:
``length_scaled_tpm`` multiplies TPM by the gene's mean length over samples, ``scaled_tpm``
takes TPM as it is, and either is scaled so each sample keeps the sum of its estimated counts.
Taking the mean length before the gaps are filled, or the arithmetic for the geometric mean,
gives other counts (``tests/test_quant_tximport.py`` plants both).

Samples stream in one at a time, so 210 samples of 644,914 transcripts never sit in memory
together: only genes by samples does. A transcript the transcript-to-gene table lacks is its
own gene while counts are scaled and is dropped after, as nf-core's ``tximport.r`` does: its
reads stay in each sample's total, so leaving it out first gives other scaled counts. Genes
come out in byte order, as R's ``rowsum(reorder = TRUE)`` under the C collation nf-core sets.
"""

from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np
import polars as pl
from numpy.typing import NDArray

from ctx.quant import CountsFrom, QuantError

# quant.sf's columns, by name (Salmon's documentation, "Quantification output").
NAME, LENGTH, EFFECTIVE, TPM, READS = "Name", "Length", "EffectiveLength", "TPM", "NumReads"
QUANT_COLUMNS = (NAME, LENGTH, EFFECTIVE, TPM, READS)

Matrix = NDArray[np.float64]


@dataclass(frozen=True)
class GeneTables:
    """tximport's gene-level tables, genes by samples.

    Every matrix holds one row per entry of ``rows``: the genes, and each transcript without
    a gene standing as its own; ``genes`` and the properties keep the genes alone.

    Attributes:
        rows: Gene identifiers and gene-less transcripts, in byte order.
        samples: Sample names, in the order given.
        summed: Summed estimated counts (countsFromAbundance "no").
        abundance: Summed TPM.
        length: TPM-weighted effective length, gaps filled.
        real: Whether each row is a gene of the transcript-to-gene table.
    """

    rows: tuple[str, ...]
    samples: tuple[str, ...]
    summed: Matrix
    abundance: Matrix
    length: Matrix
    real: NDArray[np.bool_]

    @property
    def genes(self) -> tuple[str, ...]:
        """The genes, in byte order."""
        return tuple(g for g, keep in zip(self.rows, self.real, strict=True) if keep)

    @property
    def counts(self) -> Matrix:
        """Summed estimated counts of the genes."""
        return self.summed[self.real]

    @property
    def tpm(self) -> Matrix:
        """Summed TPM of the genes."""
        return self.abundance[self.real]

    @property
    def lengths(self) -> Matrix:
        """Effective lengths of the genes."""
        return self.length[self.real]

    def counts_from(self, mode: CountsFrom) -> Matrix:
        """Gene counts as countsFromAbundance makes them; each sample keeps its counts' sum."""
        if mode is CountsFrom.NO:
            return self.counts
        made = (self.abundance if mode is CountsFrom.SCALED_TPM
                else self.abundance * self.length.mean(axis=1, keepdims=True))  # fmt: skip
        totals = made.sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            scale = np.where(totals > 0, self.summed.sum(axis=0) / totals, 0.0)
        scaled: Matrix = (made * scale)[self.real]
        return scaled

    def frame(self, values: Matrix, names: dict[str, str] | None = None) -> pl.DataFrame:
        """One table as gene_id (and gene_name, where ``names`` is given), then each sample."""
        frame = pl.DataFrame(values, schema=list(self.samples), orient="row")
        if names is not None:
            named = [names.get(g, g) for g in self.genes]
            frame = frame.insert_column(0, pl.Series("gene_name", named))
        return frame.insert_column(0, pl.Series("gene_id", self.genes))


@dataclass
class _Sums:
    """What one pass over the samples collects before any gap is filled."""

    tpm: list[Matrix]
    counts: list[Matrix]
    weighted: list[Matrix]
    length_total: Matrix


def _gene_index(names: pl.Series, tx2gene: pl.DataFrame
                ) -> tuple[NDArray[np.int64], list[str], set[str]]:  # fmt: skip
    """Each transcript's row number, the rows in byte order, and the transcripts with no gene."""
    mapped = dict(zip(tx2gene["transcript_id"], tx2gene["gene_id"], strict=True))
    extra = {n for n in names if n not in mapped}
    gene_of = [mapped.get(n, n) for n in names]
    rows = sorted(set(gene_of), key=lambda g: g.encode())
    number = {g: i for i, g in enumerate(rows)}
    return np.array([number[g] for g in gene_of], dtype=np.int64), rows, extra


def _rowsum(index: NDArray[np.int64], values: Matrix, n_rows: int) -> Matrix:
    """Sum transcript values into their rows, in transcript order as R's rowsum adds them."""
    summed: Matrix = np.bincount(index, weights=values, minlength=n_rows)
    return summed


def _fill(lengths: Matrix, mean_length: Matrix) -> Matrix:
    """Fill a gene's missing lengths as tximport's replaceMissingLength does.

    Args:
        lengths: Genes by samples, NaN where a gene had no TPM.
        mean_length: Per gene, the plain mean over its transcripts of their mean length.
    """
    filled = lengths.copy()
    for i in np.flatnonzero(np.isnan(filled).any(axis=1)):
        missing = np.isnan(filled[i])
        if missing.all():
            filled[i] = mean_length[i]
        else:
            filled[i, missing] = np.exp(np.mean(np.log(filled[i, ~missing])))
    return filled


def summarize(quants: Iterable[tuple[str, pl.DataFrame]], tx2gene: pl.DataFrame) -> GeneTables:
    """Summarise every sample's transcripts to genes, as tximport's summarizeToGene does.

    Args:
        quants: (sample, quant.sf as read by ``ctx.quant.outdir.quant``) per sample; every
            table must list the same transcripts in the same order, as one index writes them.
        tx2gene: transcript_id and gene_id.

    Returns:
        The gene tables; counts from abundance are made from them by ``counts_from``.

    Raises:
        QuantError: If no sample is given, or two samples list different transcripts.
    """
    samples: list[str] = []
    names: pl.Series | None = None
    sums = _Sums([], [], [], np.zeros(0))
    index, rows, extra = np.zeros(0, dtype=np.int64), list[str](), set[str]()
    for sample, quant in quants:
        if names is None:
            names = quant[NAME]
            index, rows, extra = _gene_index(names, tx2gene)
            sums.length_total = np.zeros(len(names))
        elif not quant[NAME].equals(names):
            raise QuantError(f"{sample}: its quant.sf lists other transcripts than {samples[0]}'s")
        samples.append(sample)
        tpm, effective = quant[TPM].to_numpy(), quant[EFFECTIVE].to_numpy()
        sums.tpm.append(_rowsum(index, tpm, len(rows)))
        sums.counts.append(_rowsum(index, quant[READS].to_numpy(), len(rows)))
        sums.weighted.append(_rowsum(index, tpm * effective, len(rows)))
        sums.length_total = sums.length_total + effective
    if not samples:
        raise QuantError("no sample to summarise")
    tpm, weighted = np.column_stack(sums.tpm), np.column_stack(sums.weighted)
    with np.errstate(invalid="ignore", divide="ignore"):
        lengths = weighted / tpm
    per_transcript = sums.length_total / len(samples)
    mean_length = _rowsum(index, per_transcript, len(rows)) / np.bincount(
        index, minlength=len(rows)
    )
    return GeneTables(
        tuple(rows),
        tuple(samples),
        np.column_stack(sums.counts),
        tpm,
        _fill(lengths, mean_length),
        np.array([r not in extra for r in rows]),
    )
