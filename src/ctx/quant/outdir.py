"""Read Salmon's outputs as nf-core/rnaseq 3.27.0 writes them: by column name, never position.

Two folders. The rnaseq outdir ``ctx gather`` reads holds the length-scaled gene matrix and
each sample's ``aux_info/meta_info.json`` and ``lib_format_counts.json``; those two are logs,
read only through MultiQC (``ctx.qc.multiqc_io``). The per-sample copy Q1 holds the salmon
folder's layout: each sample's ``quant.sf``, the merged ``salmon.merged.*.tsv`` tables and
``tx2gene.tsv``, written by nf-core's ``tximport.r`` (transcript tables keyed by ``tx`` and
``gene_id``, gene tables by ``gene_id`` and ``gene_name``). Nothing here writes into either.
"""

from collections.abc import Iterator, Sequence
from pathlib import Path

import polars as pl

from ctx.qc import QcError, multiqc_io
from ctx.quant import QuantError
from ctx.quant.tximport import QUANT_COLUMNS

MERGED = "salmon.merged.{}.tsv"
TX2GENE = "tx2gene.tsv"
TX2GENE_COLUMNS = ("transcript_id", "gene_id", "gene_name")
GENE_KEYS, TRANSCRIPT_KEYS = ("gene_id", "gene_name"), ("tx", "gene_id")
# The merged tables M4 reads from Q1, and their key columns.
GENE_TABLES = ("gene_counts", "gene_tpm", "gene_lengths", "gene_counts_length_scaled",
               "gene_counts_scaled")  # fmt: skip
TRANSCRIPT_TABLES = ("transcript_tpm",)
# MultiQC's Salmon table: the fates, the library and the run's settings, by column.
FATES = ("num_processed", "num_mapped", "num_decoy_fragments", "num_fragments_filtered_vm")
SETTINGS = ("library_types", "frag_length_mean", "frag_length_sd", "frag_length_source",
            "gc_bias_correct", "seq_bias_correct", "opt_type", "salmon_version",
            "num_valid_targets", "num_decoy_targets")  # fmt: skip
STRAND_BIAS = "salmon-strand_mapping_bias"
SHOWN = 5


def _shown(names: Sequence[str]) -> str:
    """The first names of a list, and how many more."""
    more = f" and {len(names) - SHOWN} more" if len(names) > SHOWN else ""
    return f"{list(names[:SHOWN])}{more}"


def header(path: Path) -> list[str]:
    """A tab-separated table's column names."""
    return pl.scan_csv(path, separator="\t", n_rows=0).collect_schema().names()


def table_problems(path: Path, keys: Sequence[str], samples: Sequence[str]) -> list[str]:
    """Why a merged table cannot serve: missing, its keys absent, a sample missing."""
    if not path.is_file():
        return [f"{path} is missing"]
    columns = header(path)
    problems = []
    if missing := [k for k in keys if k not in columns]:
        problems.append(f"{path.name} has no {missing} column")
    if absent := [s for s in samples if s not in columns]:
        problems.append(f"{path.name} has no column for {_shown(absent)}")
    return problems


def merged(folder: Path, what: str, keys: Sequence[str], samples: Sequence[str]) -> pl.DataFrame:
    """One merged table, its keys as text and the samples' columns as Float64, in that order."""
    path = folder / MERGED.format(what)
    schema = {**dict.fromkeys(keys, pl.String), **dict.fromkeys(samples, pl.Float64)}
    return pl.scan_csv(path, separator="\t", schema_overrides=schema).select(*keys,
                                                                          *samples).collect()


def column_sums(folder: Path, what: str, samples: Sequence[str]) -> dict[str, float]:
    """Each sample's column sum of a merged table, streamed: the table is never held whole."""
    path = folder / MERGED.format(what)
    sums = (pl.scan_csv(path, separator="\t", schema_overrides=dict.fromkeys(samples, pl.Float64))
            .select(pl.col(list(samples)).sum()).collect())  # fmt: skip
    return {s: float(sums[s][0]) for s in samples}


def tx2gene(folder: Path) -> pl.DataFrame:
    """nf-core's transcript-to-gene table: transcript_id, gene_id, gene_name.

    Raises:
        QuantError: If the table is missing or lacks one of those columns.
    """
    path = folder / TX2GENE
    if not path.is_file():
        raise QuantError(f"{path} is missing")
    if missing := [c for c in TX2GENE_COLUMNS if c not in header(path)]:
        raise QuantError(f"{path.name} has no {missing} column")
    return pl.read_csv(path, separator="\t", infer_schema=False).select(TX2GENE_COLUMNS)


def quant_problems(folder: Path, samples: Sequence[str]) -> list[str]:
    """Samples with no quant.sf, and a quant.sf without Salmon's columns, checked on the first."""
    absent = [s for s in samples if not (folder / s / "quant.sf").is_file()]
    problems = [f"{folder}: no quant.sf for {_shown(absent)}"] if absent else []
    present = [s for s in samples if s not in absent]
    if present and (missing := [c for c in QUANT_COLUMNS
                                if c not in header(folder / present[0] / "quant.sf")]):  # fmt: skip
        problems.append(f"{present[0]}/quant.sf has no {missing} column")
    return problems


def quant(path: Path) -> pl.DataFrame:
    """One quant.sf, Salmon's five columns by name, the numbers as Float64."""
    numbers = dict.fromkeys(QUANT_COLUMNS[1:], pl.Float64)
    return pl.read_csv(path, separator="\t", schema_overrides=numbers).select(QUANT_COLUMNS)


def quants(folder: Path, samples: Sequence[str]) -> Iterator[tuple[str, pl.DataFrame]]:
    """Every sample's quant.sf, one at a time, for ``tximport.summarize``."""
    for sample in samples:
        yield sample, quant(folder / sample / "quant.sf")


def salmon_run(rnaseq_dir: Path, dest: Path, samples: Sequence[str]) -> multiqc_io.MultiqcRun:
    """MultiQC's Salmon module on the outdir's salmon folder; every sample must be there.

    Raises:
        QuantError: If MultiQC fails, or a sample has no row or one it does not know.
    """
    try:
        run = multiqc_io.run([rnaseq_dir / "salmon"], dest, modules=("salmon",))
        names = multiqc_io.map_names(multiqc_io.table(run, "salmon")["Sample"],
                                     {s: s for s in samples})  # fmt: skip
    except QcError as error:
        raise QuantError(str(error)) from error
    if missing := sorted(set(samples) - set(names.values())):
        raise QuantError(f"MultiQC found no Salmon run for {_shown(missing)}")
    return run


def salmon_table(run: multiqc_io.MultiqcRun) -> pl.DataFrame:
    """One row per sample: Salmon's fates as integers, its settings as text, its strand bias.

    Raises:
        QuantError: If MultiQC's table lacks a column read here: a zero from a column the
            table does not carry is not a result.
    """
    table = multiqc_io.table(run, "salmon")
    stats = multiqc_io.table(run, "general_stats")
    for name, frame, needed in (("salmon", table, (*FATES, *SETTINGS)),
                                ("general_stats", stats, (STRAND_BIAS,))):  # fmt: skip
        if missing := [c for c in needed if c not in frame.columns]:
            raise QuantError(f"MultiQC's {name} table has no {missing}")
    fates = multiqc_io.numeric(table.select("Sample", *FATES), FATES)
    return (fates.with_columns(pl.col(list(FATES)).cast(pl.Int64))
            .join(table.select("Sample", *SETTINGS), on="Sample")
            .join(multiqc_io.numeric(stats.select("Sample", STRAND_BIAS), (STRAND_BIAS,)),
                  on="Sample", how="left")
            .rename({"Sample": "sample", STRAND_BIAS: "strand_bias"}))  # fmt: skip
