"""GTEx V11 from its release files: per-tissue TPM summaries and exon shares for a few genes.

V11 (2025-08-22, GENCODE 47) has no API, so the card reads the files ``ctx setup`` registers:
the gene TPM GCT (2.6 GB gzipped), the exon read parquet (21.8 GB) and the sample and subject
annotation tables. Genes meet by unversioned Ensembl number: NPTX2 is version 4 in GENCODE
47 and 6 in our GENCODE 50. Neither large file is ever held whole: the GCT is streamed line
by line, and the parquet, one row group of 432,367 exons by 19,788 samples, is read in
batches of sample columns, keeping only the wanted genes' rows of each.
"""

import gzip
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl

from ctx.genecard import GenecardError
from ctx.genecard.summary import DETECTED, spread

RELEASE = "V11"
RELEASE_DATE = "2025-08-22"
GENCODE = "47"
# The resource names a dataset file registers the four files under (ctx setup checks them).
RESOURCES = {"tpm": "gtex_gene_tpm", "exons": "gtex_exon_reads",
             "samples": "gtex_samples", "subjects": "gtex_subjects"}  # fmt: skip
# SubjectPhenotypesDD: SEX 1 male, 2 female.
SEXES = {1: "male", 2: "female"}
ALL = "all"
GROUPS = (ALL, "female", "male")
# Sample columns read from the exon parquet at once: 256 float32 columns of 432,367 rows are
# 442,743,808 B, measured on V11 with polars 1.44.2; 128 halves that.
BATCH = 128
_GCT_ID = re.compile(r"ENSG([0-9]{11})\.[0-9]+")
_EXON = re.compile(r"ENSG([0-9]{11})\.[0-9]+_([0-9]+)")


@dataclass(frozen=True)
class Files:
    """The four V11 files, as the dataset file registers them under the root."""

    tpm: Path
    exons: Path
    samples: Path
    subjects: Path


@dataclass(frozen=True)
class GeneTpm:
    """One gene's GCT row: its GENCODE 47 ID and TPM in every sample, in ``samples`` order."""

    gene_id: str
    tpm: np.ndarray


@dataclass(frozen=True)
class ExonReads:
    """One gene's exon read counts: exons of the collapsed model by samples."""

    gene_id: str
    exons: tuple[str, ...]
    reads: np.ndarray


def donors(files: Files) -> pl.DataFrame:
    """Every RNA-seq sample's tissue and donor sex, from the two annotation tables.

    Raises:
        GenecardError: If a sample's donor is missing from the subject table, or its sex is
            coded other than 1 or 2.
    """
    samples = (
        pl.read_csv(
            files.samples,
            separator="\t",
            columns=["SAMPID", "SMTSD", "SMAFRZE"],
            infer_schema=False,
        )  # fmt: skip
        .filter(pl.col("SMAFRZE") == "RNASEQ")
        .select(
            sample="SAMPID", tissue="SMTSD", subject=pl.col("SAMPID").str.extract(r"^([^-]+-[^-]+)")
        )  # fmt: skip
    )
    subjects = pl.read_csv(files.subjects, separator="\t", columns=["SUBJID", "SEX"],
                           schema_overrides={"SEX": pl.Int64})  # fmt: skip
    joined = samples.join(subjects, left_on="subject", right_on="SUBJID", how="left")
    if (lost := joined.filter(pl.col("SEX").is_null())).height:
        raise GenecardError(f"{files.subjects.name}: no donor for {lost.height} samples, "
                            f"such as {lost['sample'][0]}")  # fmt: skip
    if bad := sorted(set(joined["SEX"].to_list()) - set(SEXES)):
        raise GenecardError(f"{files.subjects.name}: SEX codes {bad}; the dictionary has 1, 2")
    return joined.select("sample", "tissue", sex=pl.col("SEX").replace_strict(SEXES))


def read_tpm(path: Path, numbers: Iterable[int]) -> tuple[tuple[str, ...], dict[int, GeneTpm]]:
    """The wanted genes' rows of the GCT, streamed; the file's own dimensions are checked.

    A ``_PAR_Y`` row shares its number with the X copy, so it is never a wanted gene's row.

    Returns:
        The sample IDs in column order, and each found gene's row by unversioned number.

    Raises:
        GenecardError: If the file is not GCT 1.2, its rows or columns differ from the
            dimensions it declares, or one number names two rows.
    """
    wanted, found = set(numbers), {}
    with gzip.open(path, "rt") as handle:
        if handle.readline().strip() != "#1.2":
            raise GenecardError(f"{path.name}: not a GCT 1.2 file")
        rows, columns = (int(v) for v in handle.readline().split("\t")[:2])
        samples = tuple(handle.readline().rstrip("\n").split("\t")[2:])
        seen = 0
        for line in handle:
            seen += 1
            gene_id = line[: line.index("\t")]
            match = _GCT_ID.fullmatch(gene_id)
            if match is None or (n := int(match[1])) not in wanted:
                continue
            if n in found:
                raise GenecardError(f"{path.name}: ENSG{n:011d} names two rows")
            values = np.array(line.rstrip("\n").split("\t")[2:], dtype=float)
            found[n] = GeneTpm(gene_id, values)
    if (seen, len(samples)) != (rows, columns):
        raise GenecardError(f"{path.name}: {seen} rows by {len(samples)} samples; the header "
                            f"declares {rows} by {columns}")  # fmt: skip
    return samples, found


def tissue_spread(tpm: np.ndarray, samples: tuple[str, ...], who: pl.DataFrame) -> pl.DataFrame:
    """n, median and IQR of TPM per tissue, for all donors, then female, then male.

    A sex with no sample in a tissue (prostate, uterus) has no row: absent, never zero.

    Raises:
        GenecardError: If a sample of the GCT has no annotation row.
    """
    order = pl.DataFrame({"sample": samples, "tpm": tpm})
    joined = order.join(who, on="sample", how="left", maintain_order="left")
    if (lost := joined.filter(pl.col("tissue").is_null())).height:
        raise GenecardError(f"{lost.height} GCT samples have no annotation, such as "
                            f"{lost['sample'][0]}")  # fmt: skip
    rows = []
    for (tissue,), part in sorted(joined.group_by("tissue"), key=lambda kv: str(kv[0][0])):
        for group in GROUPS:
            values = part["tpm"] if group == ALL else part.filter(pl.col("sex") == group)["tpm"]
            if values.len():
                s = spread(values.to_numpy())
                rows.append((tissue, group, s.n, s.median, s.q1, s.q3))
    schema = {"tissue": pl.String, "sex": pl.String, "n": pl.Int64, "median": pl.Float64,
              "q1": pl.Float64, "q3": pl.Float64}  # fmt: skip
    return pl.DataFrame(rows, schema=schema, orient="row")


def detected(spread: pl.DataFrame, floor: float = DETECTED) -> list[str]:
    """Tissues whose all-donor median TPM reaches ``floor``, in the spread's order.

    The exon panel reads exon use only there (H.3): a tissue that barely reads the gene gives
    shares from a handful of reads, mostly 0, and sets every reference to 0.
    """
    rows = spread.filter((pl.col("sex") == ALL) & (pl.col("median") >= floor))
    return rows["tissue"].to_list()


def read_exons(path: Path, numbers: Iterable[int], samples: tuple[str, ...],
               batch: int = BATCH) -> dict[int, ExonReads]:  # fmt: skip
    """The wanted genes' exon reads, in ``samples`` order, read a batch of columns at a time.

    Raises:
        GenecardError: If the parquet lacks a sample of ``samples``.
    """
    wanted = set(numbers)
    scan = pl.scan_parquet(path)
    if missing := set(samples) - set(scan.collect_schema().names()):
        raise GenecardError(f"{path.name}: no column for {len(missing)} GCT samples")
    names = scan.select("Name").collect()["Name"].to_list()
    rows: dict[int, list[tuple[int, int, str]]] = {}
    for row, name in enumerate(names):
        if (match := _EXON.fullmatch(name)) and (n := int(match[1])) in wanted:
            rows.setdefault(n, []).append((int(match[2]), row, name))
    order = {n: sorted(found) for n, found in rows.items()}
    picked = [row for found in order.values() for _, row, _ in found]
    reads = np.empty((len(picked), len(samples)), dtype=np.float64)
    for start in range(0, len(samples), batch):
        columns = list(samples[start : start + batch])
        block = scan.select(columns).collect()[picked]
        reads[:, start : start + len(columns)] = block.to_numpy()
    out, at = {}, 0
    for n, found in order.items():
        gene_id = found[0][2].rsplit("_", 1)[0]
        out[n] = ExonReads(gene_id, tuple(f"E{k}" for k, _, _ in found),
                           reads[at : at + len(found)])  # fmt: skip
        at += len(found)
    return out


def exon_shares(exons: ExonReads, samples: tuple[str, ...], who: pl.DataFrame
                ) -> tuple[pl.DataFrame, int]:  # fmt: skip
    """Each exon's share of the gene's exon reads: the median over a tissue's samples.

    The share is taken in each sample first, then the median across samples, so a few deep
    samples cannot set a tissue's profile. A sample with no exon read has no share.

    Returns:
        One row per exon and tissue (exon, tissue, n, share), and the samples left out for
        having no read on the gene.
    """
    totals = exons.reads.sum(axis=0)
    covered = totals > 0
    share = exons.reads[:, covered] / totals[covered]
    tissue_of: Mapping[str, str] = dict(who.select("sample", "tissue").iter_rows())
    tissues = np.array([tissue_of[s] for s, keep in zip(samples, covered, strict=True) if keep])
    rows = []
    for tissue in sorted(set(tissues.tolist())):
        here = tissues == tissue
        medians = np.median(share[:, here], axis=1)
        rows += [(e, tissue, int(here.sum()), float(m))
                 for e, m in zip(exons.exons, medians, strict=True)]  # fmt: skip
    schema = {"exon": pl.String, "tissue": pl.String, "n": pl.Int64, "share": pl.Float64}
    return pl.DataFrame(rows, schema=schema, orient="row"), int((~covered).sum())
