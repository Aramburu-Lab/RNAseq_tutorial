"""The miRNA half of a GDC project: its miRNA and isoform files, read file by file.

``ctx gather --mirna-dir --isoform-dir`` beside ``--gdc-dir`` is the only command that opens
these folders. A miRNA file has one row per miRBase hairpin; only an isoform file names the
mature miRNA a read came from (``miRNA_region`` ``mature,MIMAT...``), so the mature matrix
sums the mature rows of every hairpin of a MIMAT, each read once: a read that maps to
let-7a-1, -2 and -3 alike is listed under one of them. Each isoform file is aggregated as it
is read, so no file's rows are held beside another's. The two kinds must list the same
aliquots, and per aliquot every hairpin's isoform rows must sum to the miRNA file's count.
Files are chosen by ``gdc.choose`` without its protocol step, and selected through
``gdc.chosen``, as the STAR files are; each is joined to the chosen STAR file of its sample
barcode, never by aliquot, which differs between the two libraries of one sample.
"""

import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl

from ctx.dataset import GdcSpec
from ctx.gather import gdc

MIRNA, ISOFORM = "mirna", "isoform"
# Per kind: the fetch's metadata and manifest files and the GDC data type its query names.
KINDS = {
    MIRNA: ("gdc_files_mirna.tsv", "gdc_manifest_mirna.txt", "miRNA Expression Quantification"),
    ISOFORM: ("gdc_files_isoform.tsv", "gdc_manifest_isoform.txt",
              "Isoform Expression Quantification"),
}  # fmt: skip
MIRNA_COLUMNS = ("miRNA_ID", "read_count", "reads_per_million_miRNA_mapped", "cross-mapped")
ISOFORM_COLUMNS = ("miRNA_ID", "isoform_coords", "read_count", "reads_per_million_miRNA_mapped",
                   "cross-mapped", "miRNA_region")  # fmt: skip
# The region classes BCGSC writes; a mature row carries its MIMAT after a comma.
MATURE = "mature"
REGIONS = (MATURE, "precursor", "stemloop", "unannotated")
COORDS = re.compile(r"^[^:]+:(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+):(?P<strand>[+-])$")


@dataclass(frozen=True)
class Matrices:
    """Both kinds of one project, aligned on the aliquots of ``files``.

    Attributes:
        files: One row per aliquot: sample (the aliquot), the GDC fields of ``gdc.FIELDS`` from
            the miRNA kind, and ``isoform_file_id``.
        hairpins: Hairpin IDs, in the files' order.
        hairpin_counts: int32, hairpins by aliquots, from the miRNA files.
        mature: MIMAT accessions, sorted.
        mature_counts: int32, MIMATs by aliquots: the mature rows' reads, summed over hairpins.
        loci: Per hairpin: chrom, start, end and strand spanning every isoform row of it in any
            file (1-based, inclusive), and the MIMATs its mature rows name; null where no file
            holds a row of it.
        stats: Per aliquot: reads of each region class, cross-mapped reads of the isoform rows,
            and the distinct mature isoform rows.
        paths: The miRNA folder, then the isoform folder.
    """

    files: pl.DataFrame
    hairpins: tuple[str, ...]
    hairpin_counts: np.ndarray
    mature: tuple[str, ...]
    mature_counts: np.ndarray
    loci: pl.DataFrame
    stats: pl.DataFrame
    paths: tuple[Path, Path]


def read_kind(path: Path, spec: GdcSpec, kind: str) -> gdc.Folder:
    """One kind's folder: metadata, manifest, release, retrieval log and query checked.

    Raises:
        gdc.GdcError: Listing every problem found.
    """
    metadata, manifest, data_type = KINDS[kind]
    files, dropped = gdc._metadata(path, spec, metadata)  # noqa: SLF001 - the one metadata reader
    if files.is_empty():
        raise gdc.GdcError(f"{path / metadata}: no {kind} file of this project and site")
    problems = gdc._manifest(path, files, manifest)  # noqa: SLF001 - the one manifest check
    found = json.loads(gdc._text(path / gdc.STATUS)).get("data_release")  # noqa: SLF001
    if found != spec.release:
        problems.append(f"{gdc.STATUS}: {found!r}, not {spec.release!r}")
    log = gdc._table(path / gdc.RETRIEVAL)  # noqa: SLF001 - the one TSV reader
    problems += [f"{gdc.RETRIEVAL}: {r['file']} retrieved from {r['gdc_release']!r}"
                 for r in log if r["gdc_release"] != spec.release][:5]  # fmt: skip
    logged = {r["file"] for r in log}
    needed = [metadata, manifest, gdc.STATUS]
    needed += [f"{i}/{n}" for i, n in files.select("file_id", "file_name").rows()]
    problems += [f"{gdc.RETRIEVAL}: no row for {n}" for n in needed if n not in logged][:5]
    wanted = (("project_id", spec.project), ("workflow_type", spec.mirna_workflow),
              ("data_type", data_type))  # fmt: skip
    for field, value in wanted:
        if gdc._queried(path, metadata, field) != {value}:  # noqa: SLF001 - the one query reader
            problems.append(f"{metadata}{gdc.QUERY}: {field} is not {value!r} alone")
    if problems:
        raise gdc.GdcError(f"{path}:\n  " + "\n  ".join(problems))
    return gdc.Folder(path, files, spec.release, dropped)


def _columns(path: Path, header: list[str], wanted: tuple[str, ...]) -> None:
    """A file holds exactly the expected columns, by name, in any order."""
    if sorted(header) != sorted(wanted):
        raise gdc.GdcError(f"{path.name}: columns {header}, not {list(wanted)}")


def read_hairpins(folder: gdc.Folder) -> tuple[tuple[str, ...], np.ndarray]:
    """Every miRNA file's read counts, hairpins by files, each file's hairpins those of the first.

    Raises:
        gdc.GdcError: On the first file whose columns, hairpins or counts do not fit.
    """
    hairpins: tuple[str, ...] = ()
    matrix = np.empty(0, dtype=np.int32)
    for j, (file_id, name) in enumerate(folder.files.select("file_id", "file_name").rows()):
        path = folder.path / file_id / name
        frame = pl.read_csv(path, separator="\t", schema_overrides={"read_count": pl.Int64})
        _columns(path, frame.columns, MIRNA_COLUMNS)
        if not hairpins:
            hairpins = tuple(frame["miRNA_ID"].to_list())
            matrix = np.zeros((len(hairpins), folder.files.height), dtype=np.int32)
        if tuple(frame["miRNA_ID"].to_list()) != hairpins:
            raise gdc.GdcError(f"{name}: hairpins differ from {folder.files['file_name'][0]}")
        counts = frame["read_count"].to_numpy()
        if counts.min() < 0 or counts.max() > gdc.INT32_MAX:
            raise gdc.GdcError(f"{name}: a count outside int32")
        matrix[:, j] = counts
    return hairpins, matrix


def _one_isoform(path: Path) -> pl.DataFrame:
    """One isoform file's rows, its columns and region classes checked."""
    frame = pl.read_csv(path, separator="\t", schema_overrides={"read_count": pl.Int64})
    _columns(path, frame.columns, ISOFORM_COLUMNS)
    region = frame["miRNA_region"].str.split(",").list.first()
    if stray := sorted(set(region.to_list()) - set(REGIONS)):
        raise gdc.GdcError(f"{path.name}: miRNA_region {stray}, not one of {REGIONS}")
    coords = frame["isoform_coords"].str.extract_groups(COORDS.pattern)
    if coords.struct.field("chrom").null_count():
        raise gdc.GdcError(f"{path.name}: an isoform_coords value is not genome:chrom:a-b:strand")
    return frame.with_columns(region=region, mimat=frame["miRNA_region"].str.extract(r",(.+)$"),
                              coords=coords)  # fmt: skip


def _loci(parts: list[pl.DataFrame], hairpins: tuple[str, ...]) -> pl.DataFrame:
    """Each hairpin's span over every file, and the MIMATs its mature rows name."""
    rows = pl.concat(parts).group_by("miRNA_ID").agg(
        pl.col("chrom").unique().sort(), pl.col("strand").unique().sort(),
        pl.col("start").min(), pl.col("end").max(),
        pl.col("mimat").drop_nulls().unique().sort().str.join(";").alias("mimats"),
    )  # fmt: skip
    bad = rows.filter((pl.col("chrom").list.len() > 1) | (pl.col("strand").list.len() > 1))
    if bad.height:
        raise gdc.GdcError(f"hairpins on two chromosomes or strands: {bad['miRNA_ID'][:5]}")
    rows = rows.with_columns(pl.col("chrom").list.first(), pl.col("strand").list.first())
    order = pl.DataFrame({"miRNA_ID": list(hairpins)})
    return order.join(rows, on="miRNA_ID", how="left", maintain_order="left").select(
        "miRNA_ID", "chrom", "start", "end", "strand", "mimats"
    )


def read_isoforms(folder: gdc.Folder, hairpins: tuple[str, ...]) -> tuple[
        tuple[str, ...], np.ndarray, np.ndarray, pl.DataFrame, pl.DataFrame]:  # fmt: skip
    """Every isoform file aggregated as read: mature and hairpin sums, regions, loci.

    Returns:
        The MIMATs; mature counts, MIMATs by files; each file's isoform reads per hairpin,
        hairpins by files, for the check against the miRNA files; per file: reads per region,
        cross-mapped reads, mature isoform rows; and the hairpins' loci.

    Raises:
        gdc.GdcError: On a file whose columns or regions do not fit, or a hairpin unknown to
            the miRNA files.
    """
    index = {h: i for i, h in enumerate(hairpins)}
    per_hairpin = np.zeros((len(hairpins), folder.files.height), dtype=np.int64)
    mature: list[dict[str, int]] = []
    stats, spans = [], []
    for j, (file_id, name) in enumerate(folder.files.select("file_id", "file_name").rows()):
        frame = _one_isoform(folder.path / file_id / name)
        if stray := sorted(set(frame["miRNA_ID"]) - set(index)):
            raise gdc.GdcError(f"{name}: hairpins {stray[:5]} are in no miRNA file")
        by = frame.group_by("miRNA_ID").agg(pl.col("read_count").sum())
        per_hairpin[[index[h] for h in by["miRNA_ID"]], j] = by["read_count"].to_numpy()
        rows = frame.filter(pl.col("region") == MATURE)
        mature.append(dict(rows.group_by("mimat").agg(pl.col("read_count").sum()).rows()))
        regions = dict(frame.group_by("region").agg(pl.col("read_count").sum()).rows())
        stats.append({**{f"{r}_reads": int(regions.get(r, 0)) for r in REGIONS},
                      "cross_mapped_reads": int(frame.filter(pl.col("cross-mapped") == "Y")[
                          "read_count"].sum()),
                      "mature_isoforms": rows.height})  # fmt: skip
        span = frame.select("miRNA_ID", "mimat", pl.col("coords").struct.unnest())
        spans.append(span.with_columns(pl.col("start", "end").cast(pl.Int64)))
    names = tuple(sorted({m for found in mature for m in found}))
    matrix = np.zeros((len(names), len(mature)), dtype=np.int64)
    at = {m: i for i, m in enumerate(names)}
    for j, found in enumerate(mature):
        matrix[[at[m] for m in found], j] = list(found.values())
    if matrix.max(initial=0) > gdc.INT32_MAX:
        raise gdc.GdcError("a mature count outside int32")
    return (names, matrix.astype(np.int32), per_hairpin, pl.DataFrame(stats),
            _loci(spans, hairpins))  # fmt: skip


def _aligned(mirna: gdc.Folder, isoform: gdc.Folder) -> pl.DataFrame:
    """The miRNA kind's files with the isoform file of each aliquot; both list the same ones."""
    left, right = set(mirna.files["sample"]), set(isoform.files["sample"])
    if left != right:
        raise gdc.GdcError(f"miRNA and isoform folders list different aliquots: "
                           f"{sorted(left ^ right)[:5]}")  # fmt: skip
    ids = isoform.files.select("sample", isoform_file_id="file_id",
                               isoform_file_name="file_name")  # fmt: skip
    return mirna.files.join(ids, on="sample", maintain_order="left")


def read_matrices(mirna: gdc.Folder, isoform: gdc.Folder) -> Matrices:
    """Both kinds read, aligned on aliquots and checked against each other per hairpin.

    Raises:
        gdc.GdcError: If the folders list other aliquots, or any aliquot's isoform rows of a
            hairpin do not sum to its miRNA file's count; every such aliquot is named.
    """
    files = _aligned(mirna, isoform)
    hairpins, hairpin_counts = read_hairpins(mirna)
    ordered = gdc.Folder(isoform.path, files.select("sample", file_id="isoform_file_id",
                                                    file_name="isoform_file_name"),
                         isoform.release)  # fmt: skip
    names, mature, per_hairpin, stats, loci = read_isoforms(ordered, hairpins)
    differ = np.flatnonzero((per_hairpin != hairpin_counts).any(axis=0))
    if differ.size:
        shown = [files["sample"][int(j)] for j in differ[:5]]
        raise gdc.GdcError(f"{differ.size} aliquots' isoform rows do not sum to their miRNA "
                           f"files' hairpin counts: {shown}")  # fmt: skip
    stats = stats.with_columns(
        sample=files["sample"], mirna_reads=pl.Series(hairpin_counts.sum(axis=0, dtype=np.int64))
    )
    return Matrices(files.drop("isoform_file_name"), hairpins, hairpin_counts, names, mature,
                    loci, stats, (mirna.path, isoform.path))  # fmt: skip


def sheet(found: Matrices, notes: pl.DataFrame, rna: pl.DataFrame, spec: GdcSpec) -> pl.DataFrame:
    """One row per aliquot: identity, reads, the choice, and the STAR file of its sample.

    Args:
        found: Both kinds, read.
        notes: ``gdc.annotations`` of the STAR folder over these files; a redaction drops one.
        rna: Gather's STAR sample sheet; its chosen file of the same sample barcode is joined.
        spec: The dataset file's ``[gdc]``.
    """
    hit = notes.filter(pl.col("action") == "drops its files")["samples"]
    redacted = {s for joined in hit for s in joined.split(";")}
    condition = (pl.when(pl.col("sample_type") == spec.tumour).then(pl.lit("tumour"))
                 .when(pl.col("sample_type") == spec.normal).then(pl.lit("normal"))
                 .otherwise(pl.lit("other")))  # fmt: skip
    frame = found.files.join(found.stats, on="sample", maintain_order="left").with_columns(
        condition=condition, donor=pl.col("case"),
        run_accession=gdc.run_accession(),
        redacted=pl.col("sample").is_in(sorted(redacted)),
    )  # fmt: skip
    frame = gdc.choose(frame, None)
    taken = gdc.chosen(rna).select("sample_id", rna_sample="sample")
    return frame.join(taken, on="sample_id", how="left", maintain_order="left")
