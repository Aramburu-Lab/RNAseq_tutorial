"""A made-up ctx root for the M9 tests: gather's genes, M5 results, GTEx V11-shaped files.

The GTEx files have V11's shapes (GCT 1.2, the exon parquet's Name and Description beside one
float32 column per sample, the two annotation tables) and made-up values. NPTX2's HPA and
Ensembl responses are the real ones in ``tests/fixtures/genecard``; ACTB has none, so its
HPA and Ensembl panels are unavailable, which the degraded-card tests read.
"""

import gzip
import hashlib
import io
import json
import urllib.error
from datetime import date
from pathlib import Path

import numpy as np
import polars as pl

from ctx.genecard import annotation, hpa, junctions
from ctx.io.naming import stem
from ctx.io.readme_companion import write_table

FIXTURES = Path(__file__).parent / "fixtures" / "genecard"
SLUG = "gcsynth"
NPTX2, ACTB = 106236, 75624
# Our IDs (GENCODE 50) and GTEx's (GENCODE 47): the versions differ, as they do for real.
OURS = {NPTX2: "ENSG00000106236.6", ACTB: "ENSG00000075624.19"}
THEIRS = {NPTX2: "ENSG00000106236.4", ACTB: "ENSG00000075624.17"}
# Two IDs under one symbol, as GENCODE gives 446 of Step-1's symbols; and three called genes.
TWINS = ("ENSG00000900001.1", "ENSG00000900002.1")
# GENCODE 50's name for HPA v25's NDUFA4L2, whose search lists COXFA4L2 as a synonym.
RENAMED = ("ENSG00000185633.13", "COXFA4L2")
# The dataset file's HPA cell-line group and lines, as Non-cancerous-blood.toml names them.
LINES = ("cell_RNA_leukemia", ("HL-60", "THP-1", "U-937"))
# padj and |log2 FC| disagree on UPA and UPB, and UPA and UPC tie on padj.
CALLED = {"ENSG00000900011.1": ("UPA", 1e-9, 3.0, "up"),
          "ENSG00000900012.1": ("UPB", 1e-8, 4.0, "up"),
          "ENSG00000900014.1": ("UPC", 1e-9, 5.0, "up"),
          "ENSG00000900013.1": ("DOWNA", 1e-3, -2.0, "down")}  # fmt: skip
# Donor, sex code (1 male, 2 female); tissue and the donors sampled in it.
SUBJECTS = (("GTEX-AAAA", 2), ("GTEX-BBBB", 1), ("GTEX-CCCC", 1), ("GTEX-DDDD", 2))
TISSUES = {"Pituitary": ("GTEX-AAAA", "GTEX-BBBB", "GTEX-CCCC", "GTEX-DDDD"),
           "Uterus": ("GTEX-AAAA", "GTEX-DDDD"),
           "Whole Blood": ("GTEX-AAAA", "GTEX-BBBB", "GTEX-CCCC")}  # fmt: skip
NPTX2_TPM = {"Pituitary": (180.0, 220.0, 260.0, 300.0), "Uterus": (2.0, 4.0),
             "Whole Blood": (0.0, 0.05, 0.1)}  # fmt: skip
# Exon reads per sample: NPTX2's five exons; the last sample of Whole Blood has none.
EXON_READS = (10.0, 10.0, 20.0, 20.0, 140.0)
GTF_DAY = "##description: made-up GENCODE-shaped annotation for the M9 tests\n"


def samples() -> list[tuple[str, str, str]]:
    """(sample ID, tissue, subject) of every made-up sample, in column order."""
    out = []
    for t, (tissue, donors) in enumerate(TISSUES.items()):
        out += [(f"{d}-{t:04d}-SM-{i}", tissue, d) for i, d in enumerate(donors)]
    return out


def _gct(path: Path) -> None:
    """The gene TPM GCT: NPTX2, ACTB, a _PAR_Y row and a filler gene."""
    ids = [s for s, _, _ in samples()]
    nptx2 = [v for t in TISSUES for v in NPTX2_TPM[t]]
    rows = [(THEIRS[NPTX2], "NPTX2", nptx2), (THEIRS[ACTB], "ACTB", [5000.0] * len(ids)),
            ("ENSG00000106236.4_PAR_Y", "NPTX2", [999.0] * len(ids)),
            ("ENSG00000999999.1", "FILLER", [1.0] * len(ids))]  # fmt: skip
    lines = ["#1.2", f"{len(rows)}\t{len(ids)}", "\t".join(["Name", "Description", *ids])]
    lines += ["\t".join([g, d, *(f"{v:g}" for v in values)]) for g, d, values in rows]
    path.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode(), mtime=0))


def _parquet(path: Path) -> None:
    """The exon reads parquet: NPTX2's five exons and one of ACTB, a sample per column."""
    ids = [s for s, _, _ in samples()]
    names = [f"{THEIRS[NPTX2]}_{k}" for k in range(1, 6)] + [f"{THEIRS[ACTB]}_1"]
    columns: dict[str, list] = {"Description": ["NPTX2"] * 5 + ["ACTB"]}
    for i, sample in enumerate(ids):
        scale = 0.0 if i == len(ids) - 1 else 1.0 + i
        columns[sample] = [r * scale for r in EXON_READS] + [100.0]
    columns["Name"] = names
    frame = pl.DataFrame(columns).with_columns(pl.col(ids).cast(pl.Float32))
    frame.write_parquet(path)


def _annotations(folder: Path) -> tuple[Path, Path]:
    """SampleAttributesDS and SubjectPhenotypesDS, with a non-RNA-seq row that must be skipped."""
    rows = ["SAMPID\tSMTSD\tSMAFRZE"] + [f"{s}\t{t}\tRNASEQ" for s, t, _ in samples()]
    rows.append("GTEX-ZZZZ-0001-SM-9\tPituitary\tWGS")
    attributes = folder / "GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"
    attributes.write_text("\n".join(rows) + "\n")
    subjects = folder / "GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt"
    subjects.write_text("SUBJID\tSEX\tAGE\n" + "".join(f"{s}\t{x}\t60-69\n" for s, x in SUBJECTS))
    return attributes, subjects


def _gtf(path: Path) -> None:
    """NPTX2's five MANE Select exons as the real GTF holds them, plus a non-MANE exon."""
    mane = annotation.rest_mane(gzip.decompress(
        (FIXTURES / "ensembl_116_ENSG00000106236_lookup.json.gz").read_bytes()))  # fmt: skip
    assert mane is not None
    attrs = 'gene_id "ENSG00000106236.6"; transcript_id "ENST00000265634.4";'
    lines = [GTF_DAY.rstrip("\n")]
    for k, (a, b) in enumerate(mane.exons, start=1):
        lines.append(f"chr7\tHAVANA\texon\t{a}\t{b}\t.\t+\t.\t{attrs} exon_number {k}; "
                     'tag "MANE_Select";')  # fmt: skip
    lines.append('chr7\tHAVANA\texon\t1\t9\t.\t+\t.\tgene_id "ENSG00000106236.6"; '
                 'transcript_id "ENST00000466102.1"; exon_number 1;')  # fmt: skip
    path.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode(), mtime=0))


def _resource(name: str, path: Path, root: Path) -> str:
    """A [[resources]] entry for a file already in place, as ctx setup registers it."""
    body = path.read_bytes()
    return (f'[[resources]]\nname = "{name}"\nurls = ["http://127.0.0.1/{path.name}"]\n'
            f"bytes = {len(body)}\nsha256 = \"{hashlib.sha256(body).hexdigest()}\"\n"
            f'dest = "{path.relative_to(root).as_posix()}"\n\n')  # fmt: skip


def _tables(root: Path, day: date) -> None:
    """Gather's genes table and the one contrast's M5 results."""
    ids = [OURS[NPTX2], OURS[ACTB], *TWINS, *CALLED, RENAMED[0]]
    names = ["NPTX2", "ACTB", "TWIN", "TWIN", *(v[0] for v in CALLED.values()), RENAMED[1]]
    genes = pl.DataFrame({"gene_id": ids, "gene_name": names})
    described = {"gene_id": ("none", "Ensembl ID."), "gene_name": ("none", "Symbol.")}
    write_table(genes, root / "outputs" / "01_gather", stem("gather", f"{SLUG}_genes", 1, day=day),
                columns=described, data_releases={})  # fmt: skip
    rows = [(OURS[ACTB], 0.4, 0.3, 0.2, "ns")]
    rows += [(g, lfc, lfc * 0.9, p, call) for g, (_, p, lfc, call) in CALLED.items()]
    results = pl.DataFrame(rows, schema=["gene_id", "log2FoldChange", "log2FoldChange_shrunk",
                                         "padj", "call"], orient="row")  # fmt: skip
    columns = {c: ("none", c) for c in results.columns}
    write_table(results, root / "outputs" / "05_de",
                stem("de", f"{SLUG}_a_case_results", 1, day=day), columns=columns,
                data_releases={})  # fmt: skip


def make(root: Path, day: date | None = None) -> Path:
    """Write the made-up root and its dataset file; return the dataset file.

    ``ctx setup`` has not run: the caller registers the files, as a reader would.
    """
    day = day or date.today()
    folder = root / "data" / "raw" / "gtex_v11"
    folder.mkdir(parents=True)
    tpm = folder / "GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz"
    exons = folder / "GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_exon_reads.parquet"
    _gct(tpm)
    _parquet(exons)
    attributes, subjects = _annotations(folder)
    gtf = root / "resources" / "gencode" / "made_up.gtf.gz"
    gtf.parent.mkdir(parents=True)
    _gtf(gtf)
    _tables(root, day)
    spec = root / f"{SLUG}.toml"
    spec.write_text(
        f'[dataset]\nslug = "{SLUG}"\nhpa_cell_lines = {{ group = "{LINES[0]}", lines = '
        + json.dumps(list(LINES[1])) + " }\n\n"
        + _resource("annotation_gtf", gtf, root) + _resource("gtex_gene_tpm", tpm, root)
        + _resource("gtex_exon_reads", exons, root) + _resource("gtex_samples", attributes, root)
        + _resource("gtex_subjects", subjects, root)
        + "[guard]\nids = []\npatterns = []\n\n"
        '[[contrasts]]\nname = "a_case"\nvariable = "condition"\nreference = "control"\n'
        'tested = "case"\nn = [3, 3]\n\n[de]\nfeatured = "a_case"\nbox_genes = ["NPTX2"]\n'
    )  # fmt: skip
    return spec


def responses() -> dict[str, bytes]:
    """NPTX2's responses by URL, as a stub opener serves them; the V10 median TPM is made up."""
    gene = "ENSG00000106236"
    files = {hpa.url("v25", gene, "xml"): "hpa_v25_ENSG00000106236.xml.gz",
             hpa.url("v25", gene, "json"): "hpa_v25_ENSG00000106236.json.gz",
             annotation.lookup_url(gene): "ensembl_116_ENSG00000106236_lookup.json.gz",
             annotation.release_url(): "ensembl_116_info_data.json.gz",
             hpa.search_url("v25", gene, ("g", "eg", LINES[0])):
                 "hpa_v25_search_ENSG00000106236_cell_RNA_leukemia.json.gz",
             hpa.search_url("v25", "NDUFA4L2", ("g", "gs", "eg")):
                 "hpa_v25_search_name_NDUFA4L2.json.gz",
             junctions.gene_url(gene): "gtex_v10_ENSG00000106236_reference_v39.json.gz",
             junctions.junction_url("ENSG00000106236.4"):
                 "gtex_v10_ENSG00000106236_junctions.json.gz",
             junctions.tpm_url("ENSG00000106236.4"):
                 "gtex_v10_ENSG00000106236_median_tpm_madeup.json.gz"}  # fmt: skip
    return {u: gzip.decompress((FIXTURES / f).read_bytes()) for u, f in files.items()}


class Stub:
    """An opener serving fixed bodies by URL, answering 404 to the rest, counting each call."""

    def __init__(self, bodies: dict[str, bytes]) -> None:
        """Serve ``bodies``."""
        self.bodies, self.calls = bodies, []

    def __call__(self, url: str) -> io.BytesIO:
        """The body of ``url``, or an HTTP 404."""
        self.calls.append(url)
        if url not in self.bodies:
            raise urllib.error.HTTPError(url, 404, "Not Found", None, None)  # type: ignore[arg-type]
        return io.BytesIO(self.bodies[url])


def tpm_values() -> np.ndarray:
    """NPTX2's made-up TPM in column order, for the tests to recompute from."""
    return np.array([v for t in TISSUES for v in NPTX2_TPM[t]])


# M9's own samples (E.1c2): eight samples of four individuals in two batches; NPTX2's counts,
# ACTB at 1,000 everywhere. The size factors are planted, not median-of-ratios, so a reader
# that fits its own cannot meet them.
SAMPLES = {
    "s1": ("d1", "B1"),
    "s2": ("d2", "B1"),
    "s3": ("d3", "B1"),
    "s4": ("d1", "B1"),
    "s5": ("d2", "B2"),
    "s6": ("d3", "B2"),
    "s7": ("d4", "B2"),
    "s8": ("d4", "B2"),
}
NPTX2_COUNTS = {"s1": 10, "s2": 20, "s3": 40, "s4": 100, "s5": 400, "s6": 25, "s7": 7, "s8": 9}
# UPA, which M5 tested in a_case, as NPTX2 was not.
UPA_ID = "ENSG00000900011.1"
UPA_COUNTS = {"s1": 5, "s2": 6, "s3": 7, "s4": 50, "s5": 60, "s6": 70, "s7": 8, "s8": 80}
FACTORS = {
    "a_case": {"s1": ("reference", 0.5), "s2": ("reference", 1.0), "s3": ("reference", 2.0),
               "s4": ("tested", 1.0), "s5": ("tested", 4.0), "s6": ("tested", 0.25)},
    "b_case": {"s1": ("reference", 1.0), "s2": ("reference", 1.0), "s3": ("reference", 1.0),
               "s7": ("reference", 2.0), "s4": ("tested", 2.0), "s5": ("tested", 2.0),
               "s6": ("tested", 1.0), "s8": ("tested", 0.5)},
}  # fmt: skip


def own_files(
    root: Path, contrasts: tuple[str, ...] = ("a_case",), day: date | None = None
) -> None:
    """Gather's counts and sample sheet, and M5's size factors of ``contrasts``."""
    day = day or date.today()
    names = list(SAMPLES)
    counts = pl.DataFrame({"gene_id": [OURS[NPTX2], OURS[ACTB], UPA_ID],
                           **{s: [NPTX2_COUNTS[s], 1000, UPA_COUNTS[s]] for s in names}},
                          schema_overrides={s: pl.Int32 for s in names})  # fmt: skip
    gathered = root / "outputs" / "01_gather"
    described = {c: ("none", c) for c in counts.columns}
    write_table(counts, gathered, stem("gather", f"{SLUG}_counts", 1, day=day), columns=described,
                data_releases={})  # fmt: skip
    sheet = pl.DataFrame({"sample": names, "donor": [SAMPLES[s][0] for s in names],
                          "batch": [SAMPLES[s][1] for s in names]})  # fmt: skip
    write_table(sheet, gathered, stem("gather", f"{SLUG}_samplesheet", 1, day=day),
                columns={c: ("none", c) for c in sheet.columns}, data_releases={})  # fmt: skip
    rows = [(c, s, SAMPLES[s][0], g, f) for c in contrasts for s, (g, f) in FACTORS[c].items()]
    factors = pl.DataFrame(rows, schema=["contrast", "sample", "donor", "group", "size_factor"],
                           orient="row")  # fmt: skip
    write_table(factors, root / "outputs" / "05_de", stem("de", f"{SLUG}_size_factors", 1, day=day),
                columns={c: ("none", c) for c in factors.columns}, data_releases={})  # fmt: skip
