#!/usr/bin/env python3
# /// codebase-meta
# exempt: instructor-side builder of tests/fixtures/quant; runs Salmon and R tximport, offline
# ///
"""Build tests/fixtures/quant: Salmon on the chr22_se reads, then R tximport, as nf-core runs it.

Two parts, both offline. ``rnaseq/`` is an nf-core/rnaseq outdir in the layout ``ctx quant``
reads, the per-sample files of copy Q1 included: Salmon 2.7.0 from the pinned ``reads``
environment quantifies the five chr22_se samples with ``-l A`` against the chr22 fixture's
transcripts, its genome slice as decoys; R tximport 1.30.0 in nf-core's cached tximeta image
(never pulled) writes the merged tables with nf-core's ``build_table`` columns. ``planted/`` is
four made-up ``quant.sf`` files whose effective lengths differ by sample, with a gene at zero
in one sample, a gene at zero everywhere and a transcript the tx2gene table lacks, and R's gene
tables of them: the cases that tell tximport's order of operations apart.

Usage:
    pixi run python scripts/python/make_quant_fixture.py --out <empty directory>
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from ctx import tools  # noqa: E402 - the repository's src is put on the path first

CHR22 = ROOT / "tests" / "fixtures" / "chr22"
SE = ROOT / "tests" / "fixtures" / "chr22_se"
PIPELINE_INFO = ROOT / "tests" / "fixtures" / "rnaseq" / "pipeline_info"
SAMPLES = ("se_a", "se_b", "se_c", "se_d", "se_e")
# nf-core/rnaseq 3.27.0's tximport image, cached locally; a missing image stops the build.
IMAGE = "community.wave.seqera.io/library/bioconductor-tximeta_jq:b3a71c30effcaa9d"
MANIFEST = "MANIFEST.sha256"
# The planted set: gene G1's three isoforms are absent from s3, G2 is zero in every sample,
# G3 has one isoform, and T9 is in no gene. Effective lengths differ by sample, as they do
# between paired-end runs, so the arithmetic and geometric means of a gene's lengths differ.
PLANTED_GENES = {"T1": "G1", "T2": "G1", "T3": "G1", "T4": "G2", "T5": "G2", "T6": "G3"}
PLANTED_LENGTHS = {"T1": 3000, "T2": 1500, "T3": 800, "T4": 2000, "T5": 900, "T6": 1200, "T9": 500}
PLANTED_READS = {
    "s1": (300, 50, 20, 0, 0, 400, 7), "s2": (40, 200, 90, 0, 0, 380, 3),
    "s3": (0, 0, 0, 0, 0, 500, 5), "s4": (120, 120, 10, 0, 0, 350, 0),
}  # fmt: skip
PLANTED_FRAGMENT = {"s1": 180.0, "s2": 240.0, "s3": 210.0, "s4": 160.0}

# nf-core's tximport.r, its tables only: transcript tables keyed by tx and gene_id, gene tables
# by gene_id and gene_name; genes in C collation, as its rowsum(reorder = TRUE) leaves them.
R_SCRIPT = r"""
Sys.setlocale("LC_COLLATE", "C")
suppressMessages(library(tximport))
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]; samples <- strsplit(args[2], ",")[[1]]; prefix <- args[3]
files <- setNames(file.path(dir, samples, "quant.sf"), samples)
txi <- tximport(files, type = "salmon", txOut = TRUE)
info <- read.csv(file.path(dir, "tx2gene.tsv"), sep = "\t", header = TRUE, check.names = FALSE)
tx <- data.frame(tx = info[[1]], gene_id = info[[2]], gene_name = info[[3]])
extra <- setdiff(rownames(txi$counts), tx$tx)
tx <- rbind(tx, data.frame(tx = extra, gene_id = extra, gene_name = extra))
tx <- tx[match(rownames(txi$counts), tx$tx), ]
gi <- summarizeToGene(txi, tx2gene = tx[, 1:2])
gi.ls <- summarizeToGene(txi, tx2gene = tx[, 1:2], countsFromAbundance = "lengthScaledTPM")
gi.s <- summarizeToGene(txi, tx2gene = tx[, 1:2], countsFromAbundance = "scaledTPM")
real <- setdiff(rownames(gi$counts), extra)
genes <- unique(tx[, 2:3]); genes <- genes[match(real, genes$gene_id), ]
w <- function(rows, m, suffix) write.table(data.frame(rows, m, check.names = FALSE),
    file.path(dir, paste0(prefix, ".", suffix)), sep = "\t", quote = FALSE, row.names = FALSE)
w(tx[, 1:2], txi$abundance, "transcript_tpm.tsv")
w(tx[, 1:2], txi$counts, "transcript_counts.tsv")
w(tx[, 1:2], txi$length, "transcript_lengths.tsv")
w(genes, gi$length[real, , drop = FALSE], "gene_lengths.tsv")
w(genes, gi$abundance[real, , drop = FALSE], "gene_tpm.tsv")
w(genes, gi$counts[real, , drop = FALSE], "gene_counts.tsv")
w(genes, gi.ls$counts[real, , drop = FALSE], "gene_counts_length_scaled.tsv")
w(genes, gi.s$counts[real, , drop = FALSE], "gene_counts_scaled.tsv")
cat(as.character(packageVersion("tximport")), "\n")
"""


def require_empty(out: Path) -> None:
    """Refuse a directory that holds anything: the fixture is written whole or not at all."""
    if out.exists() and any(out.iterdir()):
        sys.exit(f"{out} is not empty")
    out.mkdir(parents=True, exist_ok=True)


def tx2gene(path: Path) -> None:
    """nf-core's tx2gene.tsv for the chr22 fixture: transcript_id, gene_id, gene_name."""
    rows = ["transcript_id\tgene_id\tgene_name"]
    for line in (CHR22 / "transcripts.fa").read_text().splitlines():
        if line.startswith(">"):
            name, *fields = line[1:].split()
            tags = dict(f.split("=", 1) for f in fields)
            rows.append(f"{name}\t{tags['gene_id']}\t{tags['gene_name']}")
    path.write_text("\n".join(rows) + "\n")


def salmon(work: Path, salmon_dir: Path) -> None:
    """Index the transcripts with the genome slice as decoys, then quantify every sample."""
    shutil.copy(CHR22 / "transcripts.fa", work / "transcripts.fa")
    gentrome = work / "gentrome.fa"
    gentrome.write_text((CHR22 / "transcripts.fa").read_text() + (CHR22 / "slice.fa").read_text())
    decoy = (CHR22 / "slice.fa").read_text().splitlines()[0][1:].split()[0]
    (work / "decoys.txt").write_text(decoy + "\n")
    tools.run("salmon", ["index", "-t", "gentrome.fa", "-d", "decoys.txt", "-i", "index",
                         "-k", "25", "-p", "1"], cwd=work, log=work / "index.log")  # fmt: skip
    for sample in SAMPLES:
        shutil.copy(SE / f"{sample}.fq", work / f"{sample}.fq")
        tools.run("salmon", ["quant", "-i", "index", "-l", "A", "-r", f"{sample}.fq", "-p", "1",
                             "-o", sample], cwd=work, log=work / f"{sample}.log")  # fmt: skip
        kept = salmon_dir / sample
        for name in ("quant.sf", "cmd_info.json", "lib_format_counts.json",
                     "aux_info/meta_info.json", "aux_info/ambig_info.tsv",
                     "libParams/flenDist.txt"):  # fmt: skip
            (kept / name).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(work / sample / name, kept / name)


def r_tximport(folder: Path, samples: tuple[str, ...], prefix: str) -> str:
    """R tximport in nf-core's image on ``folder``; return its version."""
    script = folder / "tximport.R"
    script.write_text(R_SCRIPT)
    done = subprocess.run(
        ["docker", "run", "--rm", "--pull", "never", "-v", f"{folder}:/work", IMAGE, "Rscript",
         "/work/tximport.R", "/work", ",".join(samples), prefix],
        capture_output=True, text=True, check=False,
    )  # fmt: skip
    script.unlink()
    if done.returncode:
        sys.exit(f"R tximport failed:\n{done.stderr[-2000:]}")
    return done.stdout.strip().splitlines()[-1]


def planted(folder: Path) -> None:
    """Four made-up quant.sf files and their tx2gene table."""
    folder.mkdir()
    lines = ["transcript_id\tgene_id\tgene_name", *(f"{t}\t{g}\t{g}" for t, g in
                                                   PLANTED_GENES.items())]  # fmt: skip
    (folder / "tx2gene.tsv").write_text("\n".join(lines) + "\n")
    names = list(PLANTED_LENGTHS)
    for sample, reads in PLANTED_READS.items():
        effective = [PLANTED_LENGTHS[t] - PLANTED_FRAGMENT[sample] for t in names]
        rate = [n / e for n, e in zip(reads, effective, strict=True)]
        tpm = [1e6 * r / sum(rate) for r in rate]
        rows = ["Name\tLength\tEffectiveLength\tTPM\tNumReads"]
        rows += [f"{t}\t{PLANTED_LENGTHS[t]}\t{e:.3f}\t{p:.6f}\t{n:.3f}"
                 for t, e, p, n in zip(names, effective, tpm, reads, strict=True)]  # fmt: skip
        (folder / sample).mkdir()
        (folder / sample / "quant.sf").write_text("\n".join(rows) + "\n")


def manifest(out: Path) -> None:
    """SHA-256 of every file, sorted by path."""
    lines = [f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out).as_posix()}"
             for p in sorted(out.rglob("*")) if p.is_file() and p.name != MANIFEST]  # fmt: skip
    (out / MANIFEST).write_text("\n".join(lines) + "\n", newline="\n")


def write(out: Path, work: Path) -> None:
    """Build both parts into ``out``, Salmon's scratch in ``work``."""
    require_empty(out)
    salmon_dir = out / "rnaseq" / "salmon"
    salmon_dir.mkdir(parents=True)
    shutil.copytree(PIPELINE_INFO, out / "rnaseq" / "pipeline_info")
    tx2gene(salmon_dir / "tx2gene.tsv")
    salmon(work, salmon_dir)
    version = r_tximport(salmon_dir, SAMPLES, "salmon.merged")
    planted(out / "planted")
    r_tximport(out / "planted", tuple(PLANTED_READS), "r")
    meta = json.loads((salmon_dir / SAMPLES[0] / "aux_info" / "meta_info.json").read_text())
    (out / "README.md").write_text(README.format(version=version, salmon=meta["salmon_version"],
                                                 image=IMAGE))  # fmt: skip
    manifest(out)


README = """# quant fixture: Salmon and R tximport on chr22_se

Written by `scripts/python/make_quant_fixture.py`; never edit it by hand. Salmon {salmon}
(the pinned `reads` environment) quantified the five `tests/fixtures/chr22_se` samples with
`-l A` against `tests/fixtures/chr22/transcripts.fa`, its `slice.fa` as decoys, k 25. R
tximport {version} in `{image}` wrote the merged tables with nf-core/rnaseq's columns.

- `rnaseq/` is an nf-core/rnaseq outdir as `ctx quant` reads it. `pipeline_info/` is the
  yeast fixture's (`tests/fixtures/rnaseq`), so the release check passes; its params name
  the yeast genome, not this one. `salmon/` holds what copy Q1 holds per sample: `quant.sf`,
  `cmd_info.json`, `lib_format_counts.json`, `aux_info/meta_info.json`,
  `aux_info/ambig_info.tsv` and `libParams/flenDist.txt`, and the merged tables.
- `planted/` holds four made-up `quant.sf` files, their `tx2gene.tsv` and R's gene tables
  (`r.gene_*.tsv`): effective lengths differ by sample, G1 is zero in s3, G2 is zero
  everywhere and T9 is in no gene.

A rebuild repeats every count; `meta_info.json` records run times, so its bytes differ.
"""


def main(argv: list[str] | None = None) -> int:
    """Build the fixture into an empty directory."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True, help="An empty scratch directory.")
    args = parser.parse_args(argv)
    require_empty(args.work)
    write(args.out.resolve(), args.work.resolve())
    return 0


if __name__ == "__main__":
    sys.exit(main())
