"""A cluster alignment folder as ED.6 block D copies it, built from the chr22_se local run.

``ctx process --align-dir`` reads ``data/raw/<slug>_align/``: tasks.tsv, versions.txt,
refs.txt, ``runs/<sample>/{fastp,star,featurecounts,time}/``, block E's
``featured/<sample>/{untrimmed,cut}/`` and ``reference/{main,cut}/`` (the plan,
docs/plans/2026-10-22_ed6_m2_m3_core_cluster_blocks.md). Every file comes from a real tool run
on the fixture except two, which no tool here can write:

- The pinned osx-arm64 STAR maps nothing (D-036), so each ``Log.final.out`` and
  ``ReadsPerGene.out.tab`` is written in STAR 2.7.11b's layout from the local HISAT2 BAM of the
  same reads: input, unique and multimapped reads from HISAT2's own summary, its unaligned
  reads as STAR's "too short" (STAR's other unmapped kinds and chimeric reads are written 0),
  mapped length, splices, junction motifs and indel and mismatch rates from the BAM's CIGAR
  and NM tags, and the gene rows from featureCounts at -s 0, 1 and 2 on that BAM.
- GNU time's ``-v`` files: the Mac's /usr/bin/time prints another layout, so the two lines
  ctx reads are written from each featureCounts run's own wall time and ``wait4`` peak memory.

featureCounts runs exactly as block B1 calls it, one sample and setting at a time, and the
references are ``reference.prepare``'s, the code block E runs. The fixture is single-end only.
"""

import os
import re
import shutil
import subprocess
import time
from pathlib import Path

import polars as pl
from process_root import SAMPLES, SLICE

from ctx import tools
from ctx.process import fc_compact, reference

ROOT = Path(__file__).resolve().parents[1]
GTF = ROOT / "tests" / "fixtures" / "chr22" / "slice.gtf"
FEATURED = SAMPLES[0]
INDEX, CUT_INDEX = "star_2.7.11b_gencode50_sjdb100", "star_2.7.11b_gencode50_half"
# Block B1's 12 settings: (file name, featureCounts options after -s).
VARIANTS = (("", ()), ("O_", ("-O",)), ("O_fraction_", ("-O", "--fraction")), ("M_", ("-M",)))
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
MOTIFS = {("GT", "AG"): "GT/AG", ("CT", "AC"): "GT/AG", ("GC", "AG"): "GC/AG",
          ("CT", "GC"): "GC/AG", ("AT", "AC"): "AT/AC", ("GT", "AT"): "AT/AC"}  # fmt: skip


def versions(dest: Path) -> None:
    """versions.txt as block B1 writes it, from this machine's pinned tools."""

    def out(name: str, *args: str) -> str:
        done = subprocess.run([str(tools.require(name)), *args], capture_output=True, text=True,
                              check=False)  # fmt: skip
        return done.stdout + done.stderr

    lines = [f"STAR {out('STAR', '--version').strip()}", out("fastp", "--version").strip(),
             out("fastqc", "--version").strip(), out("seqkit", "version").strip(),
             out("samtools", "--version").splitlines()[0],
             next(line for line in out("featureCounts", "-v").splitlines()
                  if "featureCounts v" in line).strip()]  # fmt: skip
    dest.write_text("\n".join(lines) + "\n")


def timed(args: list[str], cwd: Path, record: Path) -> None:
    """Run a tool and write its wall time and peak memory in GNU time -v's two lines."""
    started = time.perf_counter()
    child = subprocess.Popen(args, cwd=cwd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    _, status, usage = os.wait4(child.pid, 0)
    assert status == 0, args
    wall = time.perf_counter() - started
    # macOS reports ru_maxrss in bytes; GNU time prints kilobytes.
    record.parent.mkdir(parents=True, exist_ok=True)
    record.write_text(
        f'\tCommand being timed: "{" ".join(args)}"\n'
        f"\tElapsed (wall clock) time (h:mm:ss or m:ss): 0:{wall:05.2f}\n"
        f"\tMaximum resident set size (kbytes): {usage.ru_maxrss // 1024}\n"
    )


def feature_counts(bam: Path, gtf: Path, out: Path, options: tuple[str, ...]) -> None:
    """featureCounts on one BAM, run in its folder as block B1 runs it, timed."""
    out.parent.mkdir(parents=True, exist_ok=True)
    args = [str(tools.require("featureCounts")), "-T", "1", "-a", str(gtf), *options, "-o",
            str(out), bam.name]  # fmt: skip
    timed(args, bam.parent, out.parent.parent / "time" / f"fc_{out.name.split('.')[0]}.time")


def hisat2_reads(summary: Path) -> dict[str, int]:
    """Total, once, more than once and unaligned reads of HISAT2's --new-summary, single-end."""
    text = summary.read_text()
    number = {k: int(re.search(rf"{p}: (\d+)", text)[1]) for k, p in (  # type: ignore[index]
        ("total", "Total reads"), ("unique", "Aligned 1 time"), ("multi", r"Aligned >1 times"),
        ("unmapped", "Aligned 0 time"))}  # fmt: skip
    return number


def bam_details(bam: Path, introns: set[tuple[str, int, int]], genome: dict[str, str]
                ) -> dict[str, float]:  # fmt: skip
    """Mapped length, splices by motif and annotation, and per-base rates of unique reads."""
    view = subprocess.run([str(tools.require("samtools")), "view", "-F", "0x904", str(bam)],
                          capture_output=True, text=True, check=True).stdout  # fmt: skip
    sums = dict.fromkeys(("n", "matched", "mismatched", "deleted", "inserted", "deletions",
                          "insertions"), 0)  # fmt: skip
    splices = dict.fromkeys(("total", "annotated", "GT/AG", "GC/AG", "AT/AC", "non-canonical"), 0)
    for line in view.splitlines():
        fields = line.split("\t")
        if "NH:i:1" not in fields[11:]:
            continue
        sums["n"] += 1
        at, indels = int(fields[3]), 0
        for length, op in ((int(a), b) for a, b in CIGAR.findall(fields[5])):
            if op == "M":
                sums["matched"] += length
            elif op in "DI":
                kind = "deleted" if op == "D" else "inserted"
                sums[kind] += length
                sums["deletions" if op == "D" else "insertions"] += 1
                indels += length
            elif op == "N":
                start, end = at, at + length - 1
                seq = genome[fields[2]]
                splices["total"] += 1
                splices["annotated"] += (fields[2], start, end) in introns
                motif = (seq[start - 1 : start + 1].upper(), seq[end - 2 : end].upper())
                splices[MOTIFS.get(motif, "non-canonical")] += 1
            if op in "MDN":
                at += length
        nm = next(int(f[5:]) for f in fields[11:] if f.startswith("NM:i:"))
        sums["mismatched"] += nm - indels
    m = sums["matched"] or 1
    return {"mapped_length": (sums["matched"] + sums["inserted"]) / max(sums["n"], 1),
            "mismatch": 100 * sums["mismatched"] / m, "deletion": 100 * sums["deleted"] / m,
            "deletion_length": sums["deleted"] / max(sums["deletions"], 1),
            "insertion": 100 * sums["inserted"] / m,
            "insertion_length": sums["inserted"] / max(sums["insertions"], 1),
            **{f"splices_{k}": v for k, v in splices.items()}}  # fmt: skip


def introns_of(gtf: Path) -> set[tuple[str, int, int]]:
    """Every annotated intron, (contig, first base, last base), from each transcript's exons."""
    exons: dict[str, list[tuple[str, int, int]]] = {}
    for line in gtf.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) > 8 and fields[2] == "exon":  # noqa: PLR2004 - GTF's nine columns
            transcript = re.search(r'transcript_id "([^"]+)"', fields[8])[1]  # type: ignore[index]
            exons.setdefault(transcript, []).append((fields[0], int(fields[3]), int(fields[4])))
    found = set()
    for spans in exons.values():
        ordered = sorted(spans, key=lambda s: s[1])
        found |= {(a[0], a[2] + 1, b[1] - 1) for a, b in zip(ordered, ordered[1:], strict=False)}
    return found


def log_final(path: Path, reads: dict[str, int], details: dict[str, float],
              read_length: float) -> None:  # fmt: skip
    """A Log.final.out in STAR 2.7.11b's layout, its numbers measured as the docstring says."""
    total = reads["total"]

    def pct(n: float) -> str:
        return f"{100 * n / total:.2f}%" if total else "0.00%"

    rows = [
        ("Started job on", "Oct 22 12:00:00"), ("Started mapping on", "Oct 22 12:00:01"),
        ("Finished on", "Oct 22 12:00:02"), ("Mapping speed, Million of reads per hour", "3.60"),
        None, ("Number of input reads", total), ("Average input read length", round(read_length)),
        ("UNIQUE READS:", None), ("Uniquely mapped reads number", reads["unique"]),
        ("Uniquely mapped reads %", pct(reads["unique"])),
        ("Average mapped length", f"{details['mapped_length']:.2f}"),
        ("Number of splices: Total", details["splices_total"]),
        ("Number of splices: Annotated (sjdb)", details["splices_annotated"]),
        ("Number of splices: GT/AG", details["splices_GT/AG"]),
        ("Number of splices: GC/AG", details["splices_GC/AG"]),
        ("Number of splices: AT/AC", details["splices_AT/AC"]),
        ("Number of splices: Non-canonical", details["splices_non-canonical"]),
        ("Mismatch rate per base, %", f"{details['mismatch']:.2f}%"),
        ("Deletion rate per base", f"{details['deletion']:.2f}%"),
        ("Deletion average length", f"{details['deletion_length']:.2f}"),
        ("Insertion rate per base", f"{details['insertion']:.2f}%"),
        ("Insertion average length", f"{details['insertion_length']:.2f}"),
        ("MULTI-MAPPING READS:", None),
        ("Number of reads mapped to multiple loci", reads["multi"]),
        ("% of reads mapped to multiple loci", pct(reads["multi"])),
        ("Number of reads mapped to too many loci", 0),
        ("% of reads mapped to too many loci", pct(0)),
        ("UNMAPPED READS:", None),
        ("Number of reads unmapped: too many mismatches", 0),
        ("% of reads unmapped: too many mismatches", pct(0)),
        ("Number of reads unmapped: too short", reads["unmapped"]),
        ("% of reads unmapped: too short", pct(reads["unmapped"])),
        ("Number of reads unmapped: other", 0), ("% of reads unmapped: other", pct(0)),
        ("CHIMERIC READS:", None), ("Number of chimeric reads", 0), ("% of chimeric reads", pct(0)),
    ]  # fmt: skip
    lines = []
    for row in rows:
        if row is None:
            lines.append("")
        elif row[1] is None:
            lines.append(f"{row[0]:>48}")
        else:
            lines.append(f"{row[0]:>48} |\t{row[1]}")
    path.write_text("\n".join(lines) + "\n")


def reads_per_gene(path: Path, reads: dict[str, int], counts: Path) -> None:
    """ReadsPerGene.out.tab from featureCounts at -s 0, 1 and 2 in ``counts`` (s<k> files).

    STAR's columns 2 to 4 are -s 0, 1 and 2; each must sum to the input reads.
    """
    genes, rows = None, {}
    for k in (0, 1, 2):
        table = pl.read_csv(counts / f"s{k}.featureCounts.txt", separator="\t", comment_prefix="#")
        genes = table["Geneid"].to_list()
        rows[k] = table[table.columns[-1]].to_list()
        lines = (counts / f"s{k}.featureCounts.txt.summary").read_text().splitlines()[1:]
        summary = dict(line.split("\t") for line in lines)
        rows[k] = [reads["unmapped"], reads["multi"], int(summary["Unassigned_NoFeatures"]),
                   int(summary["Unassigned_Ambiguity"]), *rows[k]]  # fmt: skip
        assert sum(rows[k]) == reads["total"], (path, k, sum(rows[k]), reads)
    names = ["N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous", *(genes or [])]
    path.write_text("".join(f"{n}\t{rows[0][i]}\t{rows[1][i]}\t{rows[2][i]}\n"
                            for i, n in enumerate(names)))  # fmt: skip


def star_files(align: Path, sample: str, dest: Path, gtf: Path,  # noqa: PLR0913, PLR0917 - one each
               keep: tuple[tuple[str, tuple[str, ...]], ...], context: dict,
               annotation: Path | None = None) -> None:  # fmt: skip
    """``<dest>/star/`` and ``<dest>/featurecounts/`` for one sample of one local alignment.

    ``keep`` names the featureCounts files kept at each strand, by prefix and options: block
    B1's 12 (``VARIANTS``), block E's ``cut_``, or none for the untrimmed alignment. Each kept
    file is compacted as the blocks compact it (``fc_compact``), its annotation columns kept once
    in ``annotation``.
    """
    star, kept = dest / "star", dest / "featurecounts"
    star.mkdir(parents=True, exist_ok=True)
    work = context["work"] / f"{dest.name}_{sample}"
    work.mkdir(parents=True)
    bam = work / f"{sample}.bam"
    shutil.copyfile(align / f"{sample}.bam", bam)
    for k in (0, 1, 2):
        feature_counts(bam, gtf, work / "plain" / f"s{k}.featureCounts.txt", ("-s", str(k)))
        for prefix, options in keep:
            matrix = kept / f"{prefix}s{k}.featureCounts.txt"
            feature_counts(bam, gtf, matrix, ("-s", str(k), *options))
            assert annotation is not None
            fc_compact.compact(str(matrix), str(annotation))
    reads = hisat2_reads(align / f"{sample}.hisat2.txt")
    details = bam_details(bam, context["introns"], context["genome"])
    log_final(star / f"{sample}.Log.final.out", reads, details, context["length"][sample])
    reads_per_gene(star / f"{sample}.ReadsPerGene.out.tab", reads, work / "plain")
    for name in ("flagstat", "idxstats"):
        shutil.copyfile(align / f"{sample}.{name}", star / f"{sample}.{name}")


def build(local: Path, root: Path, slug: str) -> Path:
    """``<root>/data/raw/<slug>_align/`` from the local HISAT2 run under ``local``.

    Args:
        local: A ctx root where ``ctx process --aligner hisat2 --settings cut,untrimmed`` ran
            on the fixture (``conftest.process_root``).
        root: The root the cluster folder goes into.
        slug: The dataset's slug.

    Returns:
        The folder.
    """
    out = local / "outputs" / "03_process"
    folders = dict(pl.read_csv(next(out.glob("*_folders_v1.csv"))).select("step", "folder").rows())
    tools_ds = out / "tools" / "chrse"
    sheet = pl.read_parquet(next((root / "outputs" / "01_gather").glob(
        f"*_{slug}_samplesheet_v1.parquet")))  # fmt: skip
    dest = root / "data" / "raw" / f"{slug}_align"
    dest.mkdir(parents=True)
    runs = dict(sheet.select("sample", "run_accession").rows())
    pl.DataFrame(
        {
            "sample": list(runs),
            "mates": ["1"] * len(runs),
            "runs": list(runs.values()),
            "r1": [f"{s}.fastq.gz" for s in runs],
            "r2": ["-"] * len(runs),
            "bytes": ["1000"] * len(runs),
        }
    ).write_csv(dest / "tasks.tsv", separator="\t")
    versions(dest / "versions.txt")
    (dest / "refs.txt").write_text(f"REFS=/cluster/refs\nIDX=/cluster/refs/{INDEX}\n"
                                   f"GTF=/cluster/refs/unzipped/{GTF.name}\n")  # fmt: skip
    work = root / "cluster_work"
    refs = {name: reference.prepare(SLICE, GTF, work / "refs", half=half)
            for name, half in (("main", False), ("cut", True))}  # fmt: skip
    for name, ref in refs.items():
        (dest / "reference" / name).mkdir(parents=True)
        for f in ("reference.json", "ctx_manifest.json"):
            shutil.copyfile(ref.folder / f, dest / "reference" / name / f)
    (dest / "reference" / "inputs.txt").write_text(f"GENOME={SLICE.name}\nCUT_INDEX={CUT_INDEX}\n")
    genome = dict(reference.records(SLICE))
    trim = tools_ds / folders["trim"] / "fastp"
    length = {}
    for sample, run in runs.items():
        (dest / "runs" / sample / "fastp").mkdir(parents=True)
        shutil.copyfile(trim / f"{run}.fastp.json", dest / "runs" / sample / "fastp" /
                        f"{run}.fastp.json")  # fmt: skip
        report = pl.read_json(trim / f"{run}.fastp.json")
        length[sample] = report["summary"][0]["after_filtering"]["read1_mean_length"]
    context = {"work": work, "introns": introns_of(GTF), "genome": genome, "length": length}
    for sample in runs:
        star_files(tools_ds / folders["align_main"], sample, dest / "runs" / sample, GTF,
                   VARIANTS, context, dest / fc_compact.ANNOTATION)  # fmt: skip
    featured = dest / "featured" / FEATURED
    star_files(tools_ds / folders["align_untrimmed"], FEATURED, featured / "untrimmed", GTF,
               (), context)  # fmt: skip
    star_files(tools_ds / folders["align_cut"], FEATURED, featured / "cut", refs["cut"].gtf,
               (("cut_", ()),), context, featured / "cut" / fc_compact.ANNOTATION)  # fmt: skip
    shutil.rmtree(work)
    return dest
