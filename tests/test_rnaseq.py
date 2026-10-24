"""ctx.gather.rnaseq on the real fixture, on tmp copies with one fault, and through the CLI."""

import hashlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from ctx.gather.fetchngs import read as read_fetchngs
from ctx.gather.rnaseq import (
    ALIGNED_FILES,
    ALIGNED_MATRIX,
    ALIGNER,
    MATRIX_FILE,
    OPTIONAL_FILES,
    READ_FILES,
    RNASEQ_VERSION,
    VERSIONS_FILE,
    RnaseqError,
    aligned,
    read,
    to_counts,
)

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "rnaseq"
FETCHNGS = ROOT / "tests" / "fixtures" / "fetchngs" / "download"
CTX = Path(sys.executable).parent / "ctx"


def fixture_samples() -> set[str]:
    with (FIXTURE / MATRIX_FILE).open() as handle:
        return set(handle.readline().rstrip("\n").split("\t")[2:])


def copy(tmp_path: Path) -> Path:
    target = tmp_path / "rnaseq"
    shutil.copytree(FIXTURE, target)
    return target


def test_round_then_cast():
    values = np.array([[0.5, 1.5], [2.5, 2.6]])
    counts, fractional = to_counts(values)
    # Half to even: 0.5 -> 0, 1.5 -> 2, 2.5 -> 2. Truncation would give 0, 1, 2, 2.
    assert counts.tolist() == [[0, 2], [2, 3]]
    assert counts.dtype == np.int32
    assert fractional == 4
    assert to_counts(np.array([[2.0**31]]))[0].dtype == np.int64
    with pytest.raises(RnaseqError, match="negative"):
        to_counts(np.array([[-1.0]]))


def test_reads_fixture():
    counts = read(FIXTURE, fixture_samples())
    raw = pl.read_csv(FIXTURE / MATRIX_FILE, separator="\t")
    samples = raw.columns[2:]
    assert counts.matrix.columns == ["gene_id", *samples]
    assert counts.matrix.height == raw.height > 0
    assert all(counts.matrix[s].dtype == pl.Int32 for s in samples)
    assert np.array_equal(
        counts.matrix.select(samples).to_numpy(), np.rint(raw.select(samples).to_numpy())
    )
    assert counts.provenance["nf-core/rnaseq"] == RNASEQ_VERSION
    assert counts.provenance["salmon"] == "2.7.0"


def test_sample_set_mismatch_lists_all():
    samples = (fixture_samples() - {"WT_REP1"}) | {"NOT_THERE"}
    with pytest.raises(RnaseqError) as raised:
        read(FIXTURE, samples)
    assert "sample NOT_THERE: in the sheet, not in the matrix" in str(raised.value)
    assert "sample WT_REP1: in the matrix, not in the sheet" in str(raised.value)


def test_other_release_refused(tmp_path):
    outdir = copy(tmp_path)
    versions = outdir / VERSIONS_FILE
    versions.write_text(versions.read_text().replace(f"v{RNASEQ_VERSION}", "v3.26.0"))
    with pytest.raises(RnaseqError, match=r"3\.26\.0 wrote"):
        read(outdir, fixture_samples())


def test_params_disagree_fails(tmp_path):
    outdir = copy(tmp_path)
    (params,) = (outdir / "pipeline_info").glob("params_*.json")
    changed = params.read_text().replace('"gencode": false', '"gencode": true')
    assert changed != params.read_text()
    params.with_name("params_2099-01-01_00-00-00.json").write_text(changed)
    with pytest.raises(RnaseqError, match=r"disagree on \['gencode'\]"):
        read(outdir, fixture_samples())


STAR = ROOT / "tests" / "fixtures" / "rnaseq_star"


def _aligned_copy(tmp_path: Path) -> Path:
    target = tmp_path / "rnaseq_star"
    shutil.copytree(STAR, target)
    return target


def test_reads_the_aligned_fixture():
    """An outdir aligned with STAR is read as a Salmon-only one, its aligned matrix recorded."""
    counts = read(STAR, fixture_samples())
    assert counts.provenance["params aligner"] == ALIGNER
    assert counts.provenance["params skip_alignment"] == "False"
    assert f"sha256:{ALIGNED_MATRIX}" in counts.provenance
    assert aligned(STAR) and not aligned(FIXTURE)


def test_an_aligned_copy_missing_a_table_is_refused(tmp_path):
    """Mutation: a cluster copy that missed one MultiQC table passes gather and stops M3 later.
    Every missing file is named, here one table and the aligned matrix."""
    outdir = _aligned_copy(tmp_path)
    (outdir / ALIGNED_FILES[-1]).unlink()
    (outdir / ALIGNED_MATRIX).unlink()
    with pytest.raises(RnaseqError) as raised:
        read(outdir, fixture_samples())
    assert ALIGNED_FILES[-1] in str(raised.value) and ALIGNED_MATRIX in str(raised.value)


def test_a_table_multiqc_may_not_write_is_not_required(tmp_path):
    """MultiQC writes the overrepresented-sequences table only when a sample has one; refusing
    its absence stopped M3 on BRCA-immune-cells-macrophages (2026-10-24). Every other table
    stays required, so the optional list cannot grow by accident into the whole list."""
    outdir = _aligned_copy(tmp_path)
    for name in OPTIONAL_FILES:
        (outdir / name).unlink()
    assert read(outdir, fixture_samples()).provenance["params aligner"] == ALIGNER
    assert len(OPTIONAL_FILES) == 1 and set(OPTIONAL_FILES) < set(ALIGNED_FILES)


def test_another_aligner_is_refused(tmp_path):
    """Mutation: an outdir aligned by HISAT2 read as STAR's; ctx reads star_salmon/ alone."""
    outdir = _aligned_copy(tmp_path)
    for params in (outdir / "pipeline_info").glob("params_*.json"):
        params.write_text(params.read_text().replace('"aligner": "star_salmon"',
                                                     '"aligner": "hisat2"'))  # fmt: skip
    with pytest.raises(RnaseqError, match="aligner 'hisat2'"):
        read(outdir, fixture_samples())


def test_launches_disagreeing_on_the_aligner_are_refused(tmp_path):
    """Mutation: AGREE without aligner; a resume that switched aligners mixes two runs."""
    outdir = _aligned_copy(tmp_path)
    last = sorted((outdir / "pipeline_info").glob("params_*.json"))[-1]
    last.write_text(last.read_text().replace('"aligner": "star_salmon"', '"aligner": "star_rsem"'))
    with pytest.raises(RnaseqError, match=r"disagree on \['aligner'\]"):
        read(outdir, fixture_samples())


def test_reads_only_read_files(tmp_path):
    only = tmp_path / "only"
    for pattern in READ_FILES:
        for path in FIXTURE.glob(pattern):
            target = only / path.relative_to(FIXTURE)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, target)
    full, reduced = read(FIXTURE, fixture_samples()), read(only, fixture_samples())
    assert full.matrix.equals(reduced.matrix)
    assert full.genes.equals(reduced.genes)


# --- through the CLI, on the fetchngs download fixture's samples ----------------------------


def synthetic_outdir(tmp_path: Path) -> Path:
    """The fixture's provenance with a matrix whose columns are the fetchngs fixture's samples."""
    outdir = copy(tmp_path)
    samples = sorted(set(read_fetchngs(FETCHNGS).runs["sample"]))
    rng = np.random.default_rng(0)
    frame = pl.DataFrame(
        {"gene_id": [f"g{i}" for i in range(30)], "gene_name": [f"G{i}" for i in range(30)]}
    ).with_columns([pl.Series(s, rng.gamma(2.0, 50.0, 30)) for s in samples])
    frame.write_csv(outdir / MATRIX_FILE, separator="\t")
    return outdir


def _gather(root: Path, rnaseq_dir: Path) -> subprocess.CompletedProcess[str]:
    root.mkdir(exist_ok=True)
    env = {k: v for k, v in os.environ.items() if not k.startswith(("SLURM_", "CTX_"))}
    args = ["--root", str(root), "gather", "--dataset", "counted",
            "--fetchngs-dir", str(FETCHNGS), "--rnaseq-dir", str(rnaseq_dir)]  # fmt: skip
    return subprocess.run(
        [str(CTX), *args], cwd=root, env=env, capture_output=True, text=True, check=False
    )


def _hashes(directory: Path) -> dict[str, tuple[str, int]]:
    return {
        p.relative_to(directory).as_posix(): (
            hashlib.sha256(p.read_bytes()).hexdigest(),
            p.stat().st_mtime_ns,
        )
        for p in directory.rglob("*")
        if p.is_file()
    }


def test_cli_writes_counts_and_leaves_outdir_untouched(tmp_path):
    outdir = synthetic_outdir(tmp_path)
    before = _hashes(outdir)
    done = _gather(tmp_path / "root", outdir)
    assert done.returncode == 0, done.stderr
    out = tmp_path / "root" / "outputs" / "01_gather"
    (counts,) = out.glob("*_gather_counted_counts_v1.parquet")
    frame = pl.read_parquet(counts)
    assert frame.width == 1 + 4 and frame.height == 30
    readme = next(out.glob("*_counted_counts_v1_README.md")).read_text()
    assert f"sha256:{MATRIX_FILE}" in readme and "nf-core/rnaseq" in readme
    assert _hashes(outdir) == before


def test_cli_bad_outdir_writes_nothing(tmp_path):
    outdir = synthetic_outdir(tmp_path)
    frame = pl.read_csv(outdir / MATRIX_FILE, separator="\t")
    frame.drop(frame.columns[-1]).write_csv(outdir / MATRIX_FILE, separator="\t")
    done = _gather(tmp_path / "root", outdir)
    assert done.returncode == 1
    assert "in the sheet, not in the matrix" in done.stderr
    assert not (tmp_path / "root" / "outputs").exists()
