"""M1c curation on the metadata-only fetchngs fixture, with an SDRF written from its own runs.

No accession is typed here: every run name comes from the fixture's samplesheet.
"""

import hashlib
import os
import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

from ctx.dataset import DatasetError, load
from ctx.gather.curation import CurationError, curate
from ctx.gather.fetchngs import read
from ctx.gather.harmonise import DEPTH_FILTERS, _depth_standard
from ctx.io.naming import figure_files, lint
from ctx.report import text
from ctx.report.builder import unpacked

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "fetchngs" / "metadata_only"
SPEC = ROOT / "tests" / "fixtures" / "datasets" / "synthetic.toml"
CTX = Path(sys.executable).parent / "ctx"
HEADER = [
    "Source Name",
    "Characteristics[individual]",
    "Characteristics[cell type]",
    "Characteristics[condition]",
    "Characteristics[time point]",
    # Repeated but not configured, as Protocol REF is in a real SDRF: allowed.
    "Protocol REF",
    "Protocol REF",
    "Comment[ENA_RUN]",
]


@pytest.fixture(scope="module")
def runs():
    return read(FIXTURE).runs


def sdrf_rows(runs: pl.DataFrame) -> list[dict[str, str]]:
    """One SDRF row per run. Runs 0 and 1 share donor, cell type and condition, and their
    time points sort the opposite way to their accessions; run 1's has two spaces."""
    rows = []
    for i, (run, title) in enumerate(runs.select("run_accession", "sample_title").rows()):
        rows.append({
            "Source Name": f" {title} ",
            "Characteristics[individual]": "d0" if i < 2 else f"d{i}",
            "Characteristics[cell type]": "type A" if i < 2 or i % 2 else "type B",
            "Characteristics[condition]": "Control" if i < 7 else "Case",
            "Characteristics[time point]": {0: "10th follow up", 1: "2nd  follow up"}.get(i, "na"),
            "Protocol REF": "P-1",
            "Comment[ENA_RUN]": run,
        })  # fmt: skip
    return rows


def sdrf_text(rows: list[dict[str, str]], header: list[str] = HEADER) -> str:
    lines = ["\t".join(header)]
    lines += ["\t".join(row.get(name, "") for name in header) for row in rows]
    return "\n".join(lines) + "\n"


def variant(tmp_path: Path, replace: dict[str, str] | None = None, append: str = ""):
    """The synthetic dataset file with text replaced or appended, loaded."""
    text = SPEC.read_text()
    for old, new in (replace or {}).items():
        assert text.count(old) == 1, old
        text = text.replace(old, new)
    path = tmp_path / "synthetic.toml"
    path.write_text(text + append)
    return load(path)


def problems_of(runs, text, spec) -> list[tuple[str, str, str]]:
    with pytest.raises(CurationError) as raised:
        curate(runs, text, spec)
    return raised.value.problems


def test_join_fills_curation(runs):
    out = curate(runs, sdrf_text(sdrf_rows(runs)), load(SPEC))
    assert out["run_accession"].to_list() == runs["run_accession"].to_list()
    assert out.select(pl.col("donor", "cell_type", "condition").is_null().any()).row(0) == (
        False, False, False,
    )  # fmt: skip
    assert set(out["cell_type"]) == {"type_a", "type_b"}
    assert dict(out["condition"].value_counts().rows()) == {"control": 7, "case": 7}


def test_unjoined_listed_together(runs):
    rows = sdrf_rows(runs)
    dropped = rows.pop(3)["Comment[ENA_RUN]"]
    rows.append({**rows[4], "Comment[ENA_RUN]": "XRR0", "Source Name": "x"})
    found = problems_of(runs, sdrf_text(rows), load(SPEC))
    joins = {(key, detail.split(";")[0]) for check, key, detail in found if check == "join"}
    assert (dropped, "run has no SDRF row (Comment[ENA_RUN])") in joins
    assert ("XRR0", "SDRF row matches no run_accession") in joins


def test_repeated_header_refused(runs):
    header = [*HEADER, "Characteristics[condition]"]
    found = problems_of(runs, sdrf_text(sdrf_rows(runs), header), load(SPEC))
    assert found == [
        ("sdrf header", "Characteristics[condition]", "appears 2 times; expected once")
    ]


def test_cross_check_fails_naming_row(runs):
    rows = sdrf_rows(runs)
    rows[5]["Source Name"] = "edited by hand"
    found = problems_of(runs, sdrf_text(rows), load(SPEC))
    assert [(c, k) for c, k, _ in found] == [("cross-check", runs["run_accession"][5])]


def test_normalise_then_map(runs):
    out = curate(runs, sdrf_text(sdrf_rows(runs)), load(SPEC))
    assert out["time_point"][1] == "fu2", "two spaces must collapse before the value map"


def test_replicate_order(runs):
    out = curate(runs, sdrf_text(sdrf_rows(runs)), load(SPEC))
    first, second = runs["run_accession"][0], runs["run_accession"][1]
    assert int(first[3:]) < int(second[3:]), "the case needs accession and time order to differ"
    # fu2 before fu10 by the declared order, though "fu10" < "fu2" as strings and by accession.
    assert out["replicate"][:2].to_list() == [2, 1]
    assert out["replicate"][2:].to_list() == [1] * (runs.height - 2)


def _lengths(runs) -> dict[str, float]:
    """Read length per run, written here independently of curation."""
    return {
        run: base / reads / (2 if layout == "PAIRED" else 1)
        for run, base, reads, layout in runs.select(
            "run_accession", pl.col("base_count").cast(pl.Int64),
            pl.col("read_count").cast(pl.Int64), "library_layout",
        ).rows()
    }  # fmt: skip


def test_read_length_rule(runs, tmp_path):
    lengths = _lengths(runs)
    fractional = sorted(r for r, v in lengths.items() if not v.is_integer())
    # The bins also hold each fractional length rounded, so rounding would bin it silently.
    whole = sorted({round(v) for v in lengths.values()})
    assert fractional and any(layout == "PAIRED" for layout in runs["library_layout"])
    spec = variant(tmp_path, append=f"\n[derived.batch]\nbins = {{ all = {whole} }}\n")
    found = problems_of(runs, sdrf_text(sdrf_rows(runs)), spec)
    assert sorted(k for c, k, _ in found if c == "batch") == fractional
    # Without bins nothing is rounded or refused, and the paired runs are halved.
    out = curate(runs, sdrf_text(sdrf_rows(runs)), load(SPEC))
    assert dict(out.select("run_accession", "read_length").rows()) == lengths


def test_expect_mismatch_fails(runs, tmp_path):
    spec = variant(tmp_path, {"control = 7\ncase = 7": "control = 6\ncase = 8"})
    found = problems_of(runs, sdrf_text(sdrf_rows(runs)), spec)
    assert {k for c, k, _ in found if c == "expect"} == {"condition = control", "condition = case"}


def test_schema_errors_listed_together(tmp_path):
    path = tmp_path / "broken.toml"
    path.write_text(
        SPEC.read_text()
        .replace('resource = "sdrf"', 'resource = "nowhere"')
        .replace("[columns.condition]", "[columns.state]")
        + "\n[bogus]\nx = 1\n"
    )
    with pytest.raises(DatasetError) as raised:
        load(path)
    message = str(raised.value)
    assert "unknown table [bogus]" in message
    assert "no [[resources]] entry named 'nowhere'" in message
    assert "must define ['condition']" in message


# --- through the CLI ------------------------------------------------------------------------


def _gather(root: Path, spec: Path = SPEC) -> subprocess.CompletedProcess[str]:
    env = {k: v for k, v in os.environ.items() if not k.startswith(("SLURM_", "CTX_"))}
    args = ["--root", str(root), "gather", "--dataset", "synthetic",
            "--fetchngs-dir", str(FIXTURE), "--dataset-file", str(spec)]  # fmt: skip
    return subprocess.run(
        [str(CTX), *args], cwd=root, env=env, capture_output=True, text=True, check=False
    )


def _place_sdrf(root: Path, rows) -> Path:
    path = root / load(SPEC).resource("sdrf").dest
    path.parent.mkdir(parents=True)
    path.write_text(sdrf_text(rows))
    return path


def test_cli_curated_end_to_end(runs, tmp_path):
    sdrf = _place_sdrf(tmp_path, sdrf_rows(runs))
    before = (hashlib.sha256(sdrf.read_bytes()).hexdigest(), sdrf.stat().st_mtime_ns)
    done = _gather(tmp_path)
    assert done.returncode == 0, done.stderr
    out = tmp_path / "outputs" / "01_gather"
    assert lint(out)[1] == []
    (sheet,) = out.glob("*_samplesheet_v1.parquet")
    frame = pl.read_parquet(sheet)
    assert frame.height == runs.height
    assert frame.select(pl.col("donor", "cell_type", "condition") == "").sum().row(0) == (0, 0, 0)
    assert {"time_point", "read_length"} <= set(frame.columns)
    for name in ("composition", "library_sizes"):
        (csv,) = out.glob(f"*_synthetic_{name}_v1.csv")
        written = sorted(p.name for p in out.glob(f"{csv.stem}*"))
        assert written == sorted(figure_files(csv.stem)), name
    (page,) = out.glob("*_synthetic_overview_v1.html")
    html = unpacked(page.read_text())
    assert "Curation from synthetic.toml" in html
    # The depth figure: a select per filter column, and ENCODE's two lines from the registry.
    assert all(f'<select data-column="{c}">' in html for c in DEPTH_FILTERS)
    good, usable = (reads / 1e6 for reads in _depth_standard()[0])
    for key, millions in (("good", good), ("usable", usable)):
        name = text.fill("gather", f"figures.library.{key}", millions=millions)
        assert f'"name":"{name}"' in html, name
    readme = next(out.glob("*_samplesheet_v1_README.md")).read_text()
    assert "sha256:synthetic.sdrf.txt" in readme
    assert (hashlib.sha256(sdrf.read_bytes()).hexdigest(), sdrf.stat().st_mtime_ns) == before


def test_failure_writes_no_sheet(runs, tmp_path):
    rows = sdrf_rows(runs)
    rows[2]["Characteristics[condition]"] = "Unknown"
    _place_sdrf(tmp_path, rows)
    done = _gather(tmp_path)
    assert done.returncode == 1
    assert "value map: condition: 'Unknown'" in done.stderr
    out = tmp_path / "outputs" / "01_gather"
    assert not list(out.glob("*_samplesheet_*"))
    (problems,) = out.glob("*_synthetic_curation_problems_v1.csv")
    assert pl.read_csv(problems)["check"].to_list() == ["value map", "expect", "expect"]


def test_slug_must_match(tmp_path):
    other = tmp_path / "other.toml"
    other.write_text(SPEC.read_text().replace('slug = "synthetic"', 'slug = "other"'))
    done = _gather(tmp_path, other)
    assert done.returncode == 1
    assert "is for 'other', not --dataset 'synthetic'" in done.stderr


def test_label_is_the_name_pages_show(tmp_path):
    # The slug names files and stays lowercase; the label is what a reader sees (D-014).
    assert load(SPEC).name == "synthetic"
    labelled = tmp_path / "labelled.toml"
    labelled.write_text(
        SPEC.read_text().replace('slug = "synthetic"', 'slug = "synthetic"\nlabel = "Step-9"')
    )
    spec = load(labelled)
    assert (spec.slug, spec.name) == ("synthetic", "Step-9")
    wrong = tmp_path / "wrong.toml"
    wrong.write_text(
        SPEC.read_text().replace('slug = "synthetic"', 'slug = "synthetic"\nlabel = 9')
    )
    with pytest.raises(DatasetError, match="label"):
        load(wrong)


@pytest.mark.parametrize("subset", [False, True])
def test_subset_curates_only_its_runs(runs, tmp_path, subset):
    # An SDRF row for a run the outdir lacks fails a full study and is dropped from a declared
    # subset, and said so in the README (D-042).
    rows = sdrf_rows(runs)
    rows.append({**rows[4], "Comment[ENA_RUN]": "XRR0", "Source Name": "x"})
    _place_sdrf(tmp_path, rows)
    spec = tmp_path / "synthetic.toml"
    declared = 'slug = "synthetic"\nsubset_of = "parent"' if subset else 'slug = "synthetic"'
    spec.write_text(SPEC.read_text().replace('slug = "synthetic"', declared))
    done = _gather(tmp_path, spec)
    out = tmp_path / "outputs" / "01_gather"
    if not subset:
        assert done.returncode == 1
        assert "XRR0: SDRF row matches no run_accession" in done.stderr
        return
    assert done.returncode == 0, done.stderr
    readme = next(out.glob("*_samplesheet_v1_README.md")).read_text()
    assert "14 of 15, the runs of this subset of parent" in readme


def test_subset_still_needs_a_row_per_run(runs, tmp_path):
    rows = sdrf_rows(runs)
    dropped = rows.pop(3)["Comment[ENA_RUN]"]
    _place_sdrf(tmp_path, rows)
    spec = tmp_path / "synthetic.toml"
    spec.write_text(
        SPEC.read_text().replace('slug = "synthetic"', 'slug = "synthetic"\nsubset_of = "parent"')
    )
    done = _gather(tmp_path, spec)
    assert done.returncode == 1
    assert f"{dropped}: run has no SDRF row" in done.stderr


def test_a_missing_registry_is_reported_not_raised(runs, tmp_path):
    # Gather reads the depth standards from the QC registry; a missing one is a clean exit 1.
    _place_sdrf(tmp_path, sdrf_rows(runs))
    env = {k: v for k, v in os.environ.items() if not k.startswith(("SLURM_", "CTX_"))}
    args = ["--root", str(tmp_path), "gather", "--dataset", "synthetic",
            "--fetchngs-dir", str(FIXTURE), "--dataset-file", str(SPEC)]  # fmt: skip
    done = subprocess.run([str(CTX), *args], cwd=tmp_path, capture_output=True, text=True,
                          env={**env, "CTX_QC_REGISTRY": str(tmp_path / "none.toml")},
                          check=False)  # fmt: skip
    assert done.returncode == 1
    assert "no QC registry at" in done.stderr and "Traceback" not in done.stderr
