"""U.3: the cost of a pipeline run, read from Nextflow's trace by its own column names.

On the committed fixtures: the rnaseq outdir's trace holds 41 rows, 40 of them CACHED (the
fixture was relaunched with -resume), and the fetchngs download's 20, all COMPLETED. Each test
names the mutation it catches.
"""

from pathlib import Path

import pytest
from conftest import ctx_run

from ctx.gather import trace
from ctx.report import text

ROOT = Path(__file__).resolve().parents[1]
RNASEQ = ROOT / "tests" / "fixtures" / "rnaseq"
DOWNLOAD = ROOT / "tests" / "fixtures" / "fetchngs" / "download"
HEADER = "task_id\thash\tnative_id\tname\tstatus\texit\ttag\tduration\trealtime\tpeak_rss\n"


@pytest.mark.parametrize(
    ("value", "expected"),
    [("3 MB", 3 * 1024**2), ("1.1 GB", round(1.1 * 1024**3)), ("334.4 KB", round(334.4 * 1024)),
     ("0", 0), ("-", None), ("", None)],
)  # fmt: skip
def test_memory_is_1024_based(value, expected):
    """Nextflow's MemoryUnit prints powers of 1024. Catches GB read as MB, and powers of 1000."""
    assert trace.memory(value) == expected


@pytest.mark.parametrize(
    ("value", "expected"),
    [("19ms", 0.019), ("7.9s", 7.9), ("1m 19s", 79.0), ("1h 2m 3s", 3723.0), ("1d 1h", 90000.0),
     ("-", None)],
)  # fmt: skip
def test_seconds(value, expected):
    """Catches minutes read as milliseconds (the m of ms), and a part dropped."""
    found = trace.seconds(value)
    assert found is None if expected is None else found == pytest.approx(expected)


def test_unreadable_duration():
    """A duration with text it cannot account for is refused, never read in part."""
    with pytest.raises(trace.TraceError, match="unreadable duration"):
        trace.seconds("3 fortnights")


def test_cached_rows_cost_nothing():
    """Only the task the last launch ran is costed. Catches CACHED rows counted as cost: the
    fixture's 40 cached rows would add 19 processes."""
    frame = trace.read(RNASEQ)
    assert frame.height == 41 and frame.filter(frame["status"] == "CACHED").height == 40
    costs = trace.costs(frame)
    assert costs.rows() == [("MULTIQC", 1, 1, trace.memory("1.1 GB"), 7.9)]


def test_sample_from_the_name(tmp_path):
    """The sample is the text in parentheses at the end of ``name``. Catches the sample read
    from ``tag``, which a pipeline may set to anything, and the process kept with its path."""
    info = tmp_path / "pipeline_info"
    info.mkdir()
    rows = [
        "1\ta\t1\tNFCORE_RNASEQ:RNASEQ:SALMON_QUANT (s1)\tCOMPLETED\t0\tx\t2s\t1m 1s\t2 GB\n",
        "2\tb\t2\tNFCORE_RNASEQ:RNASEQ:SALMON_QUANT (s2)\tCOMPLETED\t0\tx\t2s\t1s\t3 GB\n",
        "3\tc\t3\tNFCORE_RNASEQ:RNASEQ:SALMON_QUANT (s2)\tFAILED\t1\tx\t2s\t1s\t1 GB\n",
        "4\td\t4\tNFCORE_RNASEQ:CUSTOM_DUMPSOFTWAREVERSIONS\tCOMPLETED\t0\t-\t1s\t-\t-\n",
    ]
    (info / "execution_trace_2026-10-08_00-00-00.txt").write_text(HEADER + "".join(rows))
    costs = {r[0]: r[1:] for r in trace.costs(trace.read(tmp_path)).rows()}
    assert costs["SALMON_QUANT"] == (3, 2, 3 * 1024**3, 63.0)
    assert costs["CUSTOM_DUMPSOFTWAREVERSIONS"] == (1, 0, None, 0.0)


def test_launches_add(tmp_path):
    """Two launches' traces: each task once, the relaunch's cached rows left out. Catches only
    the newest trace read."""
    info = tmp_path / "pipeline_info"
    info.mkdir()
    first = "1\ta\t1\tP:INDEX (genome)\tCOMPLETED\t0\tx\t1s\t10m\t9 GB\n"
    again = first.replace("COMPLETED", "CACHED")
    again += "2\tb\t2\tP:QUANT (s1)\tCOMPLETED\t0\tx\t1s\t1m\t1 GB\n"
    (info / "execution_trace_2026-10-08_00-00-00.txt").write_text(HEADER + first)
    (info / "execution_trace_2026-10-08_01-00-00.txt").write_text(HEADER + again)
    assert sorted(trace.costs(trace.read(tmp_path))["process"]) == ["INDEX", "QUANT"]


def test_trace_without_columns(tmp_path):
    """A trace with other columns names what it lacks."""
    (tmp_path / "pipeline_info").mkdir()
    (tmp_path / "pipeline_info" / "execution_trace_x.txt").write_text("task_id\tname\n1\tA\n")
    with pytest.raises(trace.TraceError, match="peak_rss"):
        trace.read(tmp_path)


def test_gather_page_shows_the_cost(tmp_path):
    """ctx gather puts the table on its page, a row per process that ran. Catches the table
    left off the page or put on it twice, and the cached-row count left out of its note."""
    group = trace.table((("nf-core/fetchngs", DOWNLOAD), ("nf-core/rnaseq", RNASEQ)))
    assert group is not None
    note, table = group.items
    assert note == text.fill("gather", "trace.note", cached=40, total=61)
    assert [r[:3] for r in table.rows][-1] == ("nf-core/rnaseq", "MULTIQC", "1")
    done = ctx_run("--root", str(tmp_path), "gather", "--dataset", "nfcore_test",
                   "--fetchngs-dir", str(DOWNLOAD), cwd=tmp_path)  # fmt: skip
    assert done.returncode == 0, done.stderr
    (page,) = (tmp_path / "outputs" / "01_gather").glob("*_overview_v1.html")
    html = page.read_text()
    # Once: the first build appended it twice, found by planting the table's removal.
    assert html.count(text.word("gather", "trace.title")) == 1 and "SRA_IDS_TO_RUNINFO" in html
