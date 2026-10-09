"""The teacher guide reads each module's time and memory back from the pages, never typed."""

import tomllib
from pathlib import Path

import pytest

from ctx.io.naming import MODULES
from ctx.report import builder, prose, teacher, text
from ctx.report.builder import Group, ModulePage, ParameterCard, Reproduction
from ctx.report.hub import NOT_MODULES, TITLES

ROOT = Path(__file__).resolve().parents[1]


def _page(root: Path, module: str, day: str, dataset: str, cost: tuple[float, int, int]) -> Path:
    """A module page as the builder writes it, with known Reproduce values."""
    number = dict(MODULES)[module]
    folder = root / "outputs" / f"{number:02d}_{module}"
    folder.mkdir(parents=True, exist_ok=True)
    wall, own, child = cost
    page = ModulePage(
        title=f"{TITLES[module]}: {dataset}", **text.skeleton("guide"),
        parameters=(ParameterCard("--x", "x", "1", "x"),), results=(Group("g", ("t",)),),
        reproduce=Reproduction("ctx x", wall, own, child, "laptop, 4 threads", "none", {}),
    )  # fmt: skip
    return builder.write(page, folder / f"{day}_{module}_{dataset.lower()}_overview_v1.html")


@pytest.fixture
def root(tmp_path):
    _page(tmp_path, "gather", "2026-10-01", "Old", (999.0, 1, 1))
    _page(tmp_path, "gather", "2026-10-02", "Old", (5.25, 734003200, 5242880))
    _page(tmp_path, "qc", "2026-10-02", "Old", (146.4, 484442112, 1377828864))
    return tmp_path


def test_values_are_the_pages(root):
    """Each row holds its page's own Reproduce values, the newest page per dataset, in course
    order, and a module with no page reads not measured. Catches the process and child peaks
    swapped, an older page read, a time typed, and a module dropped or reordered."""
    (table,) = teacher.page([root], "x", 0.0).results[0].items
    rows = {r[0]: r for r in table.rows}
    shown = [m for m, _ in MODULES if m not in NOT_MODULES]
    assert [r[0] for r in table.rows] == [TITLES[m] for m in shown]
    assert rows["M1 Gather"][2:5] == ("5.2 s", "700 MiB", "5 MiB")
    assert rows["M2 Read QC"][2:5] == ("146.4 s", "462 MiB", "1,314 MiB")
    assert rows["M2 Read QC"][6].startswith("Track A")
    assert rows["M1 Gather"][6].startswith("Track B")
    none = text.word("teacher", "times.none")
    assert rows["M5 Differential expression"][1:6] == (none,) * 5


def test_written_and_run_from_the_command_line(root, tmp_path):
    """`python -m ctx.report.teacher OUT ROOT` writes the page, and a folder without outputs/
    is refused. Catches a root silently read as empty."""
    out = tmp_path / "guide" / "2026-10-02_report_teacher_v1.html"
    assert teacher.main([str(out), str(root)]) == out
    assert "146.4 s" in out.read_text()
    with pytest.raises(SystemExit):
        teacher.main([str(out), str(tmp_path / "guide")])


def test_no_number_typed_in_the_teacher_guide():
    """Every number on the page comes from a page it read; none sits in its text."""
    raw = tomllib.loads((ROOT / "src/ctx/report/content/teacher.toml").read_text())

    def strings(node):
        if isinstance(node, str):
            yield node
        elif isinstance(node, dict | list):
            for value in node.values() if isinstance(node, dict) else node:
                yield from strings(value)

    assert [s for s in strings(raw) if prose.quantities(s)] == []
