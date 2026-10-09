"""Every file each end-to-end run writes obeys the output contract (V3, V4).

One test over the session's module runs (conftest.py), so a module added later is held to the
same rules by adding its fixture to ``RUNS``, not by remembering to copy a per-module check.
Only the folders a run's own commands write are read: the planted inputs beside them (a
made-up gather or de folder) are test fixtures, not outputs.
"""

import html
import re
import shutil
from pathlib import Path

import pytest

from ctx.io.naming import NATIVE_DIRS, figure_files, lint

# Session fixture of conftest.py: the output folders its commands write.
RUNS = {
    "setup_run": ("00_setup",),
    "gather_run": ("01_gather",),
    "qc_run": ("02_qc",),
    "process_root": ("03_process",),
    "quant_run": ("04_quant",),
    "de_run": ("05_de", "06_function", "09_report"),
    "proteo_run": ("07_proteo",),
    "ncrna_run": ("08_ncrna",),
    "guide_run": ("09_report",),
    "genecard_run": ("11_genecard",),
    "clinical_run": ("12_clinical",),
}
TABLES = ("csv", "tsv", "parquet")
NO_TEST = re.compile(r"\bno test\b", re.IGNORECASE)
# A caption that tests names the test and states its p, FDR or adjustment (CONVENTIONS, Figures).
TEST_NAMES = re.compile(
    r"Mann-Whitney|Wilcoxon|Wald|log-rank|Kruskal-Wallis|Dunn|Welch|t-test|\bKS\b|Kolmogorov"
    r"|hypergeometric|GSEA|PERMANOVA|Cox|Fisher|chi-squared|Spearman|Pearson|permutation",
    re.IGNORECASE,
)
P_STATED = re.compile(r"\bp\b|\bFDR\b|\bBH\b|adjusted|\bq\b", re.IGNORECASE)
# Captions that break the rule on 2026-10-09 (CI run 37898263890), by figure without date and
# version. Their text is in content files other work packages own; each fix removes its line
# here, and the test fails until it does, so the list only shrinks.
KNOWN = {
    "de_desynth_calls_batch": "D, content/de.toml",
    "de_desynth_dispersion_a_case": "D, content/de.toml",
    "de_desynth_lfc_a_case_vs_b_case": "D, content/de.toml",
    "de_desynth_lfc_a_naive_vs_a_adjusted": "D, content/de.toml",
    "de_desynth_ma_a_case": "D, content/de.toml",
    "function_desynth_nes_up_set": "D, content/function.toml",
    "clinical_fake_brca_power": "C2, content/clinical.toml",
}


def caption_problem(caption: str) -> str | None:
    """Why a figure caption breaks the rule, or None: it names its test with a p, or says
    "no test"."""
    if NO_TEST.search(caption):
        return None
    if TEST_NAMES.search(caption) and P_STATED.search(caption):
        return None
    return "names no test with its p and does not say 'no test'"


def table_problems(folder: Path) -> list[str]:
    """Names the naming lint refuses, and tables with no README beside them."""
    scanned, problems = lint(folder)
    if not scanned:
        return [f"{folder.name}: no file"]
    for table in sorted(t for ext in TABLES for t in folder.rglob(f"*.{ext}")):
        if table.relative_to(folder).parts[0] in NATIVE_DIRS:
            continue  # a tool's own tables, kept as the tool wrote them
        if not table.with_name(f"{table.name.rsplit('.', 1)[0]}_README.md").is_file():
            problems.append(f"{table.name}: no README")
    return problems


def figure_problems(folder: Path) -> dict[str, str]:
    """Figure (stem without date and version) to what is wrong with it, for every figure
    under ``folder``, found by its white PDF."""
    problems = {}
    for white in sorted(folder.rglob("*_white.pdf")):
        stem = white.name.removesuffix("_white.pdf")
        name = re.sub(r"^\d{4}-\d{2}-\d{2}_|_v\d+$", "", stem)
        missing = [n for n in figure_files(stem) if not (white.parent / n).is_file()]
        if missing:
            problems[name] = f"missing {missing}"
            continue
        page = (white.parent / f"{stem}.html").read_text()
        captions = re.findall(r"<figcaption>(.*?)</figcaption>", page, re.DOTALL)
        text = html.unescape(re.sub(r"<[^>]+>", "", " ".join(captions)))
        why = caption_problem(text) if captions else "has no caption"
        if why:
            problems[name] = f"{why}: {text!r}"
    return problems


def _folders(request: pytest.FixtureRequest, run: str) -> list[Path]:
    root = Path(request.getfixturevalue(run)[0])
    return [root / "outputs" / name for name in RUNS[run]]


@pytest.mark.parametrize("run", RUNS)
def test_every_output_named_and_described(request, run):
    """V4: the naming lint clean in every folder a run writes, and a README beside every table."""
    for folder in _folders(request, run):
        assert table_problems(folder) == [], folder.name


@pytest.mark.parametrize("run", RUNS)
def test_every_figure_has_its_variants_and_a_tested_caption(request, run):
    """V3: every figure has its seven files, and its caption names its test or says no test;
    the captions in KNOWN, and only they, still break the rule."""
    folders = _folders(request, run)
    problems = {k: v for f in folders for k, v in figure_problems(f).items()}
    modules = {f.name.split("_", 1)[1] for f in folders}
    known = {k for k in KNOWN if k.split("_", 1)[0] in modules}
    assert {k: v for k, v in problems.items() if k not in KNOWN} == {}
    assert known - set(problems) == set(), "fixed: remove these from KNOWN"


def test_the_checks_see_a_planted_fault(tmp_path, quant_run):
    """The two checks on a copy of quant's folder: clean as written, then each fault planted.

    Mutations caught: a table with no README; a caption with neither a test nor "no test"; a
    figure's black SVG not written."""
    folder = tmp_path / "04_quant"
    shutil.copytree(quant_run[0] / "outputs" / "04_quant", folder)
    assert table_problems(folder) == []
    before = figure_problems(folder)
    assert set(before) <= set(KNOWN)

    (folder / "2026-10-09_quant_chrse_planted_v1.csv").write_text("a\n1\n")
    assert table_problems(folder) == ["2026-10-09_quant_chrse_planted_v1.csv: no README"]

    clean = sorted(p for p in folder.glob("*_white.pdf") if not any(k in p.name for k in before))
    page = clean[0].with_name(clean[0].name.replace("_white.pdf", ".html"))
    page.write_text(re.sub(r"<figcaption>.*?</figcaption>", "<figcaption>Reads.</figcaption>",
                           page.read_text(), flags=re.DOTALL))  # fmt: skip
    clean[1].with_name(clean[1].name.replace("_white.pdf", "_black.svg")).unlink()
    after = figure_problems(folder)
    assert len(after) == len(before) + 2
    assert any(v.startswith("names no test") for k, v in after.items() if k not in before)
    assert any(v.startswith("missing") for k, v in after.items() if k not in before)


@pytest.mark.parametrize(
    ("caption", "passes"),
    [
        ("Counts of 4 samples; no test.", True),
        ("Mann-Whitney U, two-sided; p is BH-adjusted over 3 features; n = 4, 4.", True),
        ("Spearman rho over 12 genes.", False),
        ("Reads per sample, 4 samples.", False),
    ],
)
def test_caption_rule(caption, passes):
    """The rule itself: a test needs its p; a descriptive caption must say so."""
    assert (caption_problem(caption) is None) is passes
