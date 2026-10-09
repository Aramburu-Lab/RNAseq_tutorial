"""Every file of every end-to-end run obeys the output contract (V3, V4).

One test over the session's module runs (conftest.py), so a module added later is held to the
same rules by adding its fixture to ``RUNS``, not by remembering to copy a per-module check.
"""

import html
import re
from pathlib import Path

import pytest

from ctx.io.naming import figure_files, lint

# The session fixtures of conftest.py; the pipeline root holds de, function and report.
RUNS = ["setup_run", "gather_run", "qc_run", "process_root", "quant_run", "de_run",
        "clinical_run", "proteo_run", "genecard_run", "ncrna_run", "guide_run"]  # fmt: skip
TABLES = ("csv", "tsv", "parquet")
NO_TEST = re.compile(r"\bno test\b", re.IGNORECASE)
# A caption that tests names the test and states its p, FDR or adjustment (CONVENTIONS, Figures).
TEST_NAMES = re.compile(
    r"Mann-Whitney|Wilcoxon|Wald|log-rank|Kruskal-Wallis|Dunn|Welch|t-test|\bKS\b|Kolmogorov"
    r"|hypergeometric|GSEA|PERMANOVA|Cox|Fisher|chi-squared|Spearman|Pearson|permutation",
    re.IGNORECASE,
)
P_STATED = re.compile(r"\bp\b|\bFDR\b|\bBH\b|adjusted|\bq\b", re.IGNORECASE)


def caption_problem(caption: str) -> str | None:
    """Why a figure caption breaks the rule, or None: it names its test with a p, or says
    "no test"."""
    if NO_TEST.search(caption):
        return None
    if TEST_NAMES.search(caption) and P_STATED.search(caption):
        return None
    return "names no test with its p and does not say 'no test'"


def _root(request: pytest.FixtureRequest, run: str) -> Path:
    found = request.getfixturevalue(run)
    return Path(found[0])


def _figures(folder: Path) -> list[tuple[Path, str]]:
    """(folder, stem) of every figure under ``folder``, found by its white PDF."""
    found = folder.rglob("*_white.pdf")
    return sorted((p.parent, p.name.removesuffix("_white.pdf")) for p in found)


@pytest.mark.parametrize("run", RUNS)
def test_every_output_named_and_described(request, run):
    """V4: the naming lint clean in every output folder, and a README beside every table.

    Mutation caught: a table written without its README (one planted beside a run's outputs)."""
    folders = sorted(p for p in (_root(request, run) / "outputs").iterdir() if p.is_dir())
    assert folders
    for folder in folders:
        scanned, violations = lint(folder)
        assert scanned > 0 and violations == [], (folder.name, violations)
        for table in (t for ext in TABLES for t in folder.glob(f"*.{ext}")):
            readme = table.with_name(f"{table.name.rsplit('.', 1)[0]}_README.md")
            assert readme.is_file(), f"{folder.name}/{table.name} has no README"


@pytest.mark.parametrize("run", RUNS)
def test_every_figure_has_its_variants_and_a_tested_caption(request, run):
    """V3: every figure has its seven files, and its caption names its test or says no test.

    Mutation caught: a caption with neither (one planted in a figure's page); a figure's black
    SVG not written."""
    root = _root(request, run)
    problems = []
    for module in (p for p in (root / "outputs").iterdir() if p.is_dir()):
        for folder, stem in _figures(module):
            missing = [n for n in figure_files(stem) if not (folder / n).is_file()]
            if missing:
                problems.append(f"{stem}: missing {missing}")
                continue
            page = (folder / f"{stem}.html").read_text()
            captions = re.findall(r"<figcaption>(.*?)</figcaption>", page, re.DOTALL)
            text = html.unescape(re.sub(r"<[^>]+>", "", " ".join(captions)))
            why = caption_problem(text) if captions else "has no caption"
            if why:
                problems.append(f"{stem}: {why}: {text[:200]!r}")
    assert problems == []


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
