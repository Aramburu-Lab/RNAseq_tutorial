"""The prose lint: each rule catches its planted violation, and the page files pass (D-045)."""

import pytest

from ctx.report import prose

CLEAN = {"purpose": "Verify the download.", "method": "`ctx gather` reads `samplesheet.csv`."}


@pytest.mark.parametrize(
    ("section", "text", "found"),
    [
        ("purpose", " ".join(["word"] * 31), "31 words; at most 30"),
        ("method", " ".join(["word"] * 26) + ".", "a sentence of 26 words"),
        ("method", "One. Two. Three. Four.", "a paragraph of 4 sentences"),
        ("purpose", "Here we check reads.", "banned word 'we'"),
        ("interpretation", "The download brings files.", "banned word 'brings'"),
        ("purpose", "Check 210 runs.", "a typed number"),
    ],
)
def test_each_rule_catches_its_violation(section, text, found):
    problems = list(prose.check_page("page", {**CLEAN, section: text}))
    assert any(found in p for p in problems), problems


def test_names_with_digits_are_not_numbers():
    text = "Step-1 runs on GRCh38 chr22; M5 checks each SHA-256."
    assert list(prose.check_page("page", {**CLEAN, "purpose": text})) == []


def test_placeholders_and_code_are_not_prose():
    text = "{n_ok} of {n_files} files passed `md5sum -c` in `fastq/md5/`."
    assert list(prose.check_page("page", {**CLEAN, "purpose": text})) == []


def test_lists_and_cards_have_budgets():
    items = {**CLEAN, "pitfalls": ["One.", "Two.", "Three.", "Four."]}
    assert any("4 items" in p for p in prose.check_page("page", items))
    cards = {"options": {"--x": {"meaning": " ".join(["word"] * 13), "effect": "Moves it."}}}
    assert any("13 words; at most 12" in p for p in prose.check_cards(cards))
