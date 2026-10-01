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


def test_a_caption_has_no_backticks():
    page = {**CLEAN, "figures": {"f": {"caption": "Made up; `fastp` kept {n} reads."}}}
    assert any("a backtick" in p for p in prose.check_page("page", page))


def test_page_files_pass():
    """Every content file obeys the rules: CONVENTIONS says a violation fails CI, and this is it."""
    assert prose.lint() == []


def test_course_text_is_checked():
    """The course page and kit README have no skeleton, so every string of them is checked."""
    found = list(prose.check_all("course", {"intro": "We dive in with 12 genes."}))
    assert any("banned word 'We'" in f for f in found) and any("typed number" in f for f in found)


@pytest.mark.parametrize(
    "verb",
    ["know", "knows", "knew", "known", "knowing", "bring", "brings", "brought", "bringing",
     "want", "wants", "wanted", "wanting", "tell", "tells", "told", "telling",
     "arrive", "arrives", "arrived", "arriving"],
)  # fmt: skip
def test_figurative_verbs_in_every_form(verb):
    """Software does not know, bring, want, tell or arrive, in any inflection (D-048).

    Mutation: any one form removed from BANNED fails its row; before D-048, 'knew', 'known',
    'brought', 'want', 'wanted', 'tell', 'told' and 'arrived' all passed, and before review
    every -ing form did.
    """
    found = list(prose.check_page("page", {**CLEAN, "method": f"The tool {verb} the reads."}))
    assert any(f"banned word '{verb}'" in f for f in found), found


def test_only_the_questions_ask():
    """A question mark in a title, lead or section fails; the Questions section may ask.

    Mutations: check_asks skipping titles (a key named title) passes the planted title;
    asks() flagging only a text that ends in '?' passes the lead's question mid-text;
    check_asks skipping lists passes the question in a pitfall.
    """
    page = {
        **CLEAN,
        "groups": {"title": "Where do reads go?", "lead": "Which exons? Count them."},
        "pitfalls": ["Trimming: why trim at all? Adapters do not map."],
        "questions": [{"prompt": "Why is it so?", "answer": "Because."}],
    }
    found = list(prose.check_asks("page", page))
    assert found == [
        "page.groups.title: a question; state it, or ask it under Questions",
        "page.groups.lead: a question; state it, or ask it under Questions",
        "page.pitfalls[0]: a question; state it, or ask it under Questions",
    ]
    assert list(prose.check_asks("page", {"method": "Run `ls ?.fq` over {n?} files."})) == []


def test_a_mark_counts_as_its_display_text():
    """Budgets and banned words read a glossary mark as the words the reader sees.

    Mutation: check_text without plain(): '[[Phred quality score|Q]]' counts as three words,
    and the key's 'robust' is read as a banned word the reader never sees.
    """
    rest = " ".join(["word"] * 9) + ". " + ". ".join([" ".join(["Word"] * 10)] * 2) + "."
    mark = "[[Phred quality score|Q]]"
    assert list(prose.check_page("page", {**CLEAN, "purpose": f"{mark} {rest}"})) == []
    over = list(prose.check_page("page", {**CLEAN, "purpose": f"{mark} Q {rest}"}))
    assert any("31 words; at most 30" in f for f in over), over
    shown = {**CLEAN, "method": "Flag by [[robust z-score|MAD]]; read [[md5|very]] slowly."}
    assert [f.split(": ")[-1] for f in prose.check_page("page", shown)] == ["banned word 'very'"]


def test_a_strand_end_is_not_a_number():
    """5' and 3' name a strand end, not a quantity; a trailing stop ends the number.

    Mutation: the prime taken out of the number pattern's lookahead fails the first line.
    """
    assert prose.quantities("Reads pile up at the 3' end, not the 5' end.") == []
    assert prose.quantities("Check 210 runs.") == ["210"]


def test_a_sentence_may_start_in_lower_case():
    """A stop before a lower-case word ends a sentence (fastp, lncRNAs, nf-core, pytximport
    start sentences), unless it closes an abbreviation such as e.g., i.e., vs. or et al.

    Mutations: the splitter back to capitals only (the two sentences read as one of 30 words);
    the abbreviation guard removed (the 26 words around "e.g." read as two sentences).
    """
    two = "Built from the GTF with its transcripts. " + " ".join(["fastp"] * 23) + "."
    assert list(prose.check_page("page", {**CLEAN, "method": two})) == []
    four = "One. two. three. four."
    assert any("a paragraph of 4 sentences" in f
               for f in prose.check_page("page", {**CLEAN, "method": four}))  # fmt: skip
    abbreviated = " ".join(["word"] * 12) + " e.g. " + " ".join(["fastp"] * 12) + " et al. 2016."
    found = list(prose.check_page("page", {**CLEAN, "method": abbreviated}))
    assert any("a sentence of 28 words" in f for f in found), found
