"""The glossary (D-048): its file and lint, the marks in page text, the page, and its links.

Each test names the mutation it catches; the list was run once with each mutation applied.
"""

import copy
import dataclasses
import json
import re
import subprocess
import sys
from datetime import date
from pathlib import Path

import plotly.graph_objects as go
import pytest
from pypdf import PdfReader
from synth import counts, demo_bar
from typer.testing import CliRunner

from ctx.io.naming import check, stem
from ctx.report import course, glossary, glossary_lint, marks
from ctx.report.builder import (
    FigureBlock,
    FigureText,
    GlossaryPage,
    HubPage,
    ModulePage,
    ParameterCard,
    Reproduction,
    SelfCheck,
    render,
    unpacked,
    write,
)
from ctx.report.cap import MAX_BYTES
from ctx.viz.save import write_figure

ROOT = Path(__file__).resolve().parents[1]
BUILD_KIT = ROOT / "scripts" / "python" / "build_kit.py"
RAW = glossary.read()
TODAY = date.today().isoformat()


def _entry(**fields):
    """A valid entry with the fields given replaced; a field given as ... is left out."""
    base = {"id": "made_up", "term": "made-up term", "aliases": ["planted"], "category": "data",
            "short": "A term the tests plant.", "long": "Planted by a test.",
            "see_also": ["read"]}  # fmt: skip
    return {k: v for k, v in {**base, **fields}.items() if v is not ...}


def _lint(release=False, **fields):
    """Errors and warnings of the packaged glossary with one planted entry."""
    return glossary_lint.lint([*copy.deepcopy(RAW), _entry(**fields)], release=release)


def test_the_packaged_glossary_passes_its_lint():
    """189 entries, R2's NCBI Gene ID the latest, after M10's cumulative incidence, M7's
    proteomic ruler and M8's five (context++ score to histone mRNA), and no error; the draft's
    8 unverified entries carry sources now (L.5).

    Mutation: 'robust' written back into mad.short fails it (a banned word).
    """
    errors, warnings = glossary_lint.lint()
    assert errors == []
    assert len(glossary.load()) == len(RAW) == 189
    assert not [w for w in warnings if "confirm: true" in w]


@pytest.mark.parametrize(
    ("fields", "found"),
    [
        ({"colour": "red"}, "unknown field 'colour'"),
        ({"id": "Made-Up"}, "id is not snake_case"),
        ({"category": "misc"}, "category 'misc' is not declared"),
        ({"id": "read"}, "id 'read' is used 2 times"),
        ({"aliases": ["Salmon"]}, "'salmon' also names salmon"),
        ({"aliases": ["Phred quality score"]}, "'phred quality score' also names phred"),
        ({"aliases": [3]}, "aliases holds a value that is not text"),
        ({"aliases": ...}, "no 'aliases'"),
        ({"aliases": None}, "aliases is not list"),
        ({"short": "  "}, "short is blank"),
        ({"term": ""}, "term is blank"),
        ({"see_also": ["nowhere"]}, "see_also 'nowhere' is no other entry"),
        ({"see_also": ["made_up"]}, "see_also 'made_up' is no other entry"),
        ({"short": " ".join(["word"] * 21)}, "21 words; at most 20"),
        ({"long": " ".join(["word"] * 61)}, "61 words; at most 60"),
        ({"long": " ".join(["word"] * 26) + "."}, "a sentence of 26 words"),
        ({"long": "One. Two. Three. Four."}, "a paragraph of 4 sentences"),
        ({"long": "Salmon is very fast."}, "banned word 'very'"),
        ({"rna_vs_dna": "The aligner knew the reads."}, "banned word 'knew'"),
        ({"typical": "Is it good?"}, "a question"),
    ],
)
def test_each_rule_catches_its_planted_entry(fields, found):
    """Each schema and prose rule fails its planted entry.

    Mutations: each rule's check removed in turn (the alias-against-term check, the id count,
    the word budget of short, the question rule on definitions) fails its row; the term left
    out of the names compared (an alias equal to phred's term, not its id, passes); aliases
    optional again (the missing and the null aliases pass); blank text allowed (the blank
    short and term pass); an entry allowed to see_also itself.
    """
    errors, _ = _lint(**fields)
    assert any(found in e for e in errors), errors


def test_numbers_and_confirm_warn_now_and_fail_at_a_release():
    """A number with no source, and confirm: true, are warnings, and errors at a release tag.

    A definition may type its number (the typed-number rule is off), once a source names it.
    Mutations: release ignored in lint() (the release rows pass); confirm left a warning at a
    release tag; typical left out of the fields read for numbers.
    """
    unsourced = {"short": "Q30 means 1 error in 1,000.", "confirm": True}
    confirm = "made_up: confirm: true; verify it and give its source"
    errors, warnings = _lint(**unsourced)
    assert errors == []
    assert {"made_up.short: the number 1 has no source", confirm} <= set(warnings)
    errors, warnings = _lint(release=True, **unsourced)
    assert {"made_up.short: the number 1 has no source", confirm} <= set(errors)
    assert confirm not in warnings
    _, warnings = _lint(typical="About 5 M reads per sample.")
    assert "made_up.typical: the number 5 has no source" in warnings
    sourced = _lint(release=True, short="Q30 means 1 error in 1,000.", source="Ewing 1998")
    assert not [e for e in sourced[0] if e.startswith("made_up")]


def _git(repo: Path, *args: str) -> None:
    config = ["-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false"]
    subprocess.run(["git", "-C", str(repo), *config, *args], check=True, capture_output=True)


def test_release_tag_reads_github_or_git(monkeypatch, tmp_path):
    """At a v* tag, GitHub's ref names it; any other tag is no release; without GitHub's ref,
    git describe finds a v* tag on HEAD, and none on an untagged commit.

    Mutations: the GITHUB_REF_TYPE branch removed (the tag is read from git: None here); any
    GitHub tag taken as a release ("nightly" returned); the git describe path dropped (the
    tagged repo gives None).
    """
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "-q")
    _git(repo, "commit", "-q", "--allow-empty", "-m", "one")
    monkeypatch.setattr(glossary_lint, "REPO", repo)
    monkeypatch.setenv("GITHUB_REF_TYPE", "tag")
    monkeypatch.setenv("GITHUB_REF_NAME", "v1.0.0")
    assert glossary_lint.release_tag() == "v1.0.0"
    monkeypatch.setenv("GITHUB_REF_NAME", "nightly")
    assert glossary_lint.release_tag() is None
    monkeypatch.delenv("GITHUB_REF_TYPE")
    monkeypatch.delenv("GITHUB_REF_NAME")
    _git(repo, "tag", "v9.9.9")
    assert glossary_lint.release_tag() == "v9.9.9"


def _content(tmp_path: Path, **page) -> Path:
    """A content folder with one page file."""
    folder = tmp_path / "content"
    folder.mkdir(parents=True)
    lines = []
    for key, value in page.items():
        if isinstance(value, dict):
            lines.append(f"[{key}]")
            lines += [f"{k} = {json.dumps(v)}" for k, v in value.items()]
        else:
            lines.insert(0, f"{key} = {json.dumps(value)}")
    (folder / "page.toml").write_text("\n".join(lines) + "\n")
    return folder


def test_marks_in_page_files_resolve_and_sit_where_they_render(tmp_path):
    """An unknown term, and a mark in text the page prints raw, fail; unmarked terms warn.

    Mutations: the lookup of each mark removed (the unknown term passes); section() returning
    "purpose" for every key (the lead's mark passes); the seen-set check dropped from
    _first_mentions (the marked FASTQ is also warned); the reading order ignored, marks taken
    from the whole section first (fastp's bare first mention passes).
    """
    folder = _content(
        tmp_path,
        purpose="Check [[no such term]] and [[FASTQ]] files; FASTQ again.",
        method="Count reads with featureCounts.",
        groups={"lead": "A lead with [[md5]]."},
    )
    errors, warnings = glossary_lint.page_marks(folder)
    assert errors == [
        "page: purpose: 'no such term' is not in the glossary",
        "page: groups.lead: a mark where the page prints it raw",
    ]
    folder = _content(tmp_path / "ok", purpose="Check [[FASTQ]] files; FASTQ again.",
                      method="Count reads with featureCounts.",
                      input="Trim with fastp, then read the [[fastp]] report.")  # fmt: skip
    errors, warnings = glossary_lint.page_marks(folder)
    assert errors == []
    assert "page: method: 'featureCounts' is featurecounts, unmarked" in warnings
    assert "page: input: 'fastp' is fastp, marked only at a later mention" in warnings
    assert not [w for w in warnings if "fastq" in w]


def test_a_variants_sections_render_marks_as_their_own_section(tmp_path):
    """M7's tumour-against-normal text replaces the page's sections (text.skeleton), so a mark
    there renders, and a term the page's own section marked is still unmarked in the variant.

    Mutations: section() without the variants branch (the variant's mark is an error); the
    variant grouped with the page's own section (its bare FASTQ is not warned)."""
    variant = {"variants.design": {"interpretation": "Read FASTQ, then [[md5]]."}}
    folder = _content(tmp_path, purpose="Check [[FASTQ]] files.", interpretation="[[FASTQ]] read.",
                      **variant)  # fmt: skip
    errors, warnings = glossary_lint.page_marks(folder)
    assert errors == []
    assert "page: design: interpretation: 'FASTQ' is fastq, unmarked" in warnings


def test_every_page_file_marks_only_terms_that_exist():
    """The shipped page files: every mark resolves and sits where the page renders it.

    Mutation: a mark of an unknown term planted in gather.toml's method fails it.
    """
    errors, _ = glossary_lint.page_marks()
    assert errors == []


def test_the_m3_page_marks_each_term_where_a_section_first_names_it():
    """D-048 on process.toml: in every section the first mention of a glossary term carries
    its mark. "read" and "counts", the subject of nearly every sentence, are left bare, as on
    every other module's page; any other bare first mention is a warning this test refuses.

    Mutation caught: Interpretation's [[Salmon]] written bare again (the M3 review found it,
    with Method's Spearman and the Results captions' featureCounts, unstranded and reference).
    """
    _, warnings = glossary_lint.page_marks()
    process = [w for w in warnings if w.startswith("process: ")]
    bare = ("is read, unmarked", "is raw_counts, unmarked")
    # The two bare words are found, so the lint read the page file at all.
    assert process
    assert [w for w in process if not w.endswith(bare)] == []


def test_setup_writes_the_operating_system_out_and_os_still_names_survival():
    """L.6, D-104: setup's page says "operating system", so the lint finds no OS there with no
    exception; written as survival's alias, "OS" in any section still names the endpoint.

    Mutation caught: setup.toml's Input with "the OS reports" again (setup's 'OS' warned)."""
    _, warnings = glossary_lint.page_marks()
    assert not [w for w in warnings if w.startswith("setup: ") and "'OS'" in w]
    text = ["The memory the OS reports."]
    assert list(glossary_lint._first_mentions(text, glossary_lint._names())) == [
        "'OS' is survival_endpoints, unmarked"]  # fmt: skip


def test_lookup_by_id_term_or_alias_in_any_case(capsys):
    """One entry answers to its id, its term and each alias, a long name in any case (D-104);
    explain prints it.

    Mutation: lookup without .lower(): 'Q SCORE' raises.
    """
    assert {glossary.lookup(n).id for n in ("phred", "Phred quality score", "Q SCORE")} == {"phred"}
    with pytest.raises(glossary.GlossaryError, match="not in the glossary"):
        glossary.lookup("no such term")
    glossary.explain("q score")
    printed = capsys.readouterr().out
    assert printed.startswith("Phred quality score (Read quality)\nProbability that a base")
    assert "Source: Ewing & Green 1998" in printed and "See also: Q30," in printed


def test_a_short_or_capitalised_name_resolves_only_as_written():
    """L.6, D-104: a mark resolves "OS", "GO" and "rho" only as written, the rule the lint's
    bare-word search uses, so [[os]] is no mark of overall survival; longer names keep any case.

    Mutation caught: the any-case index holding every name again ([[os]] resolves)."""
    assert glossary.lookup("OS").id == "survival_endpoints"
    assert glossary.lookup("Overall Survival").id == "survival_endpoints"
    for written in ("os", "Os", "go", "Rho"):
        with pytest.raises(glossary.GlossaryError, match="not in the glossary"):
            glossary.lookup(written)


def test_plain_prints_the_display_text():
    """Outside a page a mark is its display text; a ! and a display both parse.

    Mutation: plain() returning the key instead of the display text fails the Q scores case.
    """
    text = "Read [[phred|Q scores]] from [[FASTQ]] files, then [[md5!]]."
    assert marks.plain(text) == "Read Q scores from FASTQ files, then md5."
    assert marks.marks(text)[2] == marks.Mark("md5", "md5", True)


def _module(**overrides) -> ModulePage:
    fields = {
        "title": "M1 Gather: Demo",
        "purpose": "Check [[FASTQ]] files, [[fastq]] again, and [[FASTQ!]] once more.",
        "input": "The [[md5]] files.",
        "method": "A [[FASTQ]] in `[[FASTQ]]` code.",
        "parameters": (ParameterCard("--root", "Where outputs go.", "cwd", "Moves outputs."),),
        "results": ("One result.",),
        "interpretation": "In practice.",
        "questions": (SelfCheck("Why?", "Because one."),),
        "pitfalls": ("A mistake.",),
        "reproduce": Reproduction("ctx gather", 1.0, 1, 0, "laptop", "B", {"python": "3"}),
    }
    return ModulePage(**{**fields, **overrides})


def _section(html: str, name: str) -> str:
    return re.search(rf'<section id="{name}">(.*?)</section>', html, re.S).group(1)


def test_first_mention_per_section_carries_the_icon():
    """Per section, the first mention gets the icon, ! forces one, and code stays literal.

    Mutations: the shown-set check removed (every mention gets an icon); a section key
    ignored (method's first mention loses its icon); the force flag ignored.
    """
    html = render(_module())
    purpose, method = _section(html, "purpose"), _section(html, "method")
    assert purpose.count('data-term="fastq"') == 2
    assert method.count('data-term="fastq"') == 1 and "<code>[[FASTQ]]</code>" in method
    button = re.search(r'<button class="info"[^>]*data-term="fastq"[^>]*>', purpose).group(0)
    assert 'type="button"' in button and 'aria-expanded="false"' in button
    assert 'aria-label="Explain: FASTQ file"' in button
    (payload,) = re.findall(r'id="ctx-glossary">(.*?)</script>', html, re.S)
    assert set(json.loads(payload)["entries"]) == {"fastq", "md5"}


def test_a_page_without_marks_carries_no_tooltips():
    """A page that marks nothing embeds no glossary data and links no glossary.

    Mutation: the tooltips payload written whatever the page marks fails it.
    """
    html = render(_module(purpose="Plain.", input="Plain.", method="Plain."))
    assert 'type="application/json" id="ctx-glossary"' not in html
    assert 'class="glossary-link"' not in html


def _root(tmp_path: Path) -> Path:
    out = tmp_path / "root" / "outputs"
    (out / "01_gather").mkdir(parents=True)
    return out


def test_writing_a_page_writes_the_glossary_it_links(tmp_path):
    """A module page links ../09_report/<today>_report_glossary_v1.html, which then exists,
    is named by the naming rule and anchors every term the page marks.

    Mutations: page_path to the page's own folder (the relative link changes); the glossary
    write removed (the link does not resolve).
    """
    out = _root(tmp_path)
    page = write(_module(), out / "01_gather" / f"{TODAY}_gather_demo_overview_v1.html")
    html = unpacked(page.read_text())
    (link,) = re.findall(r'class="glossary-link" href="([^"]+)"', html)
    assert link == f"../09_report/{TODAY}_report_glossary_v1.html"
    target = (page.parent / link).resolve()
    assert target.is_file() and check(target.name)
    (payload,) = re.findall(r'id="ctx-glossary">(.*?)</script>', html, re.S)
    assert json.loads(payload)["href"] == link
    anchors = set(re.findall(r'<article class="entry" id="([^"]+)"', unpacked(target.read_text())))
    assert {"fastq", "md5"} <= anchors


def test_a_page_outside_a_ctx_root_keeps_its_tooltips_and_links_no_glossary(tmp_path):
    """A page written outside outputs/NN_<module>/, as ``ctx figure`` will write into the
    reader's own folder (Phase 2D), has no glossary beside it: its tooltips open with no link.

    Mutation: page_path without its folder check: the page links ../09_report/ and writes a
    glossary beside tmp_path.
    """
    (tmp_path / "work").mkdir()
    page = write(_module(), tmp_path / "work" / "page.html")
    html = unpacked(page.read_text())
    (payload,) = re.findall(r'id="ctx-glossary">(.*?)</script>', html, re.S)
    assert json.loads(payload)["href"] == "" and set(json.loads(payload)["entries"]) == {
        "fastq", "md5"}  # fmt: skip
    assert 'class="glossary-link"' not in html
    assert sorted(p.name for p in tmp_path.rglob("*") if p.is_file()) == ["page.html"]


def test_the_hub_links_the_glossary(tmp_path):
    """The hub marks no term but always links the glossary beside it.

    Mutation: the HubPage condition dropped from write(): the hub links nothing.
    """
    out = _root(tmp_path)
    hub = HubPage("Course hub: Demo", "Every page.", (), (), ())
    page = write(hub, out / "09_report" / f"{TODAY}_report_demo_hub_v1.html")
    (link,) = re.findall(r'class="glossary-link" href="([^"]+)"', unpacked(page.read_text()))
    assert link == f"{TODAY}_report_glossary_v1.html" and (page.parent / link).is_file()


def test_the_glossary_page(tmp_path):
    """Every entry anchored by id, by topic, with its fields; search and theme; offline.

    Mutations: the typical block removed from the template; the confirm pill removed; the
    see_also links pointed at the term, not the id; data-search built with the page's labels.
    """
    out = _root(tmp_path)
    write(_module(), out / "01_gather" / f"{TODAY}_gather_demo_overview_v1.html")
    target = out / "09_report" / f"{TODAY}_report_glossary_v1.html"
    html = unpacked(target.read_text())
    assert target.stat().st_size < MAX_BYTES
    ids = re.findall(r'<article class="entry" id="([^"]+)"', html)
    assert ids == [e.id for _, entries in glossary.grouped() for e in entries]
    assert len(ids) == 189
    assert len(re.findall(r'<section class="glossary-group"', html)) == 14
    gc = re.search(r'<article class="entry" id="gc_content" data-search="([^"]*)">(.*?)</article>',
                   html, re.S)  # fmt: skip
    for label in ("Also written:", "RNA-seq and DNA-seq:", "Typical values:", "Source:"):
        assert label in gc[2] and label not in gc[1], label
    assert "Lander et al. 2001" in gc[1] and "GC distribution" in gc[1]
    planted = glossary.Entry("made_up", "made-up term", "data", "Short.", "Long.", confirm=True)
    unverified = render(GlossaryPage("Glossary", (("data", (planted,)),)))
    assert '<span class="pill" data-kind="unverified">Not yet verified</span>' in unverified
    see_also = set(re.findall(r'<a href="#([^"]+)">', html)) - {
        f"topic-{c}" for c in glossary.CATEGORIES
    }
    assert see_also <= set(ids)
    assert 'id="ctx-glossary-search"' in html and 'id="ctx-theme-toggle"' in html
    assert not re.search(r'\b(src|href)="(https?:)?//', html)


MARKED = "Reads per [[FASTQ]] file; no test."


def test_captions_print_their_display_text_in_pdf_and_svg(tmp_path):
    """A caption's mark is a term on the page and its display text in the PDF and the SVG;
    place_caption wraps the display text, and the plotly figure holds no mark.

    Mutation: plain() removed from place_caption (the drawn caption keeps its brackets, and
    write_figure refuses it as a mark outside the caption).
    """
    from matplotlib.figure import Figure  # noqa: PLC0415

    from ctx.viz.annotate import place_caption  # noqa: PLC0415

    drawn = Figure(figsize=(3.35, 2))
    place_caption(drawn, "Reads per [[phred|Q score]] bin; no test.", 3.35, 2)
    assert drawn.texts[-1].get_text() == "Reads per Q score bin; no test."
    out = _root(tmp_path) / "02_qc"
    data = dataclasses.replace(demo_bar(counts()), caption=MARKED)
    figure_stem = stem("qc", "marked_demo", 1, day=date.today())
    files = write_figure(data, out, figure_stem, title="Demo", data_releases={})
    for path in files:
        if path.suffix == ".svg":
            svg = path.read_text()
            assert "[[" not in svg and "Reads per FASTQ file" in svg
        if path.suffix == ".pdf":
            printed = " ".join((p.extract_text() or "") for p in PdfReader(path).pages)
            assert "[[" not in printed and "FASTQ" in printed
    page = unpacked((out / f"{figure_stem}.html").read_text())
    assert re.search(r'<figcaption>Reads per <span class="term">FASTQ<button', page)
    plot = re.search(r'<div class="fig-plot">(.*?)<figcaption>', page, re.S)[1]
    assert "FASTQ" not in plot and "Samples" in plot
    assert (out.parent / "09_report" / f"{TODAY}_report_glossary_v1.html").is_file()


class _StaticMark:
    """demo_bar with one marked text in its PDF and SVG only, none in its plotly figure."""

    def __init__(self, data):
        self.data = data

    def __getattr__(self, name):
        return getattr(self.data, name)

    def static(self, theme):
        fig = self.data.static(theme)
        fig.text(0.5, 0.5, "[[FASTQ]]")
        return fig


@pytest.mark.parametrize(
    ("change", "title", "error"),
    [
        ({"value_label": "[[FASTQ]] files"}, "Demo", "a glossary mark outside the caption"),
        ({}, "Demo [[FASTQ]]", "a glossary mark outside the caption"),
        ("static", "Demo", "a glossary mark outside the caption"),
        ({"caption": "Reads per [[no such term]] file."}, "Demo", "not in the glossary"),
    ],
)
def test_a_figure_refuses_a_mark_it_cannot_draw_before_any_file(tmp_path, change, title, error):
    """A mark in an axis label, the page title or a drawn text prints raw in the plotly figure
    or the page, and a caption mark must name an entry: each fails before any file is written.

    Mutations: the page-title check removed (the title case writes its page); the check of the
    matplotlib texts removed (the static-only text passes); the caption's marks resolved only
    when the page is written (the PDF and SVG are on disk first). The axis-label case is also
    caught by the matplotlib check; test_figure_block_refuses_a_mark_in_the_plotly_figure
    covers the plotly check alone.
    """
    base = demo_bar(counts())
    data = _StaticMark(base) if change == "static" else dataclasses.replace(base, **change)
    out = tmp_path / "work"
    with pytest.raises(ValueError, match=error):
        write_figure(data, out, stem("qc", "refused", 1, day=date.today()), title=title,
                     data_releases={})  # fmt: skip
    assert not out.exists()


def test_figure_block_refuses_a_mark_in_the_plotly_figure():
    """figure_block, the door module pages use too, reads the plotly figure's text.

    Mutation: the plotly check removed from figure_block: the axis label passes.
    """
    from ctx.viz.save import figure_block  # noqa: PLC0415

    data = dataclasses.replace(demo_bar(counts()), value_label="[[FASTQ]] files")
    with pytest.raises(ValueError, match="a glossary mark outside the caption: .*FASTQ"):
        figure_block(data, "demo")


def test_a_figure_outside_a_ctx_root_is_written_whole(tmp_path):
    """``ctx figure`` will write into the reader's own folder (Phase 2D): a figure with a
    marked caption is written there whole, its tooltips with no glossary link.

    Mutation: page_path raising outside a module folder, as it did: the PDF and SVG were
    written, then the page refused.
    """
    data = dataclasses.replace(demo_bar(counts()), caption=MARKED)
    out = tmp_path / "reader_work"
    files = write_figure(data, out, stem("qc", "marked_demo", 1, day=date.today()),
                         title="Demo", data_releases={})  # fmt: skip
    assert len(files) == 7 and all(f.is_file() for f in files)
    page = unpacked(next(f for f in files if f.suffix == ".html").read_text())
    (payload,) = re.findall(r'id="ctx-glossary">(.*?)</script>', page, re.S)
    assert json.loads(payload)["href"] == "" and 'class="glossary-link"' not in page


def test_help_prints_the_display_text():
    """--help prints a card's meaning with its marks as display text.

    Mutation: HELP built from the raw meaning: '[[Phred]]' appears in --help.
    """
    from ctx.cli import app  # noqa: PLC0415

    done = CliRunner().invoke(app, ["qc", "--help"], terminal_width=200)
    assert done.exit_code == 0
    assert "[[" not in done.output and "Phred quality" in done.output


def test_the_kit_ships_the_glossary_and_start_links_it(tmp_path):
    """A kit copies the glossary its pages link, and START.html links the newest one.

    Mutation: course.page without the glossary: START.html has no glossary link.
    """
    out = _root(tmp_path)
    write(_module(), out / "01_gather" / f"{TODAY}_gather_demo_overview_v1.html")
    older = out / "09_report" / "2026-01-01_report_glossary_v1.html"
    older.write_text("<html><head><title>Glossary</title></head><body></body></html>")
    kit = tmp_path / "kit"
    done = subprocess.run(
        [sys.executable, str(BUILD_KIT), "--out", str(kit), f"--real=step={out.parent}"],
        capture_output=True, text=True, check=False,
    )  # fmt: skip
    assert done.returncode == 0, done.stderr
    shipped = f"step/outputs/09_report/{TODAY}_report_glossary_v1.html"
    assert (kit / shipped).is_file()
    home = (kit / course.INDEX).read_text()
    assert f'href="{shipped}"' in home


def test_a_mark_in_code_is_code():
    """A mark inside backticks is literal code: no icon, no tooltip, no lookup, and plain()
    leaves it as written.

    Mutation: marks() reading the whole text, code included: phred's tooltip is embedded, and
    the unknown name in code raises.
    """
    html = render(_module(purpose="Type `[[phred]]` to mark a term.", input="Plain.",
                          method="Type `[[your term]]` too."))  # fmt: skip
    assert "<code>[[phred]]</code>" in html and "data-term=" not in html
    assert 'id="ctx-glossary"' not in html
    assert marks.plain("Type `[[phred]]` or [[phred|Q]].") == "Type `[[phred]]` or Q."


def test_render_refuses_an_unknown_mark():
    """A page built in Python (M5, M6, the hub, the preview) meets no lint: rendering it is the
    guard, and a mark that names no entry fails it.

    Mutation: Terms.html and used() both skipping an unknown mark: the page renders it as text.
    """
    with pytest.raises(glossary.GlossaryError, match="'no such term' is not in the glossary"):
        render(_module(purpose="Check [[no such term]]."))


def test_blank_render_has_no_link_but_every_tooltip_names_its_term():
    """Rendered without write(), the tooltip data holds each term and short, no href.

    Mutation: tooltips() carrying every entry, not only the marked ones.
    """
    block = FigureBlock(go.Figure(), "From [[ENA]].", [], "d.csv", text=FigureText(title="T"))
    from ctx.report.builder import FigurePage  # noqa: PLC0415

    html = render(FigurePage("Alone", block))
    payload = json.loads(re.findall(r'id="ctx-glossary">(.*?)</script>', html, re.S)[0])
    assert payload["href"] == "" and list(payload["entries"]) == ["ena"]
    assert payload["entries"]["ena"]["term"] == "ENA"
