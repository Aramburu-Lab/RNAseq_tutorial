"""F.1: M9's figures, each redrawn from its spec's inputs alone in a clean folder.

The made-up root draws every kind M9 has: the evidence view of two contrasts, a box, F1, the
tissues, exons, junctions and immune panels. Every one has a spec and writes M9's data table,
PDF and SVG byte for byte. Each test names its mutation.
"""

import html
import re
from pathlib import Path

import polars as pl
import pytest
from genecard_root import SLUG, Stub, make, own_files, responses
from test_figure_panel import panels
from test_figures import copied, drawn_alike, redrawn, specs
from test_genecard import _cards, _register, _settings

from ctx.dataset import load
from ctx.de.run import newest
from ctx.figures import kinds, kinds_m9
from ctx.figures.panel import TABS
from ctx.figures.params import ParamsError, load_params
from ctx.figures.spec import parse
from ctx.genecard import boxes
from ctx.genecard import run as module
from ctx.report import cap
from ctx.report.builder import unpacked

FOLDER = ("11_genecard",)
GENES = ("UPA", "NPTX2", "ACTB")
EVERY = {"evidence_view", "box", "contrasts", "tissues", "exons", "junctions", "immune"}


def _built(root: Path, genes: tuple[str, ...] = GENES) -> Path:
    """ctx genecard on the made-up root with two contrasts, own samples and NPTX2's responses."""
    spec = make(root)
    _register(root, spec)
    own_files(root, ("a_case", "b_case"))
    results = newest(root / "outputs" / "05_de", SLUG, "a_case_results")
    for suffix in (".parquet", ".csv"):
        source = results.with_suffix(suffix)
        source.with_name(source.name.replace("a_case", "b_case")).write_bytes(source.read_bytes())
    spec.write_text(spec.read_text() + '\n[[contrasts]]\nname = "b_case"\nvariable = "condition"'
                    '\nreference = "control"\ntested = "case"\nn = [4, 4]\npaired = true\n'
                    'label = "Case against control, paired"\n')  # fmt: skip
    module.run(_settings(root), _cards(), "ctx genecard", SLUG,
               module.Options(spec, genes=genes), Stub(responses()))  # fmt: skip
    return spec


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """Every M9 figure kind drawn on the made-up root: (root, dataset file)."""
    root = tmp_path_factory.mktemp("m9") / "root"
    return root, _built(root)


def _drawn(root: Path) -> set[tuple[str, str]]:
    """(kind, gene) of every figure M9 wrote, read from its file names."""
    found = set()
    for pdf in (root / "outputs" / FOLDER[0]).glob("*_v1_white.pdf"):
        what = pdf.name[11:].removeprefix(f"genecard_{SLUG}_").removesuffix("_v1_white.pdf")
        kind = "evidence_view" if what == "evidence_view" else what.split("_", 1)[0]
        found.add((kind, what.removeprefix(f"{kind}_") if kind != what else ""))
    return found


def _named(found) -> set[tuple[str, str]]:
    """(kind, gene) of each spec."""
    return {(parse(s.figure_id)[2], getattr(s.params, "gene", "")) for s in found}


def _redraws_alike(root: Path, found, tmp_path: Path) -> None:
    """Each spec redrawn in its own clean folder equals M9's files."""
    for spec in found:
        module_csv = next(root.glob(f"outputs/11_genecard/{spec.stem}.csv"))
        csv = redrawn(copied(spec, root, tmp_path / spec.stem), spec)
        assert drawn_alike(csv, module_csv) == [], spec.stem


def test_every_figure_has_a_spec_and_redraws_byte_for_byte(built, tmp_path):
    """Every figure M9 drew has a spec, each redrawn from M9's tables alone equals its files.
    Catches a kind whose stem no longer parses, a card rebuilt from another gene's rows, the
    MANE exons left out of the junction panel, the dropped samples read as none, the F1 values
    read under the wrong contrast, and the tissue order taken from the table's rows."""
    root, spec = built
    found = specs(root, FOLDER, spec)
    drawn = _drawn(root)
    assert {k for k, _ in drawn} == EVERY, drawn
    assert _named(found) == drawn
    assert _panels(root) == sorted(f"ctx figure {s.figure_id} --params" for s in found)
    _redraws_alike(root, found, tmp_path)


def _panels(root: Path) -> list[str]:
    """The ``ctx figure`` line of each Reproduce panel on M9's page, packed cards opened.

    Each panel prints the line once, in the steps its tabs share (Y.1): the tabs must list the
    same figures, and one tab's list is the page's.
    """
    (path,) = (root / "outputs" / FOLDER[0]).glob("*_overview_v1.html")
    shown = unpacked(path.read_text())
    each = {tab: sorted(re.search(r"ctx figure \S+ --params", p).group(0)
                        for p, _ in panels(shown, tab)) for tab in TABS}  # fmt: skip
    assert len({tuple(lines) for lines in each.values()}) == 1, each
    assert len(re.findall(r"ctx figure \S+ --params", html.unescape(shown))) == len(each[TABS[0]])
    return each[TABS[0]]


def test_a_page_short_of_cards_redraws_the_whole_list(tmp_path, monkeypatch):
    """With the page cap holding one card, the evidence view still covers all three genes and
    the kept card's box its p over both tested genes, each redrawn byte for byte. Catches the
    evidence view redrawn from the page's evidence table, which holds the kept card alone, and
    the box family written over the kept cards."""
    # Every page shares cap.MAX_BYTES, so the fit alone is held to one card here.
    monkeypatch.setattr(cap, "kept", lambda pages, parts, least=0: 1)
    root = tmp_path / "root"
    spec = _built(root)
    found = [s for s in specs(root, FOLDER, spec)
             if parse(s.figure_id)[2] in {"evidence_view", "box"}]  # fmt: skip
    assert _named(found) == {("evidence_view", ""), ("box", "upa_900011")}
    shown = pl.read_parquet(newest(root / "outputs" / FOLDER[0], SLUG, "evidence"))
    listed = pl.read_parquet(newest(root / "outputs" / FOLDER[0], SLUG, "evidence_listed"))
    assert (shown["symbol"].n_unique(), listed["symbol"].n_unique()) == (1, 3)
    family = pl.read_parquet(newest(root / "outputs" / FOLDER[0], SLUG, "own_family"))
    assert family["symbol"].to_list() == ["UPA", "ACTB"]
    _redraws_alike(root, found, tmp_path / "redraw")


def test_the_immune_panel_names_the_release_the_run_named(built, tmp_path):
    """The immune panel's release comes from M9's options table: changed there, the redrawn
    caption names the new one. Catches the release taken from hpa.RELEASE, which
    --hpa-release overrides."""
    root, spec = built
    one = next(s for s in specs(root, FOLDER, spec) if parse(s.figure_id)[2] == "immune")
    folder = copied(one, root, tmp_path / "release")
    options = folder / next(i.file for i in one.inputs if i.name == "genecard_options")
    assert pl.read_parquet(options)["hpa_release"].to_list() == ["v25"]
    pl.DataFrame({"hpa_release": ["v99"]}).write_parquet(options)
    raw = (folder / "params.toml").read_text().splitlines()
    (folder / "params.toml").write_text(
        "\n".join(x for x in raw if not x.startswith("genecard_options = ") or ".parquet" in x)
        + "\n"
    )
    html = redrawn(folder, one).with_suffix(".html").read_text()
    assert "HPA v99" in html and "HPA v25" not in html


def test_a_card_is_named_as_m9_names_it(built, tmp_path):
    """A description that is no card's name matches no kind, and a params file naming one is
    refused with the rule. Catches the card pattern left out of the check, which would draw
    from whatever number the name ends in, or stop inside the drawer."""
    root, spec_path = built
    spec = load(spec_path)
    assert kinds.matched("genecard", "tissues_nptx2_106236", spec) is not None
    for bad in ("tissues_nptx2", "tissues_NPTX2_106236", "tissue_nptx2_106236"):
        assert kinds.matched("genecard", bad, spec) is None, bad
    one = next(s for s in specs(root, FOLDER, spec_path) if parse(s.figure_id)[2] == "tissues")
    folder = copied(one, root, tmp_path / "bad")
    raw = (folder / "params.toml").read_text()
    assert raw.count(f'gene = "{one.params.gene}"') == 1
    (folder / "params.toml").write_text(raw.replace(f'gene = "{one.params.gene}"',
                                                    'gene = "nptx2"'))  # fmt: skip
    with pytest.raises(ParamsError, match="gene: 'nptx2' is not a card's name"):
        load_params(folder / "params.toml", one.figure_id)


def test_contrasts_sharing_their_words_leave_the_own_panels_without_a_spec(built, tmp_path):
    """Two contrasts given the same words: M9's values table names both alike, so the box and
    F1 have no spec, while the tissues keep theirs. Catches the words check left out, which
    would redraw one contrast's samples under the other's name."""
    _, spec_path = built
    words = boxes.names(load(spec_path))["b_case"]
    raw = spec_path.read_text()
    assert raw.count("n = [3, 3]\n") == 1
    alike = tmp_path / spec_path.name
    # Labelled alike: the dataset file's own labels, which M5 never tells apart.
    alike.write_text(raw.replace("n = [3, 3]\n", f'n = [3, 3]\nlabel = "{words}"\n'))
    spec = load(alike)
    assert set(boxes.names(spec).values()) == {words}
    assert not kinds_m9.words_unique(spec)
    for kind in kinds_m9.OWN:
        assert kinds.matched("genecard", f"{kind}_upa_900011", spec) is None, kind
    assert kinds.matched("genecard", "tissues_upa_900011", spec) is not None
