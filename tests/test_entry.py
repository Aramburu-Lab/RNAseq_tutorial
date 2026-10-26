"""GS.2: each gene-list-reading module writes ``outputs/NN_<module>/index.html``; the kit walks
it and refuses a consuming module's folder without one; a packed part's links are walked too.

Each test names the mutation it catches; every one was planted and seen to fail.
"""

import base64
import gzip
import re
import subprocess
import sys

import pytest
from test_kit import CHECK, NEW, SCRIPT, a_root, link_hubs, load, page

from ctx.io.naming import ENTRY, lint
from ctx.report import builder
from ctx.report.builder import unpacked
from ctx.report.entry import CONSUMERS, write_entry


@pytest.fixture(scope="module")
def checker():
    return load(CHECK, "check_kit")


def _overview(folder, slug, name, day=NEW):
    page(folder / f"{day}_genecard_{slug}_overview_v1.html", f"M9 Gene cards: {name}")


def test_one_dataset_forwards_and_two_are_listed(tmp_path):
    """The entry page links each dataset's newest overview; one link forwards, carrying the
    fragment; nothing is fetched and no plotly is carried.

    Mutations caught: an older overview linked; the fragment left off the forwarded link; a
    module that reads no gene list given an entry page."""
    folder = tmp_path / "outputs" / "11_genecard"
    _overview(folder, "tcga_brca", "TCGA breast", day="2026-09-01")
    _overview(folder, "tcga_brca", "TCGA breast")
    html = unpacked(write_entry(folder, "genecard").read_text())
    assert re.findall(r'<a href="([^"]+)" data-entry data-dataset', html) == [
        f"{NEW}_genecard_tcga_brca_overview_v1.html"]  # fmt: skip
    assert 'location.replace(links[0].getAttribute("href"))' in html
    assert 'a.getAttribute("href") + location.hash' in html
    assert "plotly.js v" not in html and len(html.encode()) < 500_000  # noqa: PLR2004
    _overview(folder, "tcga_kirc", "TCGA kidney")
    html = unpacked(write_entry(folder, "genecard").read_text())
    # Two datasets: each link names its own, so a list's src= picks the page (Rytis
    # 2026-10-26: the reader was asked which dataset they had selected in).
    assert re.findall(r'data-entry data-dataset="([^"]+)"', html) == ["tcga_brca", "tcga_kirc"]
    assert 'getAttribute("data-dataset") === from' in html
    with pytest.raises(ValueError, match="reads no gene list"):
        write_entry(folder, "de")


def test_the_naming_lint_accepts_index_html_in_a_module_folder_only(tmp_path):
    """Mutation: index.html accepted anywhere, so a stray one deeper down passes."""
    folder = tmp_path / "11_genecard"
    _overview(folder, "s", "S")
    (folder / ENTRY).write_text("x")
    (folder / "deeper").mkdir()
    (folder / "deeper" / ENTRY).write_text("x")
    assert lint(folder)[1] == ["deeper/index.html"]


def test_check_kit_refuses_a_consuming_folder_without_its_entry_page(tmp_path, checker):
    """Mutations caught: the rule dropped; it applied to a module that reads no list; the
    twin's list of consuming modules drifting from ctx.report.entry's."""
    assert checker.CONSUMERS == CONSUMERS
    assert checker.ENTRY == ENTRY
    kit = tmp_path / "kit"
    _overview(kit / "step" / "outputs" / "11_genecard", "s", "S")
    page(kit / "step" / "outputs" / "05_de" / f"{NEW}_de_s_overview_v1.html", "M5: S")
    assert checker.entry_problems(kit) == [
        "step/outputs/11_genecard: no index.html, which a link to this module names"]  # fmt: skip
    (kit / "step" / "outputs" / "11_genecard" / ENTRY).write_text("x")
    assert checker.entry_problems(kit) == []


def test_a_link_inside_a_packed_part_is_walked(tmp_path, checker):
    """A card drawn on demand links its figure files from inside its gzipped HTML.

    Mutations caught: the link walk reading the page's text only, so a missing file inside a
    packed part passes and a kit drops the files it links; the twin unpacking differently."""
    inner = '<a href="missing_v1.pdf">PDF</a><a href="there.csv">CSV</a>'
    blob = base64.b64encode(gzip.compress(inner.encode(), mtime=0)).decode()
    html = f'<p>x</p><script type="application/octet-stream" id="ctx-packed-2">{blob}</script>'
    assert checker.unpacked(html) == builder.unpacked(html) == f"<p>x</p>{inner}"
    root = tmp_path / "kit"
    root.mkdir()
    (root / "there.csv").write_text("a\n")
    (root / "p.html").write_text(html)
    assert checker.link_problems(root, [root / "p.html"]) == [
        "p.html: missing_v1.pdf does not resolve"]  # fmt: skip


def _kit(tmp, root):
    out = tmp / "kit"
    done = subprocess.run([sys.executable, str(SCRIPT), "--out", str(out), f"--real=step={root}"],
                          capture_output=True, text=True, check=False)  # fmt: skip
    return out, done


def test_the_kit_ships_entry_pages_and_refuses_a_missing_one(tmp_path):
    """Mutations caught: the undated entry page left out of the kit; a consuming module's
    folder without one shipped."""
    root = a_root(tmp_path / "real", "Study", "study")
    folder = root / "outputs" / "12_clinical"
    page(folder / f"{NEW}_clinical_study_overview_v1.html", "M10 Clinical: Study")
    link_hubs(root)
    out, done = _kit(tmp_path / "a", root)
    assert done.returncode != 0 and "12_clinical: no index.html" in done.stderr
    assert not out.exists()
    write_entry(folder, "clinical")
    out, done = _kit(tmp_path / "b", root)
    assert done.returncode == 0, done.stderr
    assert (out / "step" / "outputs" / "12_clinical" / ENTRY).is_file()
