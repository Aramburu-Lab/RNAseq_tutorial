"""M9.8 and C3.1: gene cards packed to the page cap, the rest named with their command.

The 50-card default list made a 29,472,416-byte TCGA breast page that the cap refused (C3).
Each test names the mutation it catches; every one was planted and seen to fail.
"""

import re

import polars as pl
import pytest
from genecard_root import SLUG, Stub, make, responses
from test_genecard import _cards, _register, _settings

from ctx.genecard import page as card_page
from ctx.genecard import run as module
from ctx.report.builder import PACKED, packed_bytes, unpacked

NPTX2, ACTB = "ENSG00000106236", "ENSG00000075624"


@pytest.fixture
def root(tmp_path):
    return tmp_path / "root", make(tmp_path / "root")


def test_fit_keeps_the_cards_whose_measured_page_fits():
    """The page over k cards is measured, never estimated: a cap that holds the first three
    keeps three, and the names of the cards left out count too.

    Mutations caught: the last card left in when one more byte is over; the excess covered by
    the first packed card's bytes instead of the dropped one's; the first card dropped."""
    costs = [500, 100, 200, 300, 400]  # the first is drawn whole; the rest packed

    def pages(k: int) -> tuple[str, list[int]]:
        named = 7 * (len(costs) - k)  # the left-out line grows as cards leave
        return "x" * (1000 + sum(costs[:k]) + named), costs[1:k]

    assert card_page.fit(pages, 5, cap=2500) == 5  # noqa: PLR2004
    assert card_page.fit(pages, 5, cap=2499) == 4  # noqa: PLR2004
    assert card_page.fit(pages, 5, cap=1814) == 3  # 1800 + 14 left-out bytes  # noqa: PLR2004
    assert card_page.fit(pages, 5, cap=1813) == 2  # noqa: PLR2004
    assert card_page.fit(pages, 5, cap=10) == 1


def _run(root, monkeypatch, cap=None):
    folder, spec = root
    _register(folder, spec)
    if cap is not None:
        monkeypatch.setattr(card_page, "MAX_BYTES", cap)
    options = module.Options(spec, genes=("NPTX2", "ACTB"))
    return folder, module.run(_settings(folder), _cards(), "ctx genecard", SLUG, options,
                              Stub(responses()))  # fmt: skip


def test_cards_after_the_first_travel_packed_and_read_back_whole(root, monkeypatch):
    """Two cards: NPTX2 drawn whole, ACTB packed, and the packed HTML reads back to ACTB's card;
    two builds pack to the same bytes.

    Mutations caught: every card packed or none; a gzip timestamp, which makes two builds
    differ."""
    _, written = _run(root, monkeypatch)
    html = written[0].read_text()
    assert len(packed_bytes(html)) == 1
    (blob,) = (m[0] for m in PACKED.finditer(html))
    assert "Exon use of ACTB" not in html.replace(blob, "")
    assert "Exon use of ACTB" in unpacked(html) and "Exon use of NPTX2" in html
    _, again = _run(root, monkeypatch)
    assert [m[2] for m in PACKED.finditer(again[0].read_text())] == [
        m[2] for m in PACKED.finditer(html)]  # fmt: skip


def test_a_card_past_the_cap_is_named_with_its_command_and_writes_no_files(root, monkeypatch):
    """With the cap just under the two-card page, ACTB leaves: no figure file, no table row,
    named in its own group with the command that builds its card by Ensembl ID.

    Mutations caught: figure files written for a card the page left out; the tables keeping
    its rows; the left-out group or its command missing."""
    folder, full = _run(root, monkeypatch)
    two = len(full[0].read_bytes())
    for path in (folder / "outputs" / "11_genecard").iterdir():
        path.unlink()
    _, written = _run(root, monkeypatch, cap=two - 1)
    html = unpacked(written[0].read_text())
    assert packed_bytes(written[0].read_text()) == []
    assert "1 of the 2 listed genes have no card here" in html
    assert "Left out, in list order: ACTB." in html
    command = re.search(r'<code id="ctx-code-1">([^<]*)</code>', html)
    assert command and f"--genes {ACTB}" in command[1] and "--offline" not in command[1]
    assert not [p for p in written if "_actb_" in p.name]
    evidence = pl.read_parquet(next(p for p in written if p.name.endswith("_evidence_v1.parquet")))
    assert set(evidence["symbol"]) == {"NPTX2"}
