"""ctx.report.genes and its page twin ctxSelect (templates/_select.js.j2), on the same cases.

A gene selection is one definition in two languages (D-049): the Python reads list files and
builds pages, the page script reads links and pastes. Each case runs through both, and each
test names the mutation it catches.
"""

import json
import random
import re
import shutil
import subprocess
import time
import tomllib
from dataclasses import asdict
from pathlib import Path
from typing import Any

import pytest

from ctx.report import genes, text

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src" / "ctx" / "report" / "templates" / "_select.js.j2"
UI = ROOT / "src" / "ctx" / "report" / "content" / "ui.toml"
SRC = "tier1:ndn_gcsf"
# Real rows of Step-1's ndn_gcsf results (2026-09-30): the registered box genes, three
# symbols that name several tested genes, a gene named by its own ID, and C1orf122.
ROWS = (
    ("ENSG00000093100.13", "ENSG00000093100", 0.10242479854602089),
    ("ENSG00000111335.15", "OAS2", -3.1362422796506406),
    ("ENSG00000115415.23", "STAT1", -1.7785845403824907),
    ("ENSG00000120217.15", "CD274", 1.5774714245251056),
    ("ENSG00000126709.17", "IFI6", -1.9681595226651722),
    ("ENSG00000154451.17", "GBP5", -2.7642207381294788),
    ("ENSG00000163568.19", "AIM2", 1.164900969975241),
    ("ENSG00000187608.11", "ISG15", -2.317508935611896),
    ("ENSG00000197982.15", "C1orf122", 0.5785698818402946),
    ("ENSG00000202019.1", "Y_RNA", -2.8682360421378954),
    ("ENSG00000206341.7", "HLA-H", -1.6772959009915387),
    ("ENSG00000206990.1", "Y_RNA", -0.09045895463957457),
    ("ENSG00000274008.1", "Metazoa_SRP", -1.3847018250355694),
    ("ENSG00000275070.1", "Metazoa_SRP", -0.8729625197531826),
    ("ENSG00000278700.1", "Metazoa_SRP", -0.8667508277735625),
    ("ENSG00000310469.1", "HLA-H", -1.397357919365731),
)
TABLE = genes.table(*(list(column) for column in zip(*ROWS, strict=True)))
TABLE_IDS = [row[0] for row in ROWS]
UNSIGNED = genes.table(TABLE_IDS, list(TABLE.symbols))

HARNESS = r"""
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
const t = input.table;
const index = ctxSelect.index(t.numbers, t.symbols, t.lfc);
const own = (u) => ctxSelect.index(u.numbers, u.symbols, u.lfc);
const run = {
  select: (src, ids) => ctxSelect.select(src, ids, index),
  selectIn: (src, ids, u) => ctxSelect.select(src, ids, own(u)),
  unsigned: (src, ids) => ctxSelect.select(src, ids, ctxSelect.index(t.numbers, t.symbols)),
  encode: (sel) => ctxSelect.encode(sel),
  decode: (fragment) => ctxSelect.decode(fragment),
  carries: (fragment) => ctxSelect.carries(fragment),
  parse: (pasted) => ctxSelect.parse(pasted, index),
  parseIn: (pasted, u) => ctxSelect.parse(pasted, own(u)),
  listText: (sel, time) => ctxSelect.listText(sel, index, time),
};
console.log(JSON.stringify(input.calls.map(([name, ...args]) => {
  try { return { ok: run[name](...args) }; } catch (e) {
    if (!e.key) return { error: true };
    return { refused: { key: e.key, values: e.values, message: ctxSelect.message(e, input.ui) } };
  }
})));
"""


def _plain(value: Any) -> Any:
    """A Python result as JSON gives the page's back: dataclasses as dicts, tuples as lists."""
    plain = asdict(value) if hasattr(value, "__dataclass_fields__") else value
    return json.loads(json.dumps(plain))


def _own(table: dict[str, Any]) -> genes.GeneTable:
    """A table given as the harness's columns, for a case the Step-1 rows cannot make."""
    ids = [genes.ensembl(n) for n in table["numbers"]]
    return genes.table(ids, table["symbols"], table.get("lfc"))


def _python(name: str, *args: Any) -> dict[str, Any]:
    """One call of the Python definition, in the harness's shape."""
    run = {
        "select": lambda src, ids: genes.select(src, ids, TABLE),
        "selectIn": lambda src, ids, own: genes.select(src, ids, _own(own)),
        "unsigned": lambda src, ids: genes.select(src, ids, UNSIGNED),
        "encode": lambda sel: genes.encode(genes.Selection(sel["src"], tuple(sel["genes"]))),
        "decode": genes.decode,
        "carries": genes.carries,
        "parse": lambda pasted: genes.parse(pasted, TABLE),
        "parseIn": lambda pasted, own: genes.parse(pasted, _own(own)),
    }
    try:
        return {"ok": _plain(run[name](*args))}
    except genes.SelectionError as e:
        return {"refused": {"key": e.key, "values": e.values, "message": str(e)}}
    except ValueError:
        return {"error": True}


def _node(calls: list[list[Any]], version: int = genes.VERSION) -> list[dict[str, Any]]:
    """The same calls through ctxSelect under node, one process for all of them; ``version``
    replaces the script's VERSION, as a format change would."""
    node = shutil.which("node")
    assert node, "nodejs is pinned in the ctx-dev environment; run tests with pixi"
    script = SCRIPT.read_text().replace("{% raw %}", "").replace("{% endraw %}", "")
    script = re.sub(r"VERSION: [0-9]+,", f"VERSION: {version},", script, count=1)
    payload = {
        "calls": calls,
        "table": {"numbers": list(TABLE.numbers), "symbols": list(TABLE.symbols),
                  "lfc": list(TABLE.lfc or ())},
        "ui": {"select": text.raw("ui", "select")},
    }  # fmt: skip
    done = subprocess.run([node, "-e", script + HARNESS], input=json.dumps(payload),
                          capture_output=True, text=True, check=True)  # fmt: skip
    return json.loads(done.stdout)


def _both(calls: list[list[Any]]) -> list[dict[str, Any]]:
    """Each call's result, asserted equal in both languages."""
    page = _node(calls)
    for call, js in zip(calls, page, strict=True):
        assert js == _python(*call), call[0]
    return page


def test_select_strips_the_version_then_signs_then_drops_repeats_then_sorts():
    """The order test. Y_RNA's ENSG00000202019 enters versioned, unversioned and as its number;
    ENSG00000115415 in lower case; CD274, up, as -120217, a signed number as decode gives one.
    Catches: the version not stripped or repeats dropped before it is (two Y_RNA entries), the
    sign dropped, a number's own sign looked up or kept (CD274 "not in the gene table", or
    down), or genes sorted by their signed value (-202019 first)."""
    ids = ["ENSG00000202019.1", "ENSG00000202019", "ENSG00000163568.19", -120217,
           "ensg00000115415.23", 202019]  # fmt: skip
    (picked, fragment) = _both([["select", SRC, ids],
                                ["encode", _plain(genes.select(SRC, ids, TABLE))]])  # fmt: skip
    assert picked["ok"] == {"src": SRC, "genes": [-115415, 120217, 163568, -202019]}
    assert fragment["ok"] == "#v=1&src=tier1:ndn_gcsf&n=4&g=-115415,120217,163568,-202019&end"


def test_select_sorts_by_number_past_2_to_the_32_and_signs_zero_up():
    """JS lists an object's integer keys in ascending order only below 2**32 - 1, so the order
    test alone cannot show that the page sorts: numbers above it, entered unsorted, do. None of
    the 139,599 rows of Step-1's 10 results tables has a log2 FC of 0 (measured for GS.1), so 0
    and -0 come from a made-up table. Catches: the sort removed from either language, the
    input sorted as text in place of the output, and a log2 FC of 0 or -0 signed down (<= 0
    for < 0) in either."""
    own = {
        "numbers": [99_999_999_999, 5_000_000_000, 10_000_000_000, 115415, 1, 2],
        "symbols": ["A", "B", "C", "D", "E", "F"],
        "lfc": [1.0, -1.0, 1.0, 1.0, 0.0, -0.0],
    }
    (picked,) = _both([["selectIn", "a:b", [*own["numbers"][:4], 2, 1], own]])  # fmt: skip
    assert picked["ok"]["genes"] == [1, 2, 115415, -5_000_000_000, 10**10, 99_999_999_999]


def test_no_selection_is_made_that_no_fragment_carries():
    """Selection refuses at construction what decode would refuse, in both languages, so the
    page cannot write a link it then refuses. Catches: the page's select or encode taking any
    source (Tier1:NDN-gcsf, empty), gene 0 or 10**11, or numbers out of order."""
    stat1 = ["ENSG00000115415"]
    calls = [["select", "Tier1:NDN-gcsf", stat1], ["select", "", stat1],
             *(["encode", {"src": src, "genes": g}] for src, g in (
                 ("Tier1:ndn-gcsf", [-115415]), (SRC, [0]), (SRC, [10**11]),
                 (SRC, [120217, -115415]), (SRC, [115415, -115415])))]  # fmt: skip
    assert all(result == {"error": True} for result in _both(calls))


def test_select_refuses_a_gene_it_cannot_sign():
    """Catches: a table without log2 FC signing every gene up, a gene outside the table
    entering with no sign, a _PAR_Y copy taken as its X copy, and Python alone folding a long
    s (U+017F) into the S of ENSG."""
    results = _both([["unsigned", SRC, ["ENSG00000115415"]], ["select", SRC, ["ENSG00000999999"]],
                     ["select", SRC, ["ENSG00000115415_PAR_Y"]],
                     ["select", SRC, ["EN" + chr(0x17F) + "G00000115415"]]])  # fmt: skip
    assert all(r == {"error": True} for r in results)


def _selections() -> list[genes.Selection]:
    """The empty list, one gene each way, the extreme numbers, and 17,830 genes: as many as
    Step-1 tested over its 10 contrasts, drawn from 1 to 400,000, where its 17,830 numbers run
    up to 310,796 (measured for GS.1)."""
    rng = random.Random(7)
    many = sorted(rng.sample(range(1, 400_001), 17_830))
    return [
        genes.Selection(SRC, ()),
        genes.Selection(SRC, (115415,)),
        genes.Selection("a:b", (-115415,)),
        genes.Selection(SRC, (-1, 99_999_999_999)),
        genes.Selection(SRC, tuple(n if rng.random() < 0.5 else -n for n in many)),
    ]


def test_encoding_and_decoding_are_exact_inverses_in_both_languages():
    """Catches: a decoder that drops the sign or reorders; an encoder whose n is not the gene
    count or that leaves out the end mark; and either language drifting from the other on any
    of these."""
    selections = _selections()
    fragments = [genes.encode(s) for s in selections]
    calls = [["encode", _plain(s)] for s in selections] + [["decode", f] for f in fragments]
    results = _both(calls)
    assert [r["ok"] for r in results[: len(fragments)]] == fragments
    for selection, fragment, decoded in zip(selections, fragments, results[len(fragments):],
                                            strict=True):  # fmt: skip
        assert genes.decode(fragment) == selection
        assert genes.encode(genes.decode(fragment)) == fragment
        assert decoded["ok"] == _plain(selection)
    assert fragments[0] == "#v=1&src=tier1:ndn_gcsf&n=0&g=&end"
    assert fragments[3] == "#v=1&src=tier1:ndn_gcsf&n=2&g=-1,99999999999&end"


def test_the_version_is_defined_once_in_each_language(monkeypatch):
    """A format change bumps VERSION alone: each language's encoder then writes what its own
    decoder reads, and a link of the old version is refused as such. Catches: a fragment
    pattern holding the version as a literal (the page's '#v=1' did, so version 2 links were
    refused as select.form) or built from VERSION when the module loads."""
    selection, old = genes.Selection(SRC, (-115415,)), "#v=1&src=tier1:ndn_gcsf&n=1&g=-115415&end"
    monkeypatch.setattr(genes, "VERSION", 2)
    new = genes.encode(selection)
    assert (new, genes.decode(new)) == (old.replace("v=1", "v=2"), selection)
    with pytest.raises(genes.SelectionError, match="version 1 gene list"):
        genes.decode(old)
    page = _node([["encode", _plain(selection)], ["decode", new], ["decode", old]], version=2)
    assert [page[0]["ok"], page[1]["ok"]] == [new, _plain(selection)]
    assert page[2]["refused"]["key"] == "select.version"


GOOD = "#v=1&src=tier1:ndn_gcsf&n="
REFUSALS = (
    ("", "select.form"),
    ("#results", "select.form"),
    ("#v=", "select.form"),
    ("v=1&src=tier1:ndn_gcsf&n=0&g=&end", "select.form"),
    ("#v=2&src=tier1:ndn_gcsf&n=0&g=&end", "select.version"),
    ("#v=2", "select.version"),
    ("#v=01&src=tier1:ndn_gcsf&n=0&g=&end", "select.version"),
    # A lone surrogate, which Python's argv holds for a byte that is not UTF-8.
    ("#v=1" + chr(0xDCFF) + "&src=tier1:ndn_gcsf&n=0&g=&end", "select.version"),
    # Cut inside the source, before the g field, at a comma, and just after one.
    ("#v=1&src=tier1:nd", "select.cut"),
    ("#v=1&src=tier1:ndn_gcsf&n=0", "select.cut"),
    (GOOD + "3&g=-115415,120217", "select.cut"),
    (GOOD + "1&g=", "select.cut"),
    (GOOD + "3&g=-115415,120217,", "select.cut"),
    # Cut inside the last number of ...,99999999999&end: the count holds and every number is
    # well formed and in order, so before the end mark this read as gene 999999.
    (GOOD + "2&g=-115415,999999", "select.cut"),
    # Cut inside the end mark, and text after it.
    (GOOD + "1&g=115415&en", "select.cut"),
    (GOOD + "0&g=&end&t=black", "select.cut"),
    ("#v=1&n=0&src=tier1:ndn_gcsf&g=&end", "select.form"),
    ("#v=1&src=tier1:ndn_gcsf&n=0&end", "select.form"),
    ("#v=1&src=tier1:ndn_gcsf&n=0&g=&t=black&end", "select.form"),
    # The field pattern is anchored at the end: a second end mark is no gene list's end.
    (GOOD + "1&g=115415&end&t=black&end", "select.form"),
    ("#v=1&src=Tier1:ndn_gcsf&n=0&g=&end", "select.source"),
    ("#v=1&src=tier1&n=0&g=&end", "select.source"),
    ("#v=1&src=tier1:&n=0&g=&end", "select.source"),
    ("#v=1&src=tier1:ndn%20gcsf&n=0&g=&end", "select.source"),
    ("#v=1&src=&n=0&g=&end", "select.source_empty"),
    (GOOD + "03&g=-115415,120217,163568&end", "select.n"),
    (GOOD + "&g=&end", "select.n_empty"),
    # The end mark is there, so a count that differs is an edit: genes removed, or added.
    (GOOD + "3&g=-115415,120217&end", "select.count"),
    (GOOD + "1&g=&end", "select.count"),
    (GOOD + "1&g=115415,120217&end", "select.count"),
    # The count before the genes: this list is both short and malformed.
    (GOOD + "3&g=115415,0&end", "select.count"),
    (GOOD + "0&g=" + chr(0xD800) + "&end", "select.count"),
    # The count is of the items that hold something, so an empty one is named as such.
    (GOOD + "1&g=115415,&end", "select.number_empty"),
    (GOOD + "2&g=115415,,120217&end", "select.number_empty"),
    (GOOD + "1&g=0&end", "select.number"),
    (GOOD + "1&g=0115415&end", "select.number"),
    (GOOD + "1&g=+115415&end", "select.number"),
    (GOOD + "1&g=1e5&end", "select.number"),
    (GOOD + "1&g=-0&end", "select.number"),
    (GOOD + "1&g=123456789012&end", "select.number"),
    # 1, then 15415 in Arabic-Indic digits: Python's \d takes them, the page's does not.
    (GOOD + "1&g=1" + "".join(chr(0x660 + int(d)) for d in "15415") + "&end", "select.number"),
    (GOOD + "2&g=-115415,115415&end", "select.duplicate"),
    (GOOD + "3&g=115415,120217,115415&end", "select.duplicate"),
    (GOOD + "2&g=120217,-115415&end", "select.order"),
    (GOOD + "1&g=" + "1" * 499_980 + "&end", "select.long"),
    # 250,000 emoji: 500,001 UTF-16 units, as location.hash counts, but 250,001 code points.
    ("#" + "\U0001f600" * 250_000, "select.long"),
)


def test_decode_refuses_with_the_reason_and_the_same_words():
    """Each refusal gives the same key, values and ui.toml sentence in both languages, and
    every [select] sentence is some refusal's. Catches: any one check removed (its case then
    decodes or fails later under another key), the end mark not checked, checked after the
    fields (a cut is then select.form) or anywhere in place of at the end, the n check
    skipped so an edited list is used, a count refused only when short so extra genes pass,
    the genes checked before the count, an empty item counted so it is called a count
    mismatch, an empty field quoted as a gap ("its source ;"), \\d in place of [0-9], a bare
    #v=2 not read as a version, the length counted in code points or by a strict UTF-16
    encoding (Python raised UnicodeEncodeError on a lone surrogate), and a sentence the page
    fills differently."""
    results = _both([["decode", fragment] for fragment, _ in REFUSALS])
    keys = [r["refused"]["key"] if "refused" in r else r for r in results]
    assert keys == [key for _, key in REFUSALS]
    said = {f: r["refused"]["message"] for (f, _), r in zip(REFUSALS, results, strict=True)}
    assert said[GOOD + "1&g="] == (
        "The link does not end with its end mark: it was cut short or has text added after it. "
        "Copy the whole link again."
    )
    assert said[GOOD + "1&g=&end"] == (
        "The link gives a gene count of 1, but its list holds 0: it was edited."
    )
    assert said["#v=1&src=&n=0&g=&end"] == (
        "The link names no source; a source is a dataset and a contrast, as dataset:contrast."
    )
    assert results[-1]["refused"]["values"] == {"length": "500,001", "limit": "500,000"}
    words = tomllib.loads(UI.read_text())["select"]
    assert {f"select.{k}" for k in words} == {key for _, key in REFUSALS}


def test_a_link_cut_anywhere_after_its_version_is_refused_as_cut():
    """Every cut of an encoded link, from just after '#v=1' to just before its last character,
    is select.cut in both languages. Catches: the end mark not written or not checked, the
    check after the fields (each cut is then select.form), and a cut inside the last number
    read as another gene."""
    fragment = genes.encode(genes.Selection(SRC, (-115415, 120217, 99_999_999_999)))
    cuts = [fragment[:end] for end in range(len("#v=1"), len(fragment))]
    results = _both([["decode", cut] for cut in cuts])
    assert [r["refused"]["key"] if "refused" in r else r for r in results] == (
        ["select.cut"] * len(cuts)
    )


def test_a_count_longer_than_any_link_holds_is_refused_unformatted():
    """A link of 500,000 characters holds fewer than 250,000 numbers, so a count of 7 digits or
    more is refused before the thousands grouping, a regex quadratic in the digits: before
    the bound, a 499,971-digit count took 94.1 s in node 26.10.0 (measured for GS.1).
    Catches: the bound on the count's digits removed from either language (the short case is
    then select.count, which fails before the long case runs)."""
    (short,) = _both([["decode", GOOD + "1" + "0" * 6 + "&g=&end"]])
    assert short["refused"]["values"] == {"found": "1000000", "digits": "6"}
    longest = GOOD + "1" * (genes.LIMIT - len(GOOD) - len("&g=&end")) + "&g=&end"
    start = time.perf_counter()
    (refused,) = _both([["decode", longest]])
    assert (len(longest), refused["refused"]["key"]) == (genes.LIMIT, "select.n")
    assert time.perf_counter() - start < 1


def _fragment(length: int) -> tuple[str, list[int]]:
    """A fragment of exactly ``length`` characters and its numbers: 11-digit genes, some
    10-digit. Built by hand, not by ``encode``, which refuses one over the limit.

    Each 11-digit gene takes 12 characters with its comma, so the shortest count that reaches
    ``length`` overshoots by at most 11, and that many genes are given 10 digits.
    """
    count = 1
    while len(f"{GOOD}{count}&g=&end") + 12 * count - 1 < length:
        count += 1
    over = len(f"{GOOD}{count}&g=&end") + 12 * count - 1 - length
    numbers = [*range(10**10 - over, 10**10), *range(10**10, 10**10 + count - over)]
    fragment = f"{GOOD}{count}&g=" + ",".join(map(str, numbers)) + "&end"
    assert len(fragment) == length
    return fragment, numbers


def test_the_limit_is_500000_characters_inclusive_both_ways():
    """A fragment of 500,000 characters decodes and encodes back to itself; one character more
    is refused by decode and, so a page never writes a link it refuses, by encode. Catches:
    > turned into >=, the page holding another limit, either encoder writing what decode
    refuses (Step-1's 78,470 gene IDs, all signed down, take 625,384 and the end mark's 4),
    and a section anchor such as #results taken for a selection."""
    (at, kept), (over, too_many) = _fragment(genes.LIMIT), _fragment(genes.LIMIT + 1)
    results = _both([["decode", at], ["encode", {"src": SRC, "genes": kept}],
                     ["decode", over], ["encode", {"src": SRC, "genes": too_many}],
                     ["carries", at], ["carries", "#results"]])  # fmt: skip
    assert (len(results[0]["ok"]["genes"]), results[1]["ok"]) == (len(kept), at)
    assert [r["refused"]["key"] for r in results[2:4]] == ["select.long"] * 2
    assert results[3]["refused"]["values"] == {"length": "500,001", "limit": "500,000"}
    assert [results[4]["ok"], results[5]["ok"]] == [True, False]


PASTED = """STAT1, isg15;Ifi6\t"OAS2"
ENSG00000120217.15 | ensg00000163568
y_rna HLA-H
NOTAGENE stat1 ENSG00000999999.1 notagene
ENSG00000202019.1\tY_RNA
ENSG00000274008 Metazoa_SRP  # a saved list's line: the symbol labels its ID
# source: tier1:ndn_gcsf
C1orf122 ENSG00000093100
(GBP5.) ENSG00000206990.1/ENSG00000275070 -ENSG00000278700.1.
"""


def test_paste_reads_names_and_ids_and_refuses_ambiguous_symbols():
    """Catches: an ambiguous symbol resolved to its first gene (prototype C's error in the GS
    plan, 4 wrong picks in 3,703); a case-sensitive match (isg15, y_rna); a version not
    stripped; a comment read as names; a saved line's symbol refused as ambiguous (no label
    rule); a symbol of another gene on an ID's line taken as a label (C1orf122); and a full
    stop, dash or slash beside a name read as part of it (GBP5., IDs joined by a slash)."""
    (result,) = _both([["parse", PASTED]])
    assert result["ok"] == {
        "found": [115415, 187608, 126709, 111335, 120217, 163568, 202019, 274008, 197982, 93100,
                  154451, 206990, 275070, 278700],
        "ambiguous": [["y_rna", ["ENSG00000202019", "ENSG00000206990"]],
                      ["HLA-H", ["ENSG00000206341", "ENSG00000310469"]]],
        "missing": ["NOTAGENE", "ENSG00000999999.1"],
    }  # fmt: skip


def test_a_slash_joins_names_unless_a_symbol_holds_it():
    """THRA1/BTR, ENSG00000235300.4, is the one symbol with a slash among GENCODE v36's 60,660
    (one GDC STAR counts file) and Step-1's 78,470 (none). Catches: every slash split, so
    THRA1/BTR is two unknown names; and no slash split, so THRA/STAT1 is one."""
    own = {"numbers": [235300, 126351], "symbols": ["THRA1/BTR", "THRA"]}
    (result,) = _both([["parseIn", "thra1/btr THRA/STAT1.", own]])
    assert result["ok"] == {"found": [235300, 126351], "ambiguous": [], "missing": ["STAT1"]}


def test_a_saved_list_reads_back_as_the_selection(tmp_path):
    """The page writes the Save-list file and read_list reads it, ambiguous symbols included.
    Catches: the label rule missing (Y_RNA, HLA-H and Metazoa_SRP refused), the count check
    skipped so a list cut short is used, a header line read as names, a '# count:' after a
    gene on its line read as a header, a byte-order mark, CRLF or CR alone breaking the read
    in either language (a CR-only file read as one line, its first '#' hiding the rest), and
    a file that is not UTF-8 raising UnicodeDecodeError in place of GeneListError."""
    selection = genes.select(SRC, TABLE_IDS, TABLE)
    (written,) = _node([["listText", _plain(selection), "2026-10-01T09:00:00.000Z"]])
    saved = written["ok"]
    assert saved.splitlines()[:4] == [
        "# source: tier1:ndn_gcsf",
        "# time: 2026-10-01T09:00:00.000Z",
        "# count: 16",
        "ENSG00000093100\tENSG00000093100",
    ]
    path = tmp_path / "list.txt"
    endings = (saved, chr(0xFEFF) + saved.replace("\n", "\r\n"), saved.replace("\n", "\r"))
    for content, pasted in zip(endings, _both([["parse", e] for e in endings]), strict=True):
        path.write_text(content, encoding="utf-8", newline="")
        read = genes.read_list(path, TABLE)
        assert (read.source, read.time) == (SRC, "2026-10-01T09:00:00.000Z")
        assert (read.genes.ambiguous, read.genes.missing) == ((), ())
        assert genes.select(read.source, read.genes.found, TABLE) == selection
        assert genes.select(SRC, pasted["ok"]["found"], TABLE) == selection
    path.write_text(saved.rsplit("ENSG", 1)[0])
    with pytest.raises(genes.GeneListError, match=r"'# count: 16', but 15 lines name genes"):
        genes.read_list(path, TABLE)
    for header, problem in (
        ("# source: Step-1\n", "not <dataset>:<contrast>"),
        ("# count: 1\n# count: 1\nSTAT1\n", "two '# count:' lines"),
    ):
        path.write_text(header)
        with pytest.raises(genes.GeneListError, match=problem):
            genes.read_list(path, TABLE)
    path.write_text("STAT1 # count: 5\nCD274\n")
    assert genes.read_list(path, TABLE).genes.found == (115415, 120217)
    # Excel's "Unicode Text" export is UTF-16; a Parquet file starts with PAR1.
    for raw in ("STAT1\nCD274\n".encode("utf-16"), b"PAR1\x15\x04\xec\x00"):
        path.write_bytes(raw)
        with pytest.raises(genes.GeneListError, match="not UTF-8 text"):
            genes.read_list(path, TABLE)


def test_a_hand_written_list_needs_no_header(tmp_path):
    """Catches: the header made compulsory, so a list typed by hand is refused."""
    path = tmp_path / "mine.txt"
    path.write_text("stat1\nCD274, AIM2\n")
    read = genes.read_list(path, TABLE)
    assert (read.source, read.time, read.genes.found) == (None, None, (115415, 120217, 163568))


def test_the_table_refuses_rows_with_no_identity_of_their_own():
    """GENCODE v36, the model of GDC's STAR counts, has 44 _PAR_Y rows in one breast tumour's
    counts file (measured for GS.1), CD99's ENSG00000002586.20_PAR_Y among them: its number is
    the X copy's. Catches: a _PAR_Y copy merged into its X copy, two versions of one gene kept
    as two rows, and a gene with no finite log2 FC signed as up."""
    with pytest.raises(ValueError, match="not an Ensembl gene ID"):
        genes.table(["ENSG00000002586.20_PAR_Y"], ["CD99"])
    with pytest.raises(ValueError, match=r"more than one row: 1, such as \['ENSG00000002586'\]"):
        genes.table(["ENSG00000002586.20", "ENSG00000002586.19"], ["CD99", "CD99"])
    for lfc in (float("nan"), None, float("inf")):
        with pytest.raises(ValueError, match="without a finite log2 FC"):
            genes.table(["ENSG00000002586.20"], ["CD99"], [lfc])
    with pytest.raises(ValueError, match="differ in length"):
        genes.table(["ENSG00000002586.20"], [])
    with pytest.raises(ValueError, match="ascend without repeats"):
        genes.Selection(SRC, (120217, -115415))
