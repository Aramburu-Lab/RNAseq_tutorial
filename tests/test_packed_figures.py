"""Packed figures (Z.1): a page's figure JSON and rows as raw DEFLATE in base64, drawn by the
same inline scripts, so the page draws what it drew unpacked. The inflate is checked against
zlib under node; the packed page against the unpacked one, in Python and in Chromium. Each test
names the mutation it catches."""

import base64
import hashlib
import json
import random
import re
import shutil
import subprocess
import zlib
from dataclasses import replace
from datetime import date
from pathlib import Path

import pytest
from synth import demo_strip
from test_export import _drawn
from test_export_browser import CHROMIUM, _browse

from ctx.io.naming import stem
from ctx.report.builder import (
    FigureText,
    Group,
    ModulePage,
    ParameterCard,
    Reproduction,
    SelfCheck,
    render,
)
from ctx.report.figures import pack_json, unpack_json
from ctx.viz.save import figure_block, write_figure

ROOT = Path(__file__).resolve().parents[1]
UNPACK = ROOT / "src" / "ctx" / "report" / "templates" / "_unpack.js.j2"
DAY = date(2026, 10, 21)
PACKED_PLOT = re.compile(r'Plotly\.newPlot\(\s*"([^"]+)",\s*\.\.\.ctxUnpack\("([^"]*)"\)')
ROWS = re.compile(r'<script type="([^"]+)" id="(ctx-data-\d+)">([^<]*)</script>')

NODE = """const fs = require("fs"), crypto = require("crypto");
eval(fs.readFileSync(process.argv[1], "utf8") + ";globalThis.ctxUnpack = ctxUnpack;");
const cases = JSON.parse(fs.readFileSync(0, "utf8"));
console.log(JSON.stringify(cases.map(c => crypto.createHash("sha256")
  .update(ctxUnpack.inflate(new Uint8Array(Buffer.from(c, "base64")))).digest("hex"))));
"""


def _inputs() -> list[bytes]:
    """Text and bytes that take every block type zlib writes: stored, fixed and dynamic codes,
    long and short matches, and the empty input."""
    rng = random.Random(42)
    figure = json.dumps(demo_strip().interactive().to_plotly_json(), default=str).encode()
    two_letters = bytes(rng.choice(b"ab") for _ in range(100_000))
    return [b"", b"a", b"[1,2,3]", rng.randbytes(70_000), two_letters, b"x" * 300_000, figure * 20]


@pytest.mark.skipif(shutil.which("node") is None, reason="nodejs is pinned in the dev environment")
def test_unpack_inflates_what_zlib_deflates():
    """Every level and strategy zlib has, through the page's inflate under node, gives the
    input back byte for byte. Mutations caught: a length or distance table off by one; the
    stored block's partial byte kept; a repeat code's count base wrong."""
    cases, wanted = [], []
    for data in _inputs():
        for level in (0, 1, 6, 9):
            for strategy in (zlib.Z_DEFAULT_STRATEGY, zlib.Z_FIXED, zlib.Z_HUFFMAN_ONLY,
                             zlib.Z_RLE):  # fmt: skip
                deflate = zlib.compressobj(level, zlib.DEFLATED, -15, 9, strategy)
                cases.append(base64.b64encode(deflate.compress(data) + deflate.flush()).decode())
                wanted.append(hashlib.sha256(data).hexdigest())
    done = subprocess.run([shutil.which("node"), "-e", NODE, str(UNPACK)], input=json.dumps(cases),
                          capture_output=True, text=True, check=True)  # fmt: skip
    assert json.loads(done.stdout) == wanted


def test_pack_json_round_trips():
    """Mutation caught: a zlib header written where the page reads raw DEFLATE."""
    text = json.dumps({"x": [1.5, None, "</script>"]})
    packed = pack_json(text)
    assert unpack_json(packed) == json.loads(text) and "<" not in packed
    assert zlib.decompress(base64.b64decode(packed), -15) == text.encode()


def _page(out: Path, pack: bool) -> str:
    """A module page of two groups: one strip, then two strips as tabs."""
    blocks = []
    for i in range(3):
        data = replace(demo_strip(), caption=f"Strip {i}; no test.")
        figure_stem = stem("clinical", f"demo_{i}", 1, day=DAY)
        text = FigureText(f"Strip {i}", reading=("Each point is one made-up sample.",))
        write_figure(data, out, figure_stem, title=f"Strip {i}", data_releases={}, text=text)
        blocks.append(figure_block(data, figure_stem, text))
    groups = (Group("One", (blocks[0],)), Group("Two", tuple(blocks[1:]), tabs=True))
    built = ModulePage(
        title="M10 Clinical: Demo", purpose="Purpose.", input="Input.", method="Method.",
        parameters=(ParameterCard("--x", "A card.", "1", "Moves it."),), results=groups,
        interpretation="Interpretation.", questions=(SelfCheck("A question.", "An answer."),),
        pitfalls=("A pitfall.",), pack_figures=pack,
        reproduce=Reproduction("ctx demo", 1.0, 1, 1, "laptop", "none", {"python": "3.12"}),
    )  # fmt: skip
    return render(built)


def unpacked_figures(html: str) -> tuple[list, dict]:
    """(div id, data, layout, config) of every packed figure, and every packed rows table."""
    plots = [(m[1], *unpack_json(m[2])) for m in PACKED_PLOT.finditer(html)]
    rows = {m[2]: unpack_json(m[3]) for m in ROWS.finditer(html)
            if m[1] == "application/octet-stream"}  # fmt: skip
    return plots, rows


@pytest.fixture(scope="module")
def pages(tmp_path_factory):
    out = tmp_path_factory.mktemp("packed_figures") / "outputs" / "12_clinical"
    out.mkdir(parents=True)
    return out, _page(out, pack=False), _page(out, pack=True)


def test_a_packed_page_carries_the_same_figures_and_rows(pages):
    """Unpacked, every figure's data, layout and config and every rows table equal the
    unpacked page's, and the packed page carries no raw figure JSON. Mutations caught: the
    config argument dropped from the packed arguments; rows packed from another figure."""
    _, plain, packed = pages
    plots, rows = unpacked_figures(packed)
    assert [p[:3] for p in plots] == _drawn(plain)
    assert len(plots) == 3  # noqa: PLR2004 - three strips
    assert not re.search(r'Plotly\.newPlot\(\s*"[^"]+",\s*\[', packed)
    configs = re.findall(r"Plotly\.newPlot\((?:.|\n)*?(\{\"displaylogo\"[^}]*\})", plain)
    assert [json.dumps(p[3], sort_keys=True) for p in plots] == [
        json.dumps(json.loads(c), sort_keys=True) for c in configs
    ]
    plain_rows = {m[2]: json.loads(m[3]) for m in ROWS.finditer(plain)}
    assert rows == plain_rows and len(rows) == 3  # noqa: PLR2004
    assert len(packed.encode()) < len(plain.encode())
    assert "var ctxUnpack" not in plain and packed.count("var ctxUnpack") == 1


DUMP = """(async () => {
  await new Promise(r => setTimeout(r, 1500));
  const plots = Array.from(document.querySelectorAll('.js-plotly-plot')).map(d => ({
    id: d.id, data: JSON.stringify(d.data), layout: JSON.stringify(d.layout),
    full: JSON.stringify(d._fullData.map(t => [t.type, t.x && t.x.length, t.name]))}));
  const rows = Array.from(document.querySelectorAll('script[id^="ctx-data-"]')).map(s =>
    JSON.stringify(s.type === 'application/octet-stream' ? ctxUnpack(s.textContent)
                                                          : JSON.parse(s.textContent)));
  return {plots, rows};
})()"""


@pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
def test_chromium_draws_the_packed_page_as_the_unpacked_one(pages, tmp_path):
    """Each plot's data and layout as plotly holds them after drawing, and every rows table,
    are equal on the two pages. Mutation caught: the inflate's output cut short."""
    out, plain, packed = pages
    paths = {}
    for name, html in (("plain", plain), ("packed", packed)):
        paths[name] = out / f"{name}.html"
        paths[name].write_text(html)
    tasks = [{"name": n, "url": p.as_uri(), "width": 1280, "height": 900, "act": DUMP}
             for n, p in paths.items()]  # fmt: skip
    seen = _browse(tasks, tmp_path)
    assert seen["plain"]["errors"] == [] and seen["packed"]["errors"] == []
    assert len(seen["plain"]["acted"]["plots"]) == 3  # noqa: PLR2004
    assert seen["packed"]["acted"] == seen["plain"]["acted"]
