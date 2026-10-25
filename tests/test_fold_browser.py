"""A folded caption and reading guide open on a click, and a term inside them still works
(D-119), in Chromium. Skipped where no Chromium is found; set ``CTX_CHROMIUM``.
"""

import dataclasses
import tempfile
from datetime import date
from pathlib import Path

import pytest
from synth import counts, demo_bar
from test_export_browser import CHROMIUM, _browse

from ctx.io.naming import stem
from ctx.report.builder import FigureText
from ctx.viz.save import write_figure

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
CAPTION = "Reads per [[FASTQ]] file. Each [[phred|Q score]] bin; no test."
# A real press and release at the centre of what the expression selects, scrolled into view.
PRESS = """(() => {{
  const el = document.querySelector({sel});
  el.scrollIntoView({{block: "center", behavior: "instant"}});
  const b = el.getBoundingClientRect();
  const x = b.left + b.width / 2, y = b.top + b.height / 2;
  return [['mousePressed', x, y], ['mouseReleased', x, y]];
}})()"""
STATE = r"""(() => ({
  caption: document.querySelector("figcaption details").open,
  reading: document.querySelector(".reading details").open,
  tip: !!document.querySelector("#ctx-tip:not([hidden])"),
  shown: document.querySelector("figcaption").innerText.replace(/\s+/g, " "),
}))()"""


@pytest.fixture(scope="module")
def clicked(tmp_path_factory) -> dict:
    """The page loaded, then a click on the term in the caption's summary, on the summary's
    text, and on the reading guide's title, each step's state read after it."""
    out = tmp_path_factory.mktemp("fold") / "root" / "outputs" / "02_qc"
    out.mkdir(parents=True)
    data = dataclasses.replace(demo_bar(counts()), caption=CAPTION)
    name = stem("qc", "fold_demo", 1, day=date.today())
    text = FigureText("Demo", "A lead.", ("Look at the bars.", "A short bar is a bad run."))
    write_figure(data, out, name, title="Demo", data_releases={}, text=text)
    steps = [{"at": PRESS.format(sel=f'"{sel}"'), "then": STATE}
             for sel in ("figcaption summary .term button", "figcaption summary",
                         ".reading summary")]  # fmt: skip
    task = {"name": "fold", "url": (out / f"{name}.html").as_uri(), "width": 1280,
            "height": 1000, "act": STATE, "mouse": steps}  # fmt: skip
    with tempfile.TemporaryDirectory() as folder:
        return _browse([task], Path(folder))["fold"]


def test_both_load_folded_and_open_in_place(clicked):
    """Loaded, only the first sentence shows; a click on it opens the rest, one on the guide's
    title opens the guide.

    Mutations: the caption's details open at load; the guide's details open at load."""
    assert clicked["errors"] == []
    loaded = clicked["acted"]
    assert not loaded["caption"] and not loaded["reading"], loaded
    assert "no test" not in loaded["shown"] and "Reads per FASTQ file." in loaded["shown"]
    _, summary, guide = clicked["moused"]
    assert summary["caption"] and "no test" in summary["shown"], summary
    assert guide["reading"], guide


def test_a_term_in_the_first_sentence_opens_its_tooltip_not_the_caption(clicked):
    """The term's button takes its own click: the tooltip opens, the caption stays folded."""
    term = clicked["moused"][0]
    assert term["tip"] and not term["caption"], term
