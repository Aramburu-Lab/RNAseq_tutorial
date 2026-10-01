"""The glossary tooltips and search rule, run under node (D-048), in the style of
test_page_script.py: templates/_glossary.js.j2 is the file the page embeds.

The tooltip wiring runs on a stand-in document: elements with attributes, listeners, focus and
bubbling, enough for the events a keyboard, a mouse and a touch screen send. Each test names
the mutation it catches.
"""

import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from ctx.report.builder import HubPage, render

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src" / "ctx" / "report" / "templates" / "_glossary.js.j2"

DOCUMENT = r"""
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
let doc;
function dispatch(target, type, props, bubbles) {
  const event = { type, target, key: props.key, pointerType: props.pointerType,
                  relatedTarget: props.relatedTarget || null, preventDefault() {} };
  for (let n = target; n; n = bubbles ? n.parentNode : null) {
    (n.on[type] || []).forEach((f) => f(event));
  }
}
class El {
  constructor(tag) {
    Object.assign(this, { tagName: tag.toUpperCase(), children: [], parentNode: null,
      attrs: {}, on: {}, style: {}, hidden: false, textContent: "", className: "", id: "",
      rect: { left: 0, top: 0, right: 0, bottom: 0 } });
  }
  get firstChild() { return this.children[0] || null; }
  get offsetWidth() { return this.hidden ? 0 : input.panel.width; }
  get offsetHeight() { return this.hidden ? 0 : input.panel.height; }
  setAttribute(k, v) { this.attrs[k] = String(v); }
  getAttribute(k) { return k in this.attrs ? this.attrs[k] : null; }
  removeAttribute(k) { delete this.attrs[k]; }
  appendChild(c) {
    if (c.parentNode) c.parentNode.removeChild(c);
    c.parentNode = this; this.children.push(c); return c;
  }
  removeChild(c) {
    this.children = this.children.filter((x) => x !== c); c.parentNode = null; return c;
  }
  contains(o) { for (; o; o = o.parentNode) if (o === this) return true; return false; }
  addEventListener(t, f) { (this.on[t] = this.on[t] || []).push(f); }
  getBoundingClientRect() { return this.rect; }
  focus() { doc.focusOn(this); }
  all() { return [this, ...this.children.flatMap((c) => c.all())]; }
}
doc = new El("#document");
doc.activeElement = null;
// The page's own box. innerWidth and innerHeight read larger while a panel is open, as on a
// phone whose layout grew with it, and as a classic scrollbar makes innerWidth read wider.
doc.documentElement = { clientWidth: input.view.width, clientHeight: input.view.height };
let panel = null;
const grown = () => (panel && !panel.hidden ? 200 : 0);
doc.createElement = (t) => new El(t);
doc.getElementById = (id) => doc.all().find((e) => e.id === id) || null;
doc.querySelector = () => null;
doc.querySelectorAll = () => doc.all().filter((e) =>
  e.tagName === "BUTTON" && e.className === "info" && e.getAttribute("data-term"));
doc.focusOn = (el) => {
  const prev = doc.activeElement;
  if (prev === el) return;
  if (prev) dispatch(prev, "focusout", { relatedTarget: el }, true);
  doc.activeElement = el;
  if (el) dispatch(el, "focus", {}, false);
};
const timers = [];
const win = { get innerWidth() { return input.view.width + grown(); },
  get innerHeight() { return input.view.height + grown(); }, on: {},
  setTimeout(f) { timers.push(f); return timers.length; },
  clearTimeout(id) { if (id) timers[id - 1] = null; },
  addEventListener(t, f) { (this.on[t] = this.on[t] || []).push(f); } };
const body = doc.appendChild(new El("body"));
const data = body.appendChild(new El("script"));
data.id = "ctx-glossary";
data.textContent = JSON.stringify(input.glossary);
const els = {};
for (const [name, term, rect] of input.terms) {
  const span = body.appendChild(new El("p")).appendChild(new El("span"));
  span.className = "term"; span.rect = rect;
  const button = span.appendChild(new El("button"));
  button.className = "info"; button.setAttribute("data-term", term);
  button.setAttribute("aria-expanded", "false");
  els["t" + name] = span; els["b" + name] = button;
}
els.out = body.appendChild(new El("button"));
panel = ctxGlossary.tips(doc, win);
const mouse = { pointerType: "mouse" }, touch = { pointerType: "touch" };
function enter(term, how) {
  dispatch(els["t" + term], "pointerenter", how, false);
  dispatch(els["b" + term], "pointerenter", how, false);
}
const target = (t) => (t === "link" ? panel.children[2] : els[t === "out" ? "out" : "b" + t]);
const act = {
  // Tab: the key goes to the focused element, then focus moves.
  tab: (t) => { dispatch(doc.activeElement || doc, "keydown", { key: "Tab" }, true);
                doc.focusOn(target(t)); },
  // Focus with no key: a screen reader or a script.
  focus: (t) => doc.focusOn(target(t)),
  // A press that ends in no click (dragged off, a right click): Chrome focuses the button.
  press: (t) => { dispatch(els["b" + t], "pointerdown", mouse, true); doc.focusOn(els["b" + t]); },
  // A press Safari gives no focus, or a touch that turned into a scroll.
  nudge: (t) => dispatch(els["b" + t], "pointerdown", touch, true),
  // A click on the open panel's text: focus goes to it if it is focusable, else nowhere.
  select: () => { dispatch(panel.children[1], "pointerdown", mouse, true);
                  doc.focusOn(panel.getAttribute("tabindex") === null ? null : panel);
                  dispatch(panel.children[1], "click", {}, true); },
  // The page scrolls by n px: every term moves up, and the document hears it.
  scroll: (n) => { for (const [name] of input.terms) {
                     const r = els["t" + name].rect;
                     els["t" + name].rect = { ...r, top: r.top - n, bottom: r.bottom - n }; }
                   dispatch(doc, "scroll", {}, false); },
  key: (k) => dispatch(doc.activeElement || doc, "keydown", { key: k }, true),
  enter: (t) => dispatch(els["b" + t], "click", {}, true),
  hover: (t) => enter(t, mouse),
  leave: (t) => { dispatch(els["b" + t], "pointerleave", mouse, false);
                  dispatch(els["t" + t], "pointerleave", mouse, false); },
  back: (t) => dispatch(els["t" + t], "pointerenter", mouse, false),
  tick: () => timers.splice(0).forEach((f) => f && f()),
  click: (t) => { enter(t, mouse); dispatch(els["b" + t], "pointerdown", mouse, true);
                  doc.focusOn(els["b" + t]); dispatch(els["b" + t], "click", {}, true); },
  tap: (t) => { enter(t, touch); dispatch(els["b" + t], "pointerdown", touch, true);
                doc.focusOn(els["b" + t]); dispatch(els["b" + t], "click", {}, true); },
  away: () => { doc.focusOn(null); dispatch(body, "click", {}, true); },
  clickout: () => dispatch(body, "click", {}, true),
  resize: () => (win.on.resize || []).forEach((f) => f()),
};
const out = input.steps.map(([verb, arg]) => {
  act[verb](verb === "scroll" ? Number(arg) : arg);
  const names = input.terms.map(([n]) => n);
  const holder = names.find((n) => els["t" + n].contains(panel));
  return {
    step: verb + " " + (arg || ""),
    open: names.filter((n) => els["b" + n].getAttribute("aria-expanded") === "true"),
    described: names.filter((n) => els["b" + n].getAttribute("aria-describedby") === panel.id),
    shown: !panel.hidden && holder !== undefined ? holder : null,
    role: panel.getAttribute("role"),
    words: panel.children.map((c) => c.textContent),
    tags: panel.children.map((c) => c.tagName),
    href: panel.children[2] ? panel.children[2].getAttribute("href") : null,
    side: panel.getAttribute("data-side"),
    at: [panel.style.left, panel.style.top],
    focus: doc.activeElement === els.out ? "out" : doc.activeElement === panel.children[2]
      ? "link" : doc.activeElement === panel ? "panel"
      : names.find((n) => doc.activeElement === els["b" + n]) || null,
  };
});
console.log(JSON.stringify(out));
"""

GLOSSARY = {
    "href": "../09_report/2026-09-30_report_glossary_v1.html",
    "link": "Glossary",
    "entries": {
        "fastq": {"term": "FASTQ file", "short": "Text file of reads."},
        "md5": {"term": "md5 checksum", "short": "A file's fingerprint."},
    },
}
# One term high on the left, one low on the right of an 800 by 800 view; a 300 by 120 panel.
TERMS = [["1", "fastq", {"left": 100, "top": 100, "right": 180, "bottom": 124}],
         ["2", "md5", {"left": 700, "top": 760, "right": 790, "bottom": 784}]]  # fmt: skip


def _node(harness: str, payload: dict) -> object:
    node = shutil.which("node")
    assert node, "nodejs is pinned in the ctx-dev environment; run tests with pixi"
    code = SCRIPT.read_text().replace("{% raw %}", "").replace("{% endraw %}", "")
    done = subprocess.run([node, "-e", code + harness], input=json.dumps(payload),
                          capture_output=True, text=True, check=True)  # fmt: skip
    return json.loads(done.stdout)


def _run(*steps: tuple[str, str]) -> list[dict]:
    payload = {"glossary": GLOSSARY, "terms": TERMS, "view": {"width": 800, "height": 800},
               "panel": {"width": 300, "height": 120},
               "steps": [list(s) for s in steps]}  # fmt: skip
    result = _node(DOCUMENT, payload)
    assert isinstance(result, list)
    return result


def test_keyboard_focus_opens_escape_closes_and_returns_focus():
    """Tab to a term opens its tooltip; Tab into the link keeps it; Escape closes it and gives
    focus back to the button without opening it again; Enter toggles.

    Mutations: the focus listener removed (Tab opens nothing); the returning flag ignored
    (Escape's refocus opens it again); the focusout check without contains() (Tab to the link
    closes it); the heading drawn as a strong, which Pitfalls' rule for strong recolours in
    the warning colour (its tags change).
    """
    tab, link, esc, on, off = _run(("tab", "1"), ("tab", "link"), ("key", "Escape"),
                                   ("enter", "1"), ("enter", "1"))  # fmt: skip
    assert (tab["open"], tab["described"], tab["shown"], tab["role"]) == (["1"], ["1"], "1",
                                                                          "tooltip")  # fmt: skip
    assert tab["words"] == ["FASTQ file", "Text file of reads.", "Glossary"]
    assert tab["tags"] == ["SPAN", "SPAN", "A"]
    assert tab["href"] == GLOSSARY["href"] + "#fastq"
    assert (link["shown"], link["focus"]) == ("1", "link")
    assert (esc["open"], esc["described"], esc["shown"], esc["focus"]) == ([], [], None, "1")
    assert (on["shown"], off["shown"]) == ("1", None)


def test_blur_closes_and_one_is_open_at_a_time():
    """Focus leaving the term closes it, whichever term takes the focus.

    Mutation: the focusout listener removed (Tab to another term leaves the first open).
    """
    one, two, three, left = _run(("tab", "1"), ("tab", "2"), ("tab", "1"), ("tab", "out"))
    assert (one["open"], two["open"], three["open"], left["open"]) == (["1"], ["2"], ["1"], [])
    assert two["described"] == ["2"] and left["shown"] is None


def test_hover_opens_and_leaving_closes_after_a_pause():
    """A mouse over the icon opens it; leaving the term closes it after the pause, and coming
    back within the pause keeps it; a click pins it open past a leave; hovering another term
    opens that one and closes the pinned one.

    Keyboard focus on a hover-opened tooltip keeps it open when the mouse leaves, as long as
    the focus stays (WCAG 1.4.13, persistent).

    Mutations: the pointerleave timer removed (it never closes); the pointerenter on the term
    not clearing the timer (coming back still closes); a click on a hover-opened tooltip
    closing it instead of pinning; hide() dropped from show() (two open at once); show() of
    an open tooltip upgrading only a click, not focus (the focused one closes on leave).
    """
    steps = _run(("hover", "1"), ("leave", "1"), ("tick", ""), ("hover", "1"), ("leave", "1"),
                 ("back", "1"), ("tick", ""), ("click", "1"), ("leave", "1"), ("tick", ""),
                 ("hover", "2"), ("tab", "2"), ("leave", "2"), ("tick", ""))  # fmt: skip
    shown = [s["shown"] for s in steps]
    assert shown == ["1", "1", None, "1", "1", "1", "1", "1", "1", "1", "2", "2", "2", "2"]
    assert steps[10]["open"] == ["2"] and steps[10]["described"] == ["2"]
    assert steps[-1]["focus"] == "2"


def test_tap_toggles_and_an_outside_tap_closes():
    """On a touch screen a tap opens and a second tap closes; hover is ignored for touch;
    a tap outside, or a resize, closes; so does a click outside a tooltip pinned with no
    focus, as Safari leaves a clicked button unfocused.

    Mutations: the pointing flag ignored (the tap's focus opens it, the click then closes);
    the outside-click listener removed (the unfocused pinned one stays); the resize listener
    removed.
    """
    steps = _run(("tap", "1"), ("tap", "1"), ("tap", "2"), ("away", ""), ("tap", "2"),
                 ("resize", ""), ("hover", "1"), ("enter", "1"), ("clickout", ""))  # fmt: skip
    shown = [s["shown"] for s in steps]
    assert shown == ["1", None, "2", None, "2", None, "1", "1", None]


def test_the_panel_flips_and_stays_in_the_view():
    """A term near the bottom opens its tooltip above; near the right edge it moves left.

    The panel is fixed, so it is placed in viewport coordinates, inside the page's own box
    (documentElement's client size) even while innerWidth and innerHeight read 200 px larger.

    Mutations: place() without the flip (the low term's panel runs below the view); view()
    from innerWidth and innerHeight (the low panel sits at left 692, below the term).
    """
    high, low = _run(("tab", "1"), ("tab", "2"))
    assert (high["side"], high["at"]) == ("below", ["100px", "130px"])
    # Above: 760 - 6 - 120 = 634; left 800 - 8 - 300 = 492.
    assert (low["side"], low["at"]) == ("above", ["492px", "634px"])


def test_the_panel_follows_its_term_when_the_page_scrolls():
    """A scroll anywhere (a card table's own included, heard in the capture phase) moves the
    fixed panel with its term.

    Mutation: the scroll listener removed: the panel stays at 130 while its term moves up 40.
    """
    opened, scrolled = _run(("tab", "1"), ("scroll", "40"))
    assert (opened["at"], scrolled["at"]) == (["100px", "130px"], ["100px", "90px"])
    assert scrolled["shown"] == "1"


def test_a_press_that_is_no_click_leaves_keyboard_focus_working():
    """A press that ends in no click leaves no trace: the next focus opens its tooltip.

    Chrome focuses the button on the press, which spends the pointing flag; Safari does not,
    and a touch can turn into a scroll, so the next key press clears it.

    Mutations: the flag not spent in the focus listener (after a drag off the icon, focus
    from a screen reader opens nothing); the keydown listener not clearing it (after a
    press Safari leaves unfocused, Tab to a term opens nothing).
    """
    steps = _run(("press", "1"), ("focus", "out"), ("focus", "2"), ("focus", "out"),
                 ("nudge", "1"), ("tab", "2"))  # fmt: skip
    assert [s["shown"] for s in steps] == [None, None, "2", None, None, "2"]


def test_a_click_on_the_open_tooltip_keeps_it():
    """Clicking the text of a pinned tooltip keeps it open, focus inside the term, so the
    reader can select its words; Escape still closes it and returns focus to the button.

    Mutation: the panel's tabindex removed (the press moves focus to the page body, and
    focusout closes the tooltip, as it did in Chrome and Firefox).
    """
    pinned, selected, esc = _run(("click", "1"), ("select", ""), ("key", "Escape"))
    assert (pinned["shown"], selected["shown"], selected["focus"]) == ("1", "1", "panel")
    assert (esc["shown"], esc["focus"]) == (None, "1")


PLACE = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
console.log(JSON.stringify({
  places: input.places.map(([a, s, v]) => ctxGlossary.place(a, s, v)),
  matches: input.matches.map(([t, q]) => ctxGlossary.matches(t, q)),
}));
"""


@pytest.mark.parametrize(
    ("anchor", "view", "want"),
    [
        ((100, 100, 124), (0, 800, 800), (130, 100, False)),  # room below
        ((100, 700, 724), (0, 800, 800), (574, 100, True)),  # no room below: above
        ((100, 700, 724), (0, 800, 780), (574, 100, True)),
        ((100, 150, 174), (60, 800, 200), (180, 100, False)),  # no room either way: below
        ((700, 100, 124), (0, 800, 800), (130, 492, False)),  # past the right edge: moved
        ((4, 100, 124), (0, 800, 800), (130, 8, False)),  # past the left margin: moved
    ],
)
def test_place(anchor, view, want):
    """The placement rule on numbers: below, above, clamped left and right, the top bar.

    Mutation: the view's top (the sticky bar) ignored: the fourth row opens above, under it.
    """
    left, top, bottom = anchor
    top_, right, bottom_ = view
    payload = {
        "places": [[{"left": left, "top": top, "bottom": bottom},
                    {"width": 300, "height": 120},
                    {"left": 0, "top": top_, "right": right, "bottom": bottom_}]],
        "matches": [],
    }  # fmt: skip
    out = _node(PLACE, payload)
    assert isinstance(out, dict)
    (spot,) = out["places"]
    assert (spot["top"], spot["left"], spot["above"]) == want


def test_search_matches_every_word_in_any_case():
    """The glossary page's search keeps an entry that holds every word of the query.

    Mutation: some() in place of every(): 'phred qqq' matches.
    """
    cases = [["Phred quality score. Q30 means", "phred q30"], ["Phred quality", "phred qqq"],
             ["Anything", ""], ["Anything", "  "]]  # fmt: skip
    out = _node(PLACE, {"places": [], "matches": cases})
    assert isinstance(out, dict)
    assert out["matches"] == [True, False, True, True]


SEARCH = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
const el = (props) => Object.assign({ hidden: false, attrs: {}, on: {},
  getAttribute(k) { return k in this.attrs ? this.attrs[k] : null; },
  addEventListener(t, f) { (this.on[t] = this.on[t] || []).push(f); } }, props);
const entries = input.entries.map(([text, words]) =>
  el({ textContent: text, attrs: { "data-search": words } }));
const group = el({ querySelector: () => entries.find((a) => !a.hidden) || null });
const box = el({ value: "" });
const count = el({ attrs: { "data-template": "{shown} of {total}" } });
ctxGlossary.search({
  getElementById: (id) => (id === "ctx-glossary-search" ? box : null),
  querySelectorAll: (s) => (s === "article.entry" ? entries : [group]),
  querySelector: () => count,
});
console.log(JSON.stringify(input.queries.map((q) => {
  box.value = q;
  box.on.input.forEach((f) => f());
  return [count.textContent, group.hidden];
})));
"""


def test_search_reads_the_entry_words_not_the_labels():
    """The glossary page's search reads each entry's data-search, its own words: the labels
    every entry prints ("Also written:", "Source:") match nothing.

    Mutation: search() reading the article's textContent: "also" and "source" match both.
    """
    entries = [["Phred quality score Also written: Q score Source: Ewing & Green 1998",
                "Phred quality score Q score Ewing & Green 1998"],
               ["Q30 Source: Illumina", "Q30 Illumina"]]  # fmt: skip
    queries = ["also", "source", "phred", "illumina q30", ""]
    out = _node(SEARCH, {"entries": entries, "queries": queries})
    assert out == [["0 of 2", True], ["0 of 2", True], ["1 of 2", False], ["1 of 2", False],
                   ["2 of 2", False]]  # fmt: skip


# What gives a fixed box a containing block other than the viewport, which would clip it again.
_TRAPS = re.compile(
    r"(?<![\w-])(transform|filter|backdrop-filter|perspective|contain|will-change|container-type)"
    r"\s*:"
)


def test_the_tooltip_escapes_every_clipping_ancestor():
    """The tooltip is position: fixed, so a scrolling card table cannot clip it; only the top
    bar and the contents rail, which hold no term, set a property that would trap it; and on
    paper a term keeps no room for its hidden icon.

    Mutations: .tip back to position: absolute (the card tables clip it); a transform on
    .table-scroll (a fixed tooltip inside is clipped again); the print rule for .term's
    padding removed (a gap after every first mention on paper).
    """
    html = render(HubPage("Course hub: Demo", "Every page.", (), (), ()))
    css = re.sub(r"/\*.*?\*/", "", re.search(r"<style>(.*?)</style>", html, re.S)[1], flags=re.S)
    rules = [(sel.strip(), body) for sel, body in re.findall(r"([^{}]+)\{([^{}]*)\}", css)]
    assert any(sel == ".tip" and "position:fixed" in body for sel, body in rules)
    assert {sel for sel, body in rules if _TRAPS.search(body)} == {".topbar", ".toc"}
    printed = re.search(r"@media print\{(.*)", css, re.S)[1]
    assert ".term{padding-right:0}" in printed
