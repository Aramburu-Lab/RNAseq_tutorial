"""P.1: a burst of keystrokes in the live-calls form redraws once (``ctxCalls.debounce``).

One redraw of a 24,375-point volcano took 897 to 933 ms headless (D-081); a redraw per key
held a burst for seconds. The rule runs under node on a test clock; each test names the
mutation it catches, each planted and seen to fail.
"""

from test_page_script import _node

CLOCK = r"""
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
/* A clock the test moves: set() and clear() as setTimeout and clearTimeout, tick() to a time. */
let now = 0, next = 1;
const due = new Map();
const timers = {
  set: (f, d) => { const h = next++; due.set(h, [now + d, f]); return h; },
  clear: (h) => { due.delete(h); },
};
function tick(to) {
  for (;;) {
    const ready = [...due].filter(([, [at]]) => at <= to).sort((a, b) => a[1][0] - b[1][0]);
    if (!ready.length) break;
    const [h, [at, f]] = ready[0];
    due.delete(h); now = at; f();
  }
  now = to;
}
const out = {};
for (const [name, steps] of Object.entries(input)) {
  const runs = [];
  now = 0; due.clear();
  const later = ctxCalls.debounce(() => runs.push(now), ctxCalls.QUIET_MS, timers);
  for (const [kind, at] of steps) {
    tick(at);
    if (kind === "key") later(); else if (kind === "cancel") later.cancel();
  }
  tick(1e9);
  out[name] = runs;
}
out.quiet_ms = ctxCalls.QUIET_MS;
console.log(JSON.stringify(out));
"""


def test_a_burst_of_keys_redraws_once_after_the_pause():
    """Five keys 100 ms apart redraw once, QUIET_MS after the last; two bursts a second apart
    redraw twice. Mutations caught: no debounce (a redraw per key); the pending run not
    cleared on a new key (a redraw per key, each late); the delay taken from the first key."""
    burst = [["key", 100 * k] for k in range(5)]
    later = [["key", 2000 + 100 * k] for k in range(3)]
    got = _node({"burst": burst, "two": burst + later}, CLOCK)
    quiet = got["quiet_ms"]
    assert quiet == 300  # noqa: PLR2004 - P.1's pause, stated apart from the code
    assert got["burst"] == [400 + quiet]
    assert got["two"] == [400 + quiet, 2200 + quiet]


def test_a_change_drops_the_redraw_typing_left_waiting():
    """A cancel, as an immediate redraw makes, drops the pending one; keys after it start
    anew. Mutations caught: cancel() that leaves the timer set; cancel() that runs it."""
    steps = [["key", 0], ["key", 100], ["cancel", 150], ["key", 1000]]
    got = _node({"cancelled": steps, "only": [["key", 0], ["cancel", 10]]}, CLOCK)
    assert got["cancelled"] == [1000 + got["quiet_ms"]]
    assert got["only"] == []
