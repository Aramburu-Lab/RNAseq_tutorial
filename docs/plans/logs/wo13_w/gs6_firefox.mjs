// Package W, P.1: GS.6 in Firefox, which speaks WebDriver BiDi and no longer CDP, so
// tests/browser.mjs cannot drive it. Headless Firefox opens the M5 page, selects the up calls,
// presses Open beside on the M9 link, reads the M9 window it opened, then selects the down
// calls on M5 and reads the M9 window again: the list by link, then by message.
// Run: node gs6_firefox.mjs <firefox> <M5 page file URL>; prints JSON.
import { spawn } from "node:child_process";
import { mkdtempSync, realpathSync, rmSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";

const [firefox, m5] = process.argv.slice(2);
const profile = realpathSync(mkdtempSync(join(process.env.GS6_PROFILE_DIR || tmpdir(), "w-gs6-firefox-")));
const proc = spawn(firefox, ["--headless", "--remote-debugging-port", "0", "-profile", profile,
                             "about:blank"], { stdio: ["ignore", "pipe", "pipe"] });
const endpoint = await new Promise((resolve, reject) => {
  let seen = "";
  const look = (chunk) => {
    seen += chunk;
    const found = seen.match(/WebDriver BiDi listening on (ws:\/\/\S+)/);
    if (found) resolve(found[1]);
  };
  proc.stdout.on("data", look); proc.stderr.on("data", look);
  proc.on("exit", (code) => reject(new Error(`firefox exited ${code}: ${seen}`)));
  setTimeout(() => reject(new Error(`no BiDi endpoint in 60 s: ${seen}`)), 60000);
});
const ws = new WebSocket(`${endpoint}/session`);
await new Promise((resolve, reject) => { ws.onopen = resolve; ws.onerror = reject; });
let lastId = 0;
const pending = new Map();
ws.onmessage = (event) => {
  const msg = JSON.parse(event.data);
  if (msg.id && pending.has(msg.id)) {
    const { resolve, reject } = pending.get(msg.id);
    pending.delete(msg.id);
    if (msg.type === "error") reject(new Error(`${msg.error}: ${msg.message}`)); else resolve(msg.result);
  }
};
function send(method, params = {}) {
  const id = ++lastId;
  ws.send(JSON.stringify({ id, method, params }));
  return new Promise((resolve, reject) => pending.set(id, { resolve, reject }));
}
const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));
async function evaluate(context, expression) {
  const out = await send("script.evaluate", { expression, target: { context }, awaitPromise: true,
                                              userActivation: true, resultOwnership: "none" });
  if (out.type === "exception") throw new Error(JSON.stringify(out.exceptionDetails));
  return JSON.parse(out.result.value);
}
const READ = `(async () => {
  let panel = null;
  for (let i = 0; i < 300 && !(panel = document.querySelector(".genelist-panel")); i++) await new Promise(r => setTimeout(r, 100));
  for (let i = 0; i < 300 && panel && !panel.ctxList; i++) await new Promise(r => setTimeout(r, 100));
  return JSON.stringify({url: location.href.slice(0, 160), opener: !!window.opener,
    list: panel && panel.ctxList ? {found: panel.ctxList.found.length, missing: panel.ctxList.missing, cards: panel.ctxList.cards} : null,
    summary: panel ? panel.querySelector(".list-summary").textContent : null});
})()`;
const out = { firefox: null };
try {
  await send("session.new", { capabilities: {} }).then((r) => { out.firefox = r.capabilities.browserVersion; });
  const tree = await send("browsingContext.getTree", {});
  const main = tree.contexts[0].context;
  await send("browsingContext.navigate", { context: main, url: m5, wait: "complete" });
  await sleep(3000);
  out.picked = await evaluate(main, `(async () => {
    document.querySelector("button[data-pick=up]").click();
    await new Promise(r => setTimeout(r, 500));
    const a = Array.from(document.querySelectorAll("a[data-onward]")).filter(x => x.getAttribute("href").indexOf("11_genecard") >= 0)[0];
    a.parentElement.querySelector("[data-beside]").click();
    return JSON.stringify({link: a.textContent});
  })()`);
  let beside = null;
  for (let i = 0; i < 100 && !beside; i++) {
    await sleep(200);
    const now = await send("browsingContext.getTree", {});
    beside = now.contexts.map((c) => c.context).filter((c) => c !== main)[0] || null;
  }
  out.opened = !!beside;
  if (beside) {
    out.byLink = await evaluate(beside, READ);
    await evaluate(main, `(async () => { document.querySelector("button[data-pick=down]").click(); return "1"; })()`);
    await sleep(3000);
    out.byMessage = await evaluate(beside, READ);
  }
} catch (err) {
  out.error = String(err);
} finally {
  try { await send("browser.close"); } catch (e) { /* closing */ }
  await sleep(2000);
  proc.kill();
  rmSync(profile, { recursive: true, force: true, maxRetries: 5, retryDelay: 200 });
}
console.log(JSON.stringify(out, null, 1));
