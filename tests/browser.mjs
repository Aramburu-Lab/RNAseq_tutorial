// Drives headless Chromium over the DevTools protocol with node's own WebSocket, no package:
// tests/test_export_browser.py hands it a plan on stdin and reads what each page did.
// Plan: {chromium, downloads, trace, tasks: [{name, url, width, height, until, act, mouse,
// others, figure, reset, resize}]}. For each task it opens the page, waits, where `until` is
// set, until a loaded page's URL holds it (a page that forwards to another, as an entry page
// does), evaluates `act` as a user gesture, then
// for each step of `mouse` in turn evaluates its `at`, a list of [type, x, y] events, and
// dispatches them as a real left-button mouse (a drag is pressed, moved, released) before
// evaluating its `then`; evaluates `others` in every other page open then (a window the page
// opened), presses the PDF and SVG buttons of the figure `figure` selects when `figure` is
// set, evaluates `reset`, then, for each step of `resize` in turn, sets the window to the
// step's `width` and evaluates its `act`; it reports the results, `moused` and `resized` in
// step order, `others` by URL, the files saved, every request the page made and every console
// error.
import { spawn } from "node:child_process";
import { mkdtempSync, readFileSync, rmSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const plan = JSON.parse(readFileSync(0, "utf8"));
// With plan.trace, each step's start goes to stderr with the seconds since launch, so a
// driver that stops names the step it stopped in.
const began = Date.now();
const trace = plan.trace
  ? (step) => process.stderr.write(`trace ${((Date.now() - began) / 1000).toFixed(1)} ${step}\n`)
  : () => {};
const profile = mkdtempSync(join(tmpdir(), "ctx-browser-"));
const proc = spawn(plan.chromium, [
  "--headless=new", "--remote-debugging-port=0", `--user-data-dir=${profile}`, "--no-first-run",
  "--no-default-browser-check", "--disable-gpu", "--disable-background-networking",
  "--disable-component-update", "--disable-sync", "about:blank",
], { stdio: ["ignore", "ignore", "pipe"] });
const endpoint = await new Promise((resolve, reject) => {
  let seen = "";
  proc.stderr.on("data", (chunk) => {
    seen += chunk;
    const found = seen.match(/DevTools listening on (ws:\/\/\S+)/);
    if (found) resolve(found[1]);
  });
  proc.on("exit", (code) => reject(new Error(`chromium exited ${code}: ${seen}`)));
});
const ws = new WebSocket(endpoint);
await new Promise((resolve, reject) => { ws.onopen = resolve; ws.onerror = reject; });
let lastId = 0;
const pending = new Map(), listeners = new Set();
ws.onmessage = (event) => {
  const msg = JSON.parse(event.data);
  if (msg.id && pending.has(msg.id)) {
    const { resolve, reject } = pending.get(msg.id);
    pending.delete(msg.id);
    if (msg.error) reject(new Error(JSON.stringify(msg.error))); else resolve(msg.result);
  } else if (msg.method) {
    listeners.forEach((listen) => listen(msg));
  }
};
function send(method, params = {}, sessionId = undefined) {
  const id = ++lastId;
  ws.send(JSON.stringify({ id, method, params, ...(sessionId ? { sessionId } : {}) }));
  return new Promise((resolve, reject) => pending.set(id, { resolve, reject }));
}
function waitFor(test, ms) {
  return new Promise((resolve, reject) => {
    const timer = setTimeout(() => { listeners.delete(listen); reject(new Error("timed out")); }, ms);
    function listen(msg) {
      const got = test(msg);
      if (got) { clearTimeout(timer); listeners.delete(listen); resolve(got); }
    }
    listeners.add(listen);
  });
}

// Every download that completed, by guid, so a fast one is not missed between two awaits.
const completed = new Set();
listeners.add((m) => {
  if (m.method === "Browser.downloadProgress" && m.params.state === "completed") {
    completed.add(m.params.guid);
  }
});
async function saved(guid) {
  for (let waited = 0; !completed.has(guid); waited += 100) {
    if (waited > 60000) throw new Error(`download ${guid} did not complete`);
    await new Promise((resolve) => setTimeout(resolve, 100));
  }
}
await send("Browser.setDownloadBehavior",
           // Chromium saves nothing to a path with ".." in it, and says so nowhere.
           { behavior: "allowAndName", downloadPath: resolve(plan.downloads), eventsEnabled: true });
const results = [];
try {
  for (const task of plan.tasks) {
    const { targetId } = await send("Target.createTarget", { url: "about:blank" });
    const { sessionId } = await send("Target.attachToTarget", { targetId, flatten: true });
    const s = (method, params) => send(method, params, sessionId);
    const requests = [], errors = [];
    const watch = (msg) => {
      if (msg.sessionId !== sessionId) return;
      if (msg.method === "Network.requestWillBeSent") requests.push(msg.params.request.url);
      if (msg.method === "Runtime.exceptionThrown") errors.push(msg.params.exceptionDetails.text);
      if (msg.method === "Runtime.consoleAPICalled" && msg.params.type === "error") {
        errors.push(msg.params.args.map((a) => a.value ?? a.description).join(" "));
      }
    };
    listeners.add(watch);
    await s("Network.enable"); await s("Runtime.enable"); await s("Page.enable");
    const size = (width, height) => s("Emulation.setDeviceMetricsOverride", {
      width: width || 1280, height: height || 900, deviceScaleFactor: 1, mobile: (width || 1280) < 768,
    });
    trace(`${task.name}: size`);
    await size(task.width, task.height);
    const loaded = waitFor((m) => m.sessionId === sessionId && m.method === "Page.loadEventFired", 60000);
    trace(`${task.name}: navigate`);
    await s("Page.navigate", { url: task.url });
    await loaded;
    const evaluate = async (expression) => {
      const out = await s("Runtime.evaluate", { expression, awaitPromise: true, returnByValue: true,
                                                userGesture: true });
      if (out.exceptionDetails) throw new Error(`${task.name}: ${JSON.stringify(out.exceptionDetails)}`);
      return out.result.value;
    };
    for (let waited = 0; task.until; waited += 100) {
      // Mid-navigation the old page's context is gone and an evaluation fails: try again.
      const at = await evaluate('document.readyState === "complete" ? location.href : ""')
        .catch(() => "");
      if (at.includes(task.until)) break;
      if (waited > 60000) throw new Error(`${task.name}: no loaded page at ${task.until}`);
      await new Promise((resolve) => setTimeout(resolve, 100));
    }
    trace(`${task.name}: act`);
    const result = { name: task.name, acted: await evaluate(task.act || "null"), files: {} };
    result.moused = [];
    // Input reaches only the page in front: a background tab drew no drag.
    if (task.mouse) await s("Page.bringToFront");
    for (const step of task.mouse || []) {
      for (const [type, x, y] of await evaluate(step.at)) {
        await s("Input.dispatchMouseEvent", { type, x, y, button: "left",
                                              buttons: type === "mouseReleased" ? 0 : 1,
                                              clickCount: type === "mouseMoved" ? 0 : 1 });
      }
      result.moused.push(await evaluate(step.then || "null"));
    }
    if (task.others) {
      result.others = {};
      const { targetInfos } = await send("Target.getTargets");
      for (const info of targetInfos.filter((t) => t.type === "page" && t.targetId !== targetId)) {
        const other = await send("Target.attachToTarget", { targetId: info.targetId, flatten: true });
        const out = await send("Runtime.evaluate", { expression: task.others, awaitPromise: true,
          returnByValue: true, userGesture: true }, other.sessionId);
        result.others[info.url] = out.exceptionDetails ? { error: out.exceptionDetails.text }
          : out.result.value;
      }
    }
    for (const ext of task.figure ? ["pdf", "svg"] : []) {
      const begun = waitFor((m) => m.method === "Browser.downloadWillBegin" && m.params, 30000);
      await evaluate(`document.querySelector(${JSON.stringify(task.figure)})` +
                     `.querySelector('a[data-export="${ext}"]').click()`);
      const { guid, suggestedFilename } = await begun;
      await saved(guid);
      result.files[ext] = { guid, name: suggestedFilename };
    }
    result.reset = await evaluate(task.reset || "null");
    result.resized = [];
    for (const step of task.resize || []) {
      trace(`${task.name}: resize ${step.width}`);
      await size(step.width, step.height);
      result.resized.push(await evaluate(step.act));
    }
    result.requests = requests;
    result.errors = errors;
    results.push(result);
    listeners.delete(watch);
    trace(`${task.name}: close`);
    await send("Target.closeTarget", { targetId });
  }
} finally {
  const exited = new Promise((resolve) => proc.on("exit", resolve));
  trace("browser close");
  try { await send("Browser.close"); } catch (e) { /* already closing */ }
  trace("browser exit");
  await Promise.race([exited, new Promise((resolve) => setTimeout(resolve, 5000))]);
  proc.kill();
  trace("profile removal");
  rmSync(profile, { recursive: true, force: true, maxRetries: 5, retryDelay: 200 });
  trace("done");
}
console.log(JSON.stringify(results));
