// Drives headless Chromium over the DevTools protocol with node's own WebSocket, no package:
// tests/test_export_browser.py hands it a plan on stdin and reads what each page did.
// Plan: {chromium, downloads, trace, throttle, tasks: [{name, url, width, height, scheme,
// throttle, until, act, mouse, others, figure, reset, resize}]}. For each task it opens the page, waits, where `until` is
// set, until a loaded page's URL holds it (a page that forwards to another, as an entry page
// does), evaluates `act` as a user gesture, then
// for each step of `mouse` in turn evaluates its `at`, a list of [type, x, y] events, and
// dispatches them as a real left-button mouse (a drag is pressed, moved, released), waits,
// where the step's `until` is set, for a loaded page whose URL holds it (a click on a link),
// then evaluates its `then`; evaluates `others` in every other page open then (a window the page
// opened), presses the PDF and SVG buttons of the figure `figure` selects when `figure` is
// set, evaluates `reset`, then, for each step of `resize` in turn, sets the window to the
// step's `width` and evaluates its `act`; it reports the results, `moused` and `resized` in
// step order, `others` by URL, the files saved, every request the page made and every console
// error. Sent SIGTERM (the caller's time limit), it pauses the page in front of it through the
// debugger, prints where its script is running to stderr and exits 3 (H.6).
import { spawn } from "node:child_process";
import { mkdtempSync, readFileSync, rmSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const plan = JSON.parse(readFileSync(0, "utf8"));
// With plan.trace, each step's start goes to stderr with the seconds since launch, so a
// driver that stops names the step it stopped in.
const began = Date.now();
// The page and step the driver is in, for the diagnosis a SIGTERM asks for.
let current = { sessionId: null, step: "launch" };
const trace = (step) => {
  current.step = step;
  if (plan.trace) process.stderr.write(`trace ${((Date.now() - began) / 1000).toFixed(1)} ${step}\n`);
};
const profile = mkdtempSync(join(tmpdir(), "ctx-browser-"));
// plan.headed shows the window, as a student's double-click does; the tests run headless.
const proc = spawn(plan.chromium, [
  ...(plan.headed ? [] : ["--headless=new"]), "--remote-debugging-port=0", `--user-data-dir=${profile}`, "--no-first-run",
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
// A browser gone (closed, crashed or killed) answers no command: each one waiting fails then,
// and each later one at once, instead of waiting for the caller's time limit (H.6).
let gone = null;
ws.onclose = () => {
  gone = new Error(`the browser closed its DevTools connection at ${((Date.now() - began) / 1000).toFixed(1)} s, in ${current.step}`);
  pending.forEach(({ reject }) => reject(gone));
  pending.clear();
};
function send(method, params = {}, sessionId = undefined) {
  if (gone) return Promise.reject(gone);
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

// A page that never yields to a timer answers no Runtime.evaluate, and no Debugger.enable
// either, but V8 runs Debugger.pause in an interrupt; so each page is opened with its debugger
// enabled, and a hang is paused a few times to name the frames it is in, with script and line.
let diagnosing = false;
const scripts = new Map(), crashed = new Map();
listeners.add((m) => {
  if (m.method === "Debugger.scriptParsed") scripts.set(m.params.scriptId, m.params.url);
  // A renderer killed under memory pressure answers nothing either, and is no loop.
  if (m.method === "Inspector.targetCrashed") crashed.set(m.sessionId, (Date.now() - began) / 1000);
  // A debugger statement in a page script must not stop a test the driver is not diagnosing.
  if (m.method === "Debugger.paused" && !diagnosing) send("Debugger.resume", {}, m.sessionId).catch(() => null);
});
async function diagnose() {
  const { sessionId, step } = current;
  if (!sessionId) return `hang diagnosis: in ${step}, no page open\n`;
  if (crashed.has(sessionId)) {
    return `hang diagnosis: in ${step}, the page's renderer crashed at ${crashed.get(sessionId).toFixed(1)} s\n`;
  }
  diagnosing = true;
  const within = (promise, ms) => Promise.race([promise, new Promise((resolve) => setTimeout(resolve, ms, null))]);
  const lines = [`hang diagnosis: in ${step}`];
  for (let sample = 1; sample <= 3; sample++) {
    const paused = waitFor((m) => m.sessionId === sessionId && m.method === "Debugger.paused" && m.params, 5000)
      .catch(() => null);
    await within(send("Debugger.pause", {}, sessionId), 5000);
    const got = await paused;
    if (!got) {
      lines.push(`sample ${sample}: no script ran within 5 s; the page is idle or blocked outside JavaScript`);
      break;
    }
    lines.push(`sample ${sample}:`);
    for (const f of got.callFrames.slice(0, 8)) {
      const url = f.url || scripts.get(f.location.scriptId) || "?";
      lines.push(`  at ${f.functionName || "(anonymous)"} ${url}:${f.location.lineNumber + 1}:${f.location.columnNumber + 1}`);
    }
    await within(send("Debugger.resume", {}, sessionId), 5000);
    await new Promise((resolve) => setTimeout(resolve, 200));
  }
  return lines.join("\n") + "\n";
}
process.on("SIGTERM", async () => {
  process.stderr.write(await diagnose().catch((e) => `hang diagnosis failed: ${e.message}\n`));
  proc.kill("SIGKILL");
  rmSync(profile, { recursive: true, force: true, maxRetries: 5, retryDelay: 200 });
  process.exit(3);
});

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
    current.sessionId = sessionId;
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
    await s("Debugger.enable"); await s("Inspector.enable");
    // Pages follow prefers-color-scheme when no theme is stored, and macOS switches it by the
    // clock: the light tasks drew dark after 20:00 (2026-10-17). A task asks for dark itself.
    await s("Emulation.setEmulatedMedia", {
      features: [{ name: "prefers-color-scheme", value: task.scheme || "light" }],
    });
    const size = (width, height) => s("Emulation.setDeviceMetricsOverride", {
      width: width || 1280, height: height || 900, deviceScaleFactor: 1, mobile: (width || 1280) < 768,
    });
    // task.throttle slows this page's CPU by that factor (or plan.throttle every page's): a
    // deterministic stand-in for a loaded machine, where timing races show (H.6).
    const rate = task.throttle || plan.throttle;
    if (rate) await s("Emulation.setCPUThrottlingRate", { rate });
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
    const arrived = async (until) => {
      for (let waited = 0; until; waited += 100) {
        // Mid-navigation the old page's context is gone and an evaluation fails: try again.
        const at = await evaluate('document.readyState === "complete" ? location.href : ""')
          .catch(() => "");
        if (at.includes(until)) return;
        if (waited > 60000) throw new Error(`${task.name}: no loaded page at ${until}`);
        await new Promise((resolve) => setTimeout(resolve, 100));
      }
    };
    await arrived(task.until);
    trace(`${task.name}: act`);
    const result = { name: task.name, acted: await evaluate(task.act || "null"), files: {} };
    result.moused = [];
    // Input reaches only the page in front: a background tab drew no drag.
    if (task.mouse) await s("Page.bringToFront");
    for (const [i, step] of (task.mouse || []).entries()) {
      trace(`${task.name}: mouse ${i}`);
      for (const [type, x, y] of await evaluate(step.at)) {
        const sent = s("Input.dispatchMouseEvent", { type, x, y, button: "left",
                                                     buttons: type === "mouseReleased" ? 0 : 1,
                                                     clickCount: type === "mouseMoved" ? 0 : 1 });
        // A click that follows a link navigates before the release is acknowledged.
        await (step.until ? sent.catch(() => null) : sent);
      }
      await arrived(step.until);
      result.moused.push(await evaluate(step.then || "null"));
    }
    if (task.others) {
      trace(`${task.name}: others`);
      result.others = {};
      const { targetInfos } = await send("Target.getTargets");
      for (const info of targetInfos.filter((t) => t.type === "page" && t.targetId !== targetId)) {
        // A shown browser opens and closes pages of its own (a sign-in page) meanwhile.
        const other = await send("Target.attachToTarget", { targetId: info.targetId, flatten: true })
          .catch(() => null);
        if (!other) continue;
        const out = await send("Runtime.evaluate", { expression: task.others, awaitPromise: true,
          returnByValue: true, userGesture: true }, other.sessionId);
        result.others[info.url] = out.exceptionDetails ? { error: out.exceptionDetails.text }
          : out.result.value;
      }
    }
    for (const ext of task.figure ? ["pdf", "svg"] : []) {
      trace(`${task.name}: figure ${ext}`);
      const begun = waitFor((m) => m.method === "Browser.downloadWillBegin" && m.params, 30000);
      await evaluate(`document.querySelector(${JSON.stringify(task.figure)})` +
                     `.querySelector('a[data-export="${ext}"]').click()`);
      const { guid, suggestedFilename } = await begun;
      await saved(guid);
      result.files[ext] = { guid, name: suggestedFilename };
    }
    trace(`${task.name}: reset`);
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
    current.sessionId = null;
    await send("Target.closeTarget", { targetId });
  }
} finally {
  const exited = new Promise((resolve) => proc.on("exit", resolve));
  trace("browser close");
  try { await send("Browser.close"); } catch (e) { /* already closing */ }
  trace("browser exit");
  await Promise.race([exited, new Promise((resolve) => setTimeout(resolve, 5000))]);
  proc.kill();
  /* Chromium's stderr is one socket shared with every helper it started, and with the two
     crashpad handlers it detaches from itself: while any of them lives, reading it kept this
     driver's event loop, and so the driver, alive after "done". On a loaded Mac they outlived
     the browser by minutes, and drivers stood idle past their 900 s limit (H.6, 2026-10-21). */
  proc.stderr.destroy();
  trace("profile removal");
  rmSync(profile, { recursive: true, force: true, maxRetries: 5, retryDelay: 200 });
  trace("done");
}
console.log(JSON.stringify(results));
