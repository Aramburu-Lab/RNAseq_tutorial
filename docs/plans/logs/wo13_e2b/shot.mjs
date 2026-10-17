// Screenshots of page parts in headless Chromium, for an eye check (lessons: check interactive
// figures in the browser, in both themes). Plan on stdin:
// {chromium, shots: [{url, width, scheme, act, selector, out}]}: opens url at width and
// prefers-color-scheme, evaluates act (may return a promise), scrolls selector into view and
// saves a PNG of its box to out.
import { spawn } from "node:child_process";
import { mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";

const plan = JSON.parse(readFileSync(0, "utf8"));
const profile = mkdtempSync(join(tmpdir(), "ctx-shot-"));
const proc = spawn(plan.chromium, ["--headless=new", "--remote-debugging-port=0",
  `--user-data-dir=${profile}`, "--no-first-run", "--disable-gpu", "about:blank"],
  { stdio: ["ignore", "ignore", "pipe"] });
const endpoint = await new Promise((resolve, reject) => {
  let seen = "";
  proc.stderr.on("data", (c) => { seen += c; const f = seen.match(/DevTools listening on (ws:\/\/\S+)/); if (f) resolve(f[1]); });
  proc.on("exit", (code) => reject(new Error(`chromium exited ${code}`)));
});
const ws = new WebSocket(endpoint);
await new Promise((r, j) => { ws.onopen = r; ws.onerror = j; });
let id = 0;
const pending = new Map(), events = [];
ws.onmessage = (e) => {
  const m = JSON.parse(e.data);
  if (m.id && pending.has(m.id)) { const p = pending.get(m.id); pending.delete(m.id); m.error ? p.j(new Error(JSON.stringify(m.error))) : p.r(m.result); }
  else events.push(m);
};
const send = (method, params = {}, sessionId) => new Promise((r, j) => {
  const n = ++id; pending.set(n, { r, j });
  ws.send(JSON.stringify({ id: n, method, params, ...(sessionId ? { sessionId } : {}) }));
});
const done = [];
try {
  for (const shot of plan.shots) {
    const { targetId } = await send("Target.createTarget", { url: "about:blank" });
    const { sessionId } = await send("Target.attachToTarget", { targetId, flatten: true });
    const s = (m, p) => send(m, p, sessionId);
    await s("Page.enable"); await s("Runtime.enable");
    await s("Emulation.setEmulatedMedia", { features: [{ name: "prefers-color-scheme", value: shot.scheme || "light" }] });
    await s("Emulation.setDeviceMetricsOverride", { width: shot.width, height: 900, deviceScaleFactor: 1, mobile: shot.width < 768 });
    await s("Page.navigate", { url: shot.url });
    await new Promise((r) => setTimeout(r, 2500));
    const ev = async (expression) => {
      const out = await s("Runtime.evaluate", { expression, awaitPromise: true, returnByValue: true, userGesture: true });
      if (out.exceptionDetails) throw new Error(JSON.stringify(out.exceptionDetails));
      return out.result.value;
    };
    await ev(shot.act || "null");
    await new Promise((r) => setTimeout(r, 600));
    const box = await ev(`(() => { const e = document.querySelector(${JSON.stringify(shot.selector)});
      e.scrollIntoView(); const b = e.getBoundingClientRect();
      return {x: b.x + window.scrollX, y: b.y + window.scrollY, width: b.width, height: b.height}; })()`);
    const { data } = await s("Page.captureScreenshot", { format: "png", captureBeyondViewport: true,
      clip: { ...box, scale: 1 } });
    writeFileSync(shot.out, Buffer.from(data, "base64"));
    done.push({ out: shot.out, box });
    await send("Target.closeTarget", { targetId });
  }
} finally {
  try { await send("Browser.close"); } catch (e) { /* closing */ }
  proc.kill();
  rmSync(profile, { recursive: true, force: true, maxRetries: 5, retryDelay: 200 });
}
console.log(JSON.stringify(done));
